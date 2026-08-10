//! Empirical exact-set tripwire over the output-collision pre-flight
//! ([#414](https://github.com/FelixKrueger/TrimGalore/issues/414)).
//!
//! The same defect has shipped five times — #383, #388, #391, #409, and the report
//! hole #398 closed. Every instance was a writer whose path never entered the
//! pre-flight's candidate list, so the pre-flight silently stopped covering that
//! output; #409 destroyed live input data. `OutputSource` typing (#397) cannot catch
//! it, because a type describes candidates that exist.
//!
//! Each case here snapshots a tree, runs one dispatch arm, diffs, and asserts every
//! created file was in the candidate list **the binary itself dumped** — never in a
//! list written by hand here. Hand-enumerating would relocate the original knowledge
//! problem into this file, where it would rot the same way.
//!
//! ## If a case goes red
//!
//! The fix is to add the missing path to that arm's candidate list in `src/main.rs`,
//! at the site the failure names. Do **not** add it to the FastQC allowance and do
//! **not** narrow a candidate list — those are how this test dies. #398 hit the same
//! shape from the other side, where narrowing was the tempting wrong fix for a red
//! `..._do_not_over_reject`.
//!
//! ## Known-unplanned writes
//!
//! `--fastqc` output (`<primary-stem>_fastqc.{html,zip}`) is deliberately absent from
//! every candidate list, so it is allowed here by rule — but only on arms that pass a
//! FastQC flag, and on those arms the allowed set must be non-empty. `--fastqc_args
//! "-o DIR"` relocates both artefacts into `DIR`, which must already exist.
//!
//! Those two suffixes are the whole set: `--extract` is not in the curated `--fastqc_args`
//! subset, so no `_fastqc/` directory appears, and `--svg` sets a rendering flag without
//! writing an `.svg` file.
//!
//! `--fastqc` is a silent no-op on `--hardtrim5/3`, `--clock`, `--implicon` and on
//! paired FASTQ output from a single interleaved uBAM ([#421]). Those arms are marked
//! `fastqc_incapable`; do not pair them with a FastQC flag expecting artefacts.
//!
//! [#421]: https://github.com/FelixKrueger/TrimGalore/issues/421
//!
//! ## Limits
//!
//! Comparison uses the production `io::collision_key`, which folds case. A candidate
//! that differs from the written path only in case is therefore invisible here, which
//! on a case-sensitive filesystem is a real defect (writer writes `X`, pre-flight
//! guards `x`). Using the pre-flight's own relation is still right: a stricter one
//! would disagree with the code under audit about APFS aliasing (#216).
//!
//! A write that escapes the case root is invisible: `snapshot` only walks under it. The
//! planned side is guarded in `compare`, the created side cannot be. No arm can escape
//! while every fixture path stays relative and the child runs with the root as its CWD.
//!
//! A new dispatch arm with **no** pre-flight is not detected. `the_site_inventory_is_complete`
//! catches a 12th call site; an arm that never calls the pre-flight looks like the two
//! `Arm::Unplanned` cases, which are on record only because someone listed them.

use std::collections::{BTreeMap, BTreeSet};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use trim_galore::io::{DUMP_PLANNED_ENV, DUMP_PLANNED_PREFIX, collision_key};

// ── harness ──────────────────────────────────────────────────────────────────

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_trim_galore"))
}

/// Committed fixtures resolve from the manifest dir, not the CWD: every case runs the
/// binary with `current_dir` set to its own root, under which `test_files/…` would not
/// resolve.
fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("test_files")
        .join(name)
}

/// Canonicalised, because on macOS `temp_dir()` yields `/tmp/…` while the child's
/// `getcwd` yields `/private/tmp/…`; comparing the two spellings would pass vacuously.
///
/// `tag` must be unique across this file. All cases share one PID under
/// `cargo test --release`, and this opens with `remove_dir_all`, so a duplicate tag
/// deletes a sibling's tree mid-run — and a deletion that removes an unplanned file
/// would turn a real miss into a pass.
fn case_root(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tg_trip_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    std::fs::canonicalize(&d).unwrap_or(d)
}

/// 40 reads of 80 bp. Mates must share `id_prefix`: `--passthrough` aborts on
/// out-of-sync read IDs, and a case that cannot complete its run tests nothing.
fn fastq_body(id_prefix: &str) -> String {
    let seq = "ACGT".repeat(20);
    let qual = "I".repeat(80);
    (0..40)
        .map(|i| format!("@{id_prefix}read{i}\n{seq}\n+\n{qual}\n"))
        .collect()
}

fn write_fastq(path: &Path, id_prefix: &str) {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    std::fs::write(path, fastq_body(id_prefix)).unwrap();
}

/// Real gzip bytes, not a plain file named `.gz`: `is_gzipped` is extension-based, so
/// a fake would exercise the naming branch but never the gzip writer.
fn write_fastq_gz(path: &Path, id_prefix: &str) {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    let f = std::fs::File::create(path).unwrap();
    let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::fast());
    enc.write_all(fastq_body(id_prefix).as_bytes()).unwrap();
    enc.finish().unwrap();
}

/// Short R2 so reads genuinely strand under `--retain_unpaired`.
fn write_fastq_short(path: &Path, id_prefix: &str) {
    let body: String = (0..40)
        .map(|i| format!("@{id_prefix}read{i}\nACGTA\n+\nIIIII\n"))
        .collect();
    std::fs::write(path, body).unwrap();
}

fn copy_fixture(name: &str, dest: &Path) {
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    std::fs::copy(fixture(name), dest).unwrap();
}

fn write_barcodes(path: &Path) {
    // Real tabs. The first barcode matches every read's 3' end; the second matches
    // none, and its output file is planned regardless.
    std::fs::write(path, "sampleA\tACGTACGT\nsampleB\tTTTTGGGG\n").unwrap();
}

/// Files only, relative to `root`. `ensure_output_dir` creates directories before
/// dispatch; a directory is not a write.
fn snapshot(root: &Path) -> BTreeSet<PathBuf> {
    let mut out = BTreeSet::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let entries =
            std::fs::read_dir(&dir).unwrap_or_else(|e| panic!("read_dir {}: {e}", dir.display()));
        for entry in entries {
            let entry = entry.unwrap();
            let path = entry.path();
            if entry.file_type().unwrap().is_dir() {
                stack.push(path);
            } else {
                out.insert(path.strip_prefix(root).unwrap().to_path_buf());
            }
        }
    }
    out
}

fn hash_file(path: &Path) -> u64 {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut h = DefaultHasher::new();
    bytes.hash(&mut h);
    h.finish()
}

fn run_arm(root: &Path, args: &[&str]) -> (bool, String) {
    let out = Command::new(binary())
        .current_dir(root)
        .env(DUMP_PLANNED_ENV, "1")
        .args(args)
        .output()
        .expect("failed to run trim_galore");
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stderr).to_string(),
    )
}

/// `(site, path)` per dumped candidate. The site is repeated on every line, so this
/// needs no state; `compare` then requires all of them to agree.
fn parse_planned(stderr: &str) -> Vec<(String, PathBuf)> {
    stderr
        .lines()
        .filter_map(|l| l.strip_prefix(DUMP_PLANNED_PREFIX))
        .map(|rest| {
            let (site, path) = rest
                .split_once('\t')
                .unwrap_or_else(|| panic!("malformed dump line: {rest:?}"));
            (site.to_string(), PathBuf::from(path))
        })
        .collect()
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Arm {
    /// Reaches a pre-flight; the full assertion set applies.
    Planned,
    /// No pre-flight exists on this dispatch path (#414 §2.3). T1 and T5 cannot
    /// apply; `planned` must be empty, so the day one is added this goes red.
    ///
    /// Safe today only because both such arms derive their stem with
    /// `Path::file_stem()` and never read `cli.basename`: `file_stem` strips a
    /// `.`-prefixed extension while every appended suffix starts with `_`, so an output
    /// can never equal its input. Honouring `--basename` there would break that.
    Unplanned,
}

fn is_fastqc_artifact(path: &Path) -> bool {
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    name.ends_with("_fastqc.html") || name.ends_with("_fastqc.zip")
}

/// Production's gate is `cli.fastqc || cli.fastqc_args.is_some()` (`main.rs:1458` and
/// four siblings). The disjunction is load-bearing: `--fastqc_args`' "Implies --fastqc"
/// is doc text with no clap `requires`, so `cli.fastqc` is false on an args-only run
/// while FastQC still runs. Deriving the predicate here, once, from the args keeps one
/// source of truth for it.
fn fastqc_active(args: &[&str]) -> bool {
    args.iter()
        .any(|a| *a == "--fastqc" || *a == "--fastqc_args")
}

struct Tripwire<'a> {
    tag: &'a str,
    args: &'a [&'a str],
    /// Relative paths this case must have staged, asserted before the run. T3 hashes
    /// every pre-existing file regardless, so this is a staging check rather than the
    /// definition of the guarded set. List `guarded_inputs` (`main.rs:74-83`) — the
    /// positional inputs plus the `--passthrough` carrier plus the `--demux` barcode
    /// file, the last two being what #389 and #409 were about.
    guarded: &'a [&'a str],
    arm: Arm,
    fastqc_capable: bool,
}

/// What a case observed, for the few tests that assert something extra.
struct Outcome {
    stderr: String,
    created: BTreeSet<PathBuf>,
}

impl<'a> Tripwire<'a> {
    fn new(tag: &'a str, args: &'a [&'a str], guarded: &'a [&'a str]) -> Self {
        Self {
            tag,
            args,
            guarded,
            arm: Arm::Planned,
            fastqc_capable: true,
        }
    }

    fn unplanned(mut self) -> Self {
        self.arm = Arm::Unplanned;
        self
    }

    fn fastqc_incapable(mut self) -> Self {
        self.fastqc_capable = false;
        self
    }

    fn run(self, stage: impl FnOnce(&Path)) -> Outcome {
        let root = case_root(self.tag);
        stage(&root);

        let before = snapshot(&root);
        for rel in self.guarded {
            assert!(
                before.contains(Path::new(rel)),
                "{}: guarded input {rel} was not staged",
                self.tag
            );
        }
        // T3 covers every pre-existing file, so a case that stages one without declaring
        // it in `guarded` is still guarded.
        let hashes: BTreeMap<PathBuf, u64> = before
            .iter()
            .map(|rel| {
                let p = root.join(rel);
                (p.clone(), hash_file(&p))
            })
            .collect();

        let (ok, stderr) = run_arm(&root, self.args);
        assert!(
            ok,
            "{}: expected exit 0, args {:?}\nstderr:\n{stderr}",
            self.tag, self.args
        );

        let after = snapshot(&root);
        let created: BTreeSet<PathBuf> = after.difference(&before).cloned().collect();
        let planned = parse_planned(&stderr);

        if let Err(msg) = unchanged(&hashes) {
            panic!("{}: {msg}\nstderr:\n{stderr}", self.tag);
        }

        if let Err(msg) = compare(
            &root,
            &created,
            &planned,
            self.args,
            self.arm,
            self.fastqc_capable,
        ) {
            panic!("{}: {msg}\nstderr:\n{stderr}", self.tag);
        }

        Outcome { stderr, created }
    }
}

/// T3 — no pre-existing file's bytes changed, which is #409's harm asserted directly and
/// the substantive assertion on the `Arm::Unplanned` cases.
///
/// Pure and separate from `compare` so `comparison_detects_a_modified_input` can prove it
/// fires. Without that, the only proof was a planted mutation in `main.rs` that gets
/// reverted, and every other assertion here has a standing one.
///
/// `hash_file` panics on a missing path, so a *deleted* pre-existing file fails too.
fn unchanged(before: &BTreeMap<PathBuf, u64>) -> Result<(), String> {
    for (path, before_hash) in before {
        if hash_file(path) != *before_hash {
            return Err(format!(
                "pre-existing file {} was modified by the run",
                path.display()
            ));
        }
    }
    Ok(())
}

/// The observed-vs-planned comparison. Pure, so `comparison_*` can drive it directly
/// and prove it goes red without running the binary.
fn compare(
    root: &Path,
    created: &BTreeSet<PathBuf>,
    planned: &[(String, PathBuf)],
    args: &[&str],
    arm: Arm,
    fastqc_capable: bool,
) -> Result<(), String> {
    // T4 — liveness. `main.rs:215-230` routes clap's
    // DisplayHelpOnMissingArgumentOrSubcommand to stdout with exit 0, so a case whose
    // positional list is swallowed by a flag would satisfy every subset relation below
    // vacuously, exit-0 gate included.
    if created.is_empty() {
        return Err("the run created no files at all — it never reached its arm".into());
    }

    let active = fastqc_active(args);
    let (allowed, rest): (Vec<&PathBuf>, Vec<&PathBuf>) = created
        .iter()
        .partition(|p| active && is_fastqc_artifact(p));

    // The allowance is self-policing, so it cannot quietly widen into the drain this
    // test dies in. Scoped to arms that can reach `fastqc::run` at all: unscoped, it
    // would go red on a legitimate specialty case (#421), and a provably-wrong red is
    // the strongest argument for deleting a guard.
    match (active, fastqc_capable) {
        (true, true) if allowed.is_empty() => {
            return Err("a FastQC flag is set but no FastQC artefact was written; \
                        the allowance is masking a broken case"
                .into());
        }
        (true, false) if !allowed.is_empty() => {
            return Err(format!(
                "arm is marked fastqc_incapable but produced {allowed:?}; \
                 #421 may have been fixed — drop the marker"
            ));
        }
        _ => {}
    }
    // An allowance that swallowed the whole write set would leave T1 comparing empty sets.
    if rest.is_empty() {
        return Err(format!(
            "the FastQC allowance took every created file ({allowed:?}), so T1 has nothing \
             left to check — narrow `is_fastqc_artifact`, do not widen it"
        ));
    }
    if !active {
        let stray: Vec<&&PathBuf> = rest.iter().filter(|p| is_fastqc_artifact(p)).collect();
        if !stray.is_empty() {
            return Err(format!(
                "FastQC artefacts {stray:?} appeared with no FastQC flag in args"
            ));
        }
    }

    if arm == Arm::Unplanned {
        if !planned.is_empty() {
            return Err(format!(
                "this arm had no pre-flight, but one dumped {} candidates from {} — \
                 promote the case to Arm::Planned",
                planned.len(),
                planned[0].0
            ));
        }
        return Ok(());
    }

    if planned.is_empty() {
        return Err(
            "no candidates were dumped — the hook did not fire, or the run \
                    never reached a pre-flight"
                .into(),
        );
    }

    // T5 — no run reaches two DISTINCT pre-flight sites, so T1 compares against one arm's
    // list. It does not detect a single site invoked repeatedly: the dump carries no
    // invocation boundary, so one call of 12 candidates and two of 6 are byte-identical.
    // That shape is harmless anyway — the union of per-pair lists equals the single list.
    let site = &planned[0].0;
    if let Some((other, path)) = planned.iter().find(|(s, _)| s != site) {
        return Err(format!(
            "two pre-flights ran in one invocation ({site} and {other}, e.g. {}); \
             T1 would be comparing against a union of two arms' lists",
            path.display()
        ));
    }

    let root_key = collision_key(root);
    let planned_keys: BTreeSet<String> = planned
        .iter()
        .map(|(_, p)| collision_key(&root.join(p)))
        .collect();

    // A candidate resolving outside the case root would be invisible to `snapshot`,
    // making T1 vacuous for it.
    let root_prefix = format!("{root_key}{}", std::path::MAIN_SEPARATOR);
    for (_, p) in planned {
        let key = collision_key(&root.join(p));
        if !key.starts_with(&root_prefix) {
            return Err(format!(
                "planned path {} resolves outside the case root, where the snapshot \
                 cannot see it",
                p.display()
            ));
        }
    }

    // T1 — planned coverage. Only this direction: over-planning is legitimate
    // (`--retain_unpaired` plans `_unpaired_N` whether or not anything strands, and
    // `--demux` plans a file per barcode), so `planned ⊆ created` is false by design.
    let unplanned: Vec<String> = rest
        .iter()
        .filter(|p| !planned_keys.contains(&collision_key(&root.join(p))))
        .map(|p| p.display().to_string())
        .collect();

    if !unplanned.is_empty() {
        let mut list: Vec<String> = planned
            .iter()
            .map(|(_, p)| p.display().to_string())
            .collect();
        list.sort();
        return Err(format!(
            "created but never planned: {unplanned:?}\n\
             Add each to the candidate list at {site}. Do NOT add them to the FastQC \
             allowance and do NOT narrow the candidate list — those are how this \
             tripwire dies.\n\
             The {} candidates that arm did plan: {list:?}",
            planned.len()
        ));
    }

    Ok(())
}

// ── staging helpers ──────────────────────────────────────────────────────────

fn stage_se(root: &Path) {
    write_fastq(&root.join("a.fastq"), "A_");
}

fn stage_se2(root: &Path) {
    write_fastq(&root.join("a.fastq"), "A_");
    write_fastq(&root.join("b.fastq"), "B_");
}

fn stage_pair(root: &Path) {
    write_fastq(&root.join("s_R1.fastq"), "S_");
    write_fastq(&root.join("s_R2.fastq"), "S_");
}

/// R1 and R2 in *different* directories, so a report written beside its own mate
/// instead of beside R1 fails T1. This is #398's shape.
fn stage_pair_split(root: &Path) {
    write_fastq(&root.join("r1/s_R1.fastq"), "S_");
    write_fastq(&root.join("r2/s_R2.fastq"), "S_");
}

// ── site 1006 — single-end FASTQ trim ────────────────────────────────────────

#[test]
fn se_trim_two_inputs() {
    Tripwire::new(
        "se_trim_two_inputs",
        &["a.fastq", "b.fastq"],
        &["a.fastq", "b.fastq"],
    )
    .run(stage_se2);
}

#[test]
fn se_trim_output_dir() {
    Tripwire::new(
        "se_trim_output_dir",
        &["-o", "out", "in/a.fastq", "in/b.fastq"],
        &["in/a.fastq", "in/b.fastq"],
    )
    .run(|root| {
        write_fastq(&root.join("in/a.fastq"), "A_");
        write_fastq(&root.join("in/b.fastq"), "B_");
    });
}

#[test]
fn se_trim_fastqc() {
    Tripwire::new("se_trim_fastqc", &["--fastqc", "a.fastq"], &["a.fastq"]).run(stage_se);
}

#[test]
fn se_trim_fastqc_args_outdir() {
    // `qc` is pre-created on purpose: fastqc-rust does not create its `-o` directory,
    // and without it the run trims, writes both reports, then fails at FastQC.
    Tripwire::new(
        "se_trim_fastqc_args_outdir",
        &["--fastqc", "--fastqc_args", "-o qc", "a.fastq"],
        &["a.fastq"],
    )
    .run(|root| {
        stage_se(root);
        std::fs::create_dir_all(root.join("qc")).unwrap();
    });
}

/// `--fastqc_args` alone activates FastQC even though `cli.fastqc` stays false, which
/// is why `fastqc_active` tests for either flag.
#[test]
fn se_trim_fastqc_args_only() {
    Tripwire::new(
        "se_trim_fastqc_args_only",
        &["--fastqc_args", "--quiet", "a.fastq"],
        &["a.fastq"],
    )
    .run(stage_se);
}

#[test]
fn se_trim_demux() {
    Tripwire::new(
        "se_trim_demux",
        &["--demux", "bc.txt", "a.fastq"],
        &["a.fastq", "bc.txt"],
    )
    .run(|root| {
        stage_se(root);
        write_barcodes(&root.join("bc.txt"));
    });
}

#[test]
fn se_trim_demux_fastqc() {
    Tripwire::new(
        "se_trim_demux_fastqc",
        &["--demux", "bc.txt", "--fastqc", "a.fastq"],
        &["a.fastq", "bc.txt"],
    )
    .run(|root| {
        stage_se(root);
        write_barcodes(&root.join("bc.txt"));
    });
}

#[test]
fn se_trim_no_report_file() {
    Tripwire::new(
        "se_trim_no_report_file",
        &["--no_report_file", "a.fastq"],
        &["a.fastq"],
    )
    .run(stage_se);
}

#[test]
fn se_trim_gzipped_input() {
    Tripwire::new("se_trim_gzipped_input", &["a.fastq.gz"], &["a.fastq.gz"])
        .run(|root| write_fastq_gz(&root.join("a.fastq.gz"), "A_"));
}

#[test]
fn se_trim_clumpify() {
    let out = Tripwire::new(
        "se_trim_clumpify",
        &["--clumpify", "--cores", "2", "a.fastq"],
        &["a.fastq"],
    )
    .run(stage_se);

    // `resolve_clump_layout` falls back to plain mode with a warning and exit 0 when
    // --memory is below the floor, which would make this case observe nothing.
    assert!(
        out.stderr.contains("clumpify:"),
        "clumpify did not engage, so this case proves nothing about its writers:\n{}",
        out.stderr
    );
}

#[test]
fn se_trim_from_ubam_input() {
    Tripwire::new("se_trim_from_ubam_input", &["s.bam"], &["s.bam"])
        .run(|root| copy_fixture("ubam_test.bam", &root.join("s.bam")));
}

// ── site 927 — paired FASTQ trim ─────────────────────────────────────────────

#[test]
fn paired_trim_split_dirs() {
    let out = Tripwire::new(
        "paired_trim_split_dirs",
        &["--paired", "r1/s_R1.fastq", "r2/s_R2.fastq"],
        &["r1/s_R1.fastq", "r2/s_R2.fastq"],
    )
    .run(stage_pair_split);

    // Everything follows R1 (#398), so R2's directory gains nothing.
    let in_r2: Vec<&PathBuf> = out.created.iter().filter(|p| p.starts_with("r2")).collect();
    assert!(in_r2.is_empty(), "R2's directory gained {in_r2:?}");
}

#[test]
fn paired_trim_multipair() {
    Tripwire::new(
        "paired_trim_multipair",
        &[
            "--paired",
            "p1_R1.fastq",
            "p1_R2.fastq",
            "p2_R1.fastq",
            "p2_R2.fastq",
        ],
        &["p1_R1.fastq", "p1_R2.fastq", "p2_R1.fastq", "p2_R2.fastq"],
    )
    .run(|root| {
        for p in ["p1", "p2"] {
            write_fastq(&root.join(format!("{p}_R1.fastq")), "S_");
            write_fastq(&root.join(format!("{p}_R2.fastq")), "S_");
        }
    });
}

/// The only shape where `pair_output_dir`'s `--output_dir` arm and the writers can
/// diverge — i.e. #398's exact failure, with `-o` instead of R1's parent.
#[test]
fn paired_trim_output_dir() {
    Tripwire::new(
        "paired_trim_output_dir",
        &["--paired", "-o", "out", "r1/s_R1.fastq", "r2/s_R2.fastq"],
        &["r1/s_R1.fastq", "r2/s_R2.fastq"],
    )
    .run(stage_pair_split);
}

#[test]
fn paired_trim_retain_unpaired() {
    Tripwire::new(
        "paired_trim_retain_unpaired",
        &["--paired", "--retain_unpaired", "s_R1.fastq", "s_R2.fastq"],
        &["s_R1.fastq", "s_R2.fastq"],
    )
    .run(|root| {
        write_fastq(&root.join("s_R1.fastq"), "S_");
        write_fastq_short(&root.join("s_R2.fastq"), "S_");
    });
}

/// `--cores > 1` swaps in `parallel.rs`'s writer block, which opens its own five
/// handles — a second implementation behind the same candidate list.
#[test]
fn paired_trim_retain_unpaired_cores() {
    Tripwire::new(
        "paired_trim_retain_unpaired_cores",
        &[
            "--paired",
            "--retain_unpaired",
            "--cores",
            "2",
            "s_R1.fastq",
            "s_R2.fastq",
        ],
        &["s_R1.fastq", "s_R2.fastq"],
    )
    .run(|root| {
        write_fastq(&root.join("s_R1.fastq"), "S_");
        write_fastq_short(&root.join("s_R2.fastq"), "S_");
    });
}

#[test]
fn paired_trim_passthrough() {
    Tripwire::new(
        "paired_trim_passthrough",
        &[
            "--paired",
            "--passthrough",
            "s_I1.fastq",
            "s_R1.fastq",
            "s_R2.fastq",
        ],
        &["s_R1.fastq", "s_R2.fastq", "s_I1.fastq"],
    )
    .run(|root| {
        stage_pair(root);
        write_fastq(&root.join("s_I1.fastq"), "S_");
    });
}

#[test]
fn paired_trim_passthrough_cores() {
    Tripwire::new(
        "paired_trim_passthrough_cores",
        &[
            "--paired",
            "--passthrough",
            "s_I1.fastq",
            "--cores",
            "2",
            "s_R1.fastq",
            "s_R2.fastq",
        ],
        &["s_R1.fastq", "s_R2.fastq", "s_I1.fastq"],
    )
    .run(|root| {
        stage_pair(root);
        write_fastq(&root.join("s_I1.fastq"), "S_");
    });
}

#[test]
fn paired_trim_fastqc() {
    Tripwire::new(
        "paired_trim_fastqc",
        &["--paired", "--fastqc", "s_R1.fastq", "s_R2.fastq"],
        &["s_R1.fastq", "s_R2.fastq"],
    )
    .run(stage_pair);
}

#[test]
fn paired_trim_passthrough_fastqc() {
    Tripwire::new(
        "paired_trim_passthrough_fastqc",
        &[
            "--paired",
            "--fastqc",
            "--passthrough",
            "s_I1.fastq",
            "s_R1.fastq",
            "s_R2.fastq",
        ],
        &["s_R1.fastq", "s_R2.fastq", "s_I1.fastq"],
    )
    .run(|root| {
        stage_pair(root);
        write_fastq(&root.join("s_I1.fastq"), "S_");
    });
}

#[test]
fn paired_trim_basename() {
    Tripwire::new(
        "paired_trim_basename",
        &[
            "--paired",
            "--basename",
            "fixed",
            "s_R1.fastq",
            "s_R2.fastq",
        ],
        &["s_R1.fastq", "s_R2.fastq"],
    )
    .run(stage_pair);
}

#[test]
fn paired_trim_no_report_file() {
    Tripwire::new(
        "paired_trim_no_report_file",
        &["--paired", "--no_report_file", "s_R1.fastq", "s_R2.fastq"],
        &["s_R1.fastq", "s_R2.fastq"],
    )
    .run(stage_pair);
}

// ── sites 466 / 497 — hardtrim (outputs are CWD-anchored) ────────────────────

#[test]
fn hardtrim5_fastq() {
    Tripwire::new(
        "hardtrim5_fastq",
        &["--hardtrim5", "20", "in/a.fastq"],
        &["in/a.fastq"],
    )
    .fastqc_incapable()
    .run(|root| write_fastq(&root.join("in/a.fastq"), "A_"));
}

#[test]
fn hardtrim5_ubam() {
    Tripwire::new(
        "hardtrim5_ubam",
        &["--hardtrim5", "20", "--output-format", "ubam", "in/a.fastq"],
        &["in/a.fastq"],
    )
    .fastqc_incapable()
    .run(|root| write_fastq(&root.join("in/a.fastq"), "A_"));
}

#[test]
fn hardtrim3_fastq() {
    Tripwire::new(
        "hardtrim3_fastq",
        &["--hardtrim3", "20", "a.fastq", "b.fastq"],
        &["a.fastq", "b.fastq"],
    )
    .fastqc_incapable()
    .run(stage_se2);
}

#[test]
fn hardtrim3_ubam() {
    Tripwire::new(
        "hardtrim3_ubam",
        &["--hardtrim3", "20", "--output-format", "ubam", "a.fastq"],
        &["a.fastq"],
    )
    .fastqc_incapable()
    .run(stage_se);
}

// ── site 2670 — run_specialty_paired, serving three modes ────────────────────

#[test]
fn clock_multipair() {
    Tripwire::new(
        "clock_multipair",
        &[
            "--clock",
            "--paired",
            "c1_R1.fastq",
            "c1_R2.fastq",
            "c2_R1.fastq",
            "c2_R2.fastq",
        ],
        &["c1_R1.fastq", "c1_R2.fastq", "c2_R1.fastq", "c2_R2.fastq"],
    )
    .fastqc_incapable()
    .run(|root| {
        for p in ["c1", "c2"] {
            write_fastq(&root.join(format!("{p}_R1.fastq")), "S_");
            write_fastq(&root.join(format!("{p}_R2.fastq")), "S_");
        }
    });
}

#[test]
fn implicon_pair() {
    Tripwire::new(
        "implicon_pair",
        &["--implicon", "--paired", "i_R1.fastq", "i_R2.fastq"],
        &["i_R1.fastq", "i_R2.fastq"],
    )
    .fastqc_incapable()
    .run(|root| {
        write_fastq(&root.join("i_R1.fastq"), "S_");
        write_fastq(&root.join("i_R2.fastq"), "S_");
    });
}

#[test]
fn clump_only_paired_split_dirs() {
    Tripwire::new(
        "clump_only_paired_split_dirs",
        &["--clump_only", "--paired", "r1/s_R1.fastq", "r2/s_R2.fastq"],
        &["r1/s_R1.fastq", "r2/s_R2.fastq"],
    )
    .run(stage_pair_split);
}

#[test]
fn clump_only_paired_no_report_file() {
    Tripwire::new(
        "clump_only_paired_no_report_file",
        &[
            "--clump_only",
            "--paired",
            "--no_report_file",
            "s_R1.fastq",
            "s_R2.fastq",
        ],
        &["s_R1.fastq", "s_R2.fastq"],
    )
    .run(stage_pair);
}

// ── sites 661 / 706 / 752 / 809 — clump_only ─────────────────────────────────

#[test]
fn clump_only_se() {
    Tripwire::new(
        "clump_only_se",
        &["--clump_only", "a.fastq", "b.fastq"],
        &["a.fastq", "b.fastq"],
    )
    .run(stage_se2);
}

#[test]
fn clump_only_se_fastqc() {
    Tripwire::new(
        "clump_only_se_fastqc",
        &["--clump_only", "--fastqc", "a.fastq"],
        &["a.fastq"],
    )
    .run(stage_se);
}

#[test]
fn clump_only_ubam_se() {
    Tripwire::new(
        "clump_only_ubam_se",
        &["--clump_only", "--output-format", "ubam", "s.bam"],
        &["s.bam"],
    )
    .run(|root| copy_fixture("ubam_test.bam", &root.join("s.bam")));
}

#[test]
fn clump_only_ubam_paired_interleaved() {
    Tripwire::new(
        "clump_only_ubam_paired_interleaved",
        &[
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "ip.bam",
        ],
        &["ip.bam"],
    )
    .run(|root| copy_fixture("ubam_paired_test.bam", &root.join("ip.bam")));
}

#[test]
fn clump_only_ubam_paired_two_files() {
    Tripwire::new(
        "clump_only_ubam_paired_two_files",
        &[
            "--clump_only",
            "--paired",
            "--output-format",
            "ubam",
            "s_R1.fastq",
            "s_R2.fastq",
        ],
        &["s_R1.fastq", "s_R2.fastq"],
    )
    .run(stage_pair);
}

// ── sites 2047 / 2109 — uBAM output ──────────────────────────────────────────

#[test]
fn ubam_out_se() {
    Tripwire::new(
        "ubam_out_se",
        &["--output-format", "ubam", "a.fastq", "b.fastq"],
        &["a.fastq", "b.fastq"],
    )
    .run(stage_se2);
}

#[test]
fn ubam_out_se_fastqc() {
    Tripwire::new(
        "ubam_out_se_fastqc",
        &["--output-format", "ubam", "--fastqc", "a.fastq"],
        &["a.fastq"],
    )
    .run(stage_se);
}

#[test]
fn ubam_out_paired_two_files() {
    Tripwire::new(
        "ubam_out_paired_two_files",
        &[
            "--paired",
            "--output-format",
            "ubam",
            "r1/s_R1.fastq",
            "r2/s_R2.fastq",
        ],
        &["r1/s_R1.fastq", "r2/s_R2.fastq"],
    )
    .run(stage_pair_split);
}

// ── the two arms with no pre-flight (#414 §2.3) ──────────────────────────────

#[test]
fn orphan_paired_fastq_from_interleaved_bam() {
    Tripwire::new(
        "orphan_paired_fastq_from_interleaved_bam",
        &["--paired", "ip.bam"],
        &["ip.bam"],
    )
    .unplanned()
    .fastqc_incapable()
    .run(|root| copy_fixture("ubam_paired_test.bam", &root.join("ip.bam")));
}

/// Legal here — `--retain_unpaired` is refused only for uBAM *output* — and it adds
/// two more writers to an arm that has no pre-flight, so the full write set is on
/// record for the follow-up issue.
#[test]
fn orphan_paired_fastq_retain_unpaired() {
    Tripwire::new(
        "orphan_paired_fastq_retain_unpaired",
        &["--paired", "--retain_unpaired", "ip.bam"],
        &["ip.bam"],
    )
    .unplanned()
    .fastqc_incapable()
    .run(|root| copy_fixture("ubam_paired_test.bam", &root.join("ip.bam")));
}

#[test]
fn orphan_paired_ubam_from_interleaved_bam() {
    Tripwire::new(
        "orphan_paired_ubam_from_interleaved_bam",
        &["--paired", "--output-format", "ubam", "ip.bam"],
        &["ip.bam"],
    )
    .unplanned()
    .run(|root| copy_fixture("ubam_paired_test.bam", &root.join("ip.bam")));
}

// ── self-tests: the tripwire must be able to fail ────────────────────────────

fn synthetic_planned(site: &str, paths: &[&str]) -> Vec<(String, PathBuf)> {
    paths
        .iter()
        .map(|p| (site.to_string(), PathBuf::from(p)))
        .collect()
}

#[test]
fn comparison_reports_an_unplanned_path() {
    let root = case_root("cmp_unplanned");
    let created: BTreeSet<PathBuf> = ["a_trimmed.fq", "a.fastq_trimming_report.txt"]
        .iter()
        .map(PathBuf::from)
        .collect();
    let planned = synthetic_planned("src/main.rs:1006", &["a_trimmed.fq"]);

    let err = compare(&root, &created, &planned, &["a.fastq"], Arm::Planned, true)
        .expect_err("an unplanned created path must fail");

    assert!(err.contains("a.fastq_trimming_report.txt"), "{err}");
    assert!(err.contains("src/main.rs:1006"), "{err}");
    // The message must steer away from the fix that kills this test.
    assert!(
        err.contains("do NOT narrow") || err.contains("DoNOT") || err.contains("NOT narrow"),
        "{err}"
    );
}

#[test]
fn comparison_accepts_overplanning() {
    let root = case_root("cmp_overplan");
    let created: BTreeSet<PathBuf> = [PathBuf::from("s_R1_val_1.fq")].into_iter().collect();
    let planned = synthetic_planned(
        "src/main.rs:927",
        &["s_R1_val_1.fq", "s_R1_unpaired_1.fq", "s_R2_unpaired_2.fq"],
    );

    compare(&root, &created, &planned, &["--paired"], Arm::Planned, true)
        .expect("planned may legitimately exceed created");
}

#[test]
fn comparison_rejects_an_empty_run() {
    let root = case_root("cmp_empty");
    let err = compare(
        &root,
        &BTreeSet::new(),
        &synthetic_planned("src/main.rs:1006", &["a_trimmed.fq"]),
        &["a.fastq"],
        Arm::Planned,
        true,
    )
    .expect_err("a run that created nothing must fail, not pass vacuously");
    assert!(err.contains("never reached its arm"), "{err}");
}

#[test]
fn comparison_rejects_a_missing_dump() {
    let root = case_root("cmp_nodump");
    let created: BTreeSet<PathBuf> = [PathBuf::from("a_trimmed.fq")].into_iter().collect();
    let err = compare(&root, &created, &[], &["a.fastq"], Arm::Planned, true)
        .expect_err("an empty candidate list must fail on a planned arm");
    assert!(err.contains("hook did not fire"), "{err}");
}

/// T3's standing proof. It is the substantive assertion on the three `Arm::Unplanned`
/// cases, where T1 and T5 are skipped, so it must not be the one assertion whose ability
/// to fail rests on a mutation that gets reverted.
#[test]
fn comparison_detects_a_modified_input() {
    let root = case_root("cmp_clobber");
    let input = root.join("a.fastq");
    write_fastq(&input, "A_");
    let before: BTreeMap<PathBuf, u64> = [(input.clone(), hash_file(&input))].into_iter().collect();

    unchanged(&before).expect("an untouched file must pass");

    write_fastq(&input, "CLOBBERED_");
    let err = unchanged(&before).expect_err("a rewritten pre-existing file must fail");
    assert!(err.contains("a.fastq"), "{err}");
    assert!(err.contains("was modified by the run"), "{err}");

    // A deleted file fails too, rather than being silently skipped.
    std::fs::remove_file(&input).unwrap();
    assert!(
        std::panic::catch_unwind(|| unchanged(&before)).is_err(),
        "a deleted pre-existing file must not pass"
    );
}

/// `every_site_is_reached_by_some_case` proves 11 hand-listed shapes reach 11 *distinct*
/// sites. It cannot prove 11 is *all* of them — a 12th dispatch arm would leave that test
/// green while the new arm has no case and no tripwire. That is this bug family one level
/// up: coverage silently stops extending to a writer.
///
/// Every mention of the pre-flight in `main.rs` is a call, so counting them is a sound
/// oracle. A new arm with **no** pre-flight is still invisible to this — see the module
/// doc's Limits.
#[test]
fn the_site_inventory_is_complete() {
    let src = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("main.rs"),
    )
    .unwrap();
    let calls = src.matches("preflight_output_collisions(").count();
    assert_eq!(
        calls, 11,
        "main.rs has {calls} pre-flight call sites, not 11. If one was added, give it a \
         case and add its shape to every_site_is_reached_by_some_case, then update this \
         count — do not just update the count."
    );
}

#[test]
fn comparison_rejects_two_sites_in_one_run() {
    let root = case_root("cmp_twosites");
    let created: BTreeSet<PathBuf> = [PathBuf::from("a_trimmed.fq")].into_iter().collect();
    let mut planned = synthetic_planned("src/main.rs:1006", &["a_trimmed.fq"]);
    planned.push(("src/main.rs:927".into(), PathBuf::from("b_val_1.fq")));

    let err = compare(&root, &created, &planned, &["a.fastq"], Arm::Planned, true)
        .expect_err("two dumped sites in one run must fail");
    assert!(err.contains("two pre-flights"), "{err}");
}

#[test]
fn comparison_rejects_stray_fastqc_artifacts() {
    let root = case_root("cmp_strayfqc");
    let created: BTreeSet<PathBuf> = ["a_trimmed.fq", "a_trimmed_fastqc.zip"]
        .iter()
        .map(PathBuf::from)
        .collect();
    let planned = synthetic_planned("src/main.rs:1006", &["a_trimmed.fq"]);

    // No FastQC flag, so the allowance must not apply and the artefact is unplanned.
    let err = compare(&root, &created, &planned, &["a.fastq"], Arm::Planned, true)
        .expect_err("a _fastqc artefact with no FastQC flag must fail");
    assert!(err.contains("no FastQC flag"), "{err}");
}

#[test]
fn comparison_rejects_an_allowance_that_takes_everything() {
    let root = case_root("cmp_allswallowed");
    let created: BTreeSet<PathBuf> = [PathBuf::from("a_trimmed_fastqc.zip")]
        .into_iter()
        .collect();
    let planned = synthetic_planned("src/main.rs:1006", &["a_trimmed.fq"]);

    let err = compare(
        &root,
        &created,
        &planned,
        &["--fastqc", "a.fastq"],
        Arm::Planned,
        true,
    )
    .expect_err("an allowance covering every created file leaves T1 vacuous");
    assert!(err.contains("nothing left to check"), "{err}");
}

/// Drives hook → stderr → parse → compare and proves the chain goes red when the dump
/// is absent. Without this, the only proof the wiring is live is a manual edit that
/// gets reverted immediately.
#[test]
fn hook_is_inert_when_env_unset_and_that_goes_red() {
    let root = case_root("hook_gate");
    stage_se(&root);

    let (ok, with_dump) = run_arm(&root, &["a.fastq"]);
    assert!(ok, "run with the dump on failed:\n{with_dump}");
    assert!(
        !parse_planned(&with_dump).is_empty(),
        "the dump produced no candidate lines:\n{with_dump}"
    );

    // Same arm, fresh root, variable unset.
    let root2 = case_root("hook_gate_off");
    stage_se(&root2);
    let out = Command::new(binary())
        .current_dir(&root2)
        .env_remove(DUMP_PLANNED_ENV)
        .arg("a.fastq")
        .output()
        .expect("failed to run trim_galore");
    assert!(out.status.success());
    let without = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(
        parse_planned(&without).is_empty(),
        "the hook printed candidates with the variable unset:\n{without}"
    );

    // And with no dump to compare against, the comparison must fail rather than pass.
    let created = snapshot(&root2);
    compare(&root2, &created, &[], &["a.fastq"], Arm::Planned, true)
        .expect_err("no dump must fail the comparison, not pass it");
}

/// Turns "all 11 sites are covered" from prose into an assertion. Collects the site
/// each representative arm dumps and requires 11 distinct values — no line numbers are
/// pinned, so unrelated edits to `main.rs` cannot break it.
#[test]
fn every_site_is_reached_by_some_case() {
    let root = case_root("site_coverage");

    let shapes: [(&str, &[&str]); 11] = [
        ("se", &["a.fastq"]),
        ("pe", &["--paired", "s_R1.fastq", "s_R2.fastq"]),
        ("ht5", &["--hardtrim5", "20", "a.fastq"]),
        ("ht3", &["--hardtrim3", "20", "a.fastq"]),
        (
            "clock",
            &["--clock", "--paired", "s_R1.fastq", "s_R2.fastq"],
        ),
        ("co_se", &["--clump_only", "a.fastq"]),
        (
            "co_ubam_se",
            &["--clump_only", "--output-format", "ubam", "a.fastq"],
        ),
        (
            "co_ubam_pe1",
            &[
                "--clump_only",
                "--paired",
                "--output-format",
                "ubam",
                "ip.bam",
            ],
        ),
        (
            "co_ubam_pe2",
            &[
                "--clump_only",
                "--paired",
                "--output-format",
                "ubam",
                "s_R1.fastq",
                "s_R2.fastq",
            ],
        ),
        ("ubam_se", &["--output-format", "ubam", "a.fastq"]),
        (
            "ubam_pe2",
            &[
                "--paired",
                "--output-format",
                "ubam",
                "s_R1.fastq",
                "s_R2.fastq",
            ],
        ),
    ];

    let mut sites: BTreeMap<String, &str> = BTreeMap::new();
    for (name, args) in shapes {
        let dir = root.join(name);
        std::fs::create_dir_all(&dir).unwrap();
        write_fastq(&dir.join("a.fastq"), "A_");
        write_fastq(&dir.join("s_R1.fastq"), "S_");
        write_fastq(&dir.join("s_R2.fastq"), "S_");
        copy_fixture("ubam_paired_test.bam", &dir.join("ip.bam"));

        let (ok, stderr) = run_arm(&dir, args);
        assert!(ok, "{name}: exit non-zero\n{stderr}");
        let planned = parse_planned(&stderr);
        assert!(
            !planned.is_empty(),
            "{name}: dumped no candidates\n{stderr}"
        );
        if let Some(prev) = sites.insert(planned[0].0.clone(), name) {
            panic!(
                "{name} and {prev} both dispatched to {} — one pre-flight site is \
                 unexercised, so the coverage claim is untrue",
                planned[0].0
            );
        }
    }

    assert_eq!(
        sites.len(),
        11,
        "expected 11 distinct pre-flight sites, got {sites:?}"
    );
}

/// The tripwire tests the pre-flight's *coverage*; this pins its *purpose*. Two inputs
/// collapsing onto one output under `-o` must be refused before anything is written —
/// the contract #409 violated.
#[test]
fn a_refused_run_writes_nothing() {
    let root = case_root("refused");
    write_fastq(&root.join("d1/x.fastq"), "D1_");
    write_fastq(&root.join("d2/x.fastq"), "D2_");
    let before = snapshot(&root);

    let (ok, stderr) = run_arm(&root, &["-o", "out", "d1/x.fastq", "d2/x.fastq"]);

    assert!(
        !ok,
        "two inputs sharing an output stem must be refused\n{stderr}"
    );
    assert!(
        stderr.contains("Output path collision"),
        "expected a collision refusal:\n{stderr}"
    );
    let created: Vec<PathBuf> = snapshot(&root).difference(&before).cloned().collect();
    assert!(
        created.is_empty(),
        "a refused run must write nothing, found {created:?}"
    );
}
