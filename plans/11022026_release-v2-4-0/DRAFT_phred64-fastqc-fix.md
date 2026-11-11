# Draft — wire `--phred64` into the bundled FastQC config

**Status:** drafted, NOT applied. Needs an implementation trigger.
**Base:** `chore/fastqc-rust-1.1.1` @ `d5f41c5`, after the doc sweep.
**Type:** behaviour fix. Restores the pre-bump report for one flag combination.

## The defect

`fastqc-rust 1.1.1` adds a `phred64: bool` to `FastQCConfig`. `src/fastqc.rs`
builds its config with `..FastQCConfig::default()`, so the field stays `false`,
and nothing maps TrimGalore's own `--phred64` onto it.

FastQC 0.13.0 also changed the encoding heuristic. Measured on
`test_files/phred64_test.fastq` with `--phred64 --fastqc -q 0 --length 0`:

| Build | Encoding reported | Per-base mean, base 1 / base 25 |
|---|---|---|
| 1.0.1 (FastQC 0.12.1) | `Illumina 1.5` | 40.0 / 2.0 — correct |
| 1.1.1 as committed | `Sanger / Illumina 1.9` | **71.0 / 33.0 — 31 too high** |
| 1.1.1 + `phred64: true` | `Illumina 1.5` | 40.0 / 2.0 — correct |

The third row was measured by building `fastqc-rust 1.1.1`'s own `fastqc` binary
from the vendored crate source and running it with and without `--phred64` on the
same trimmed file, so the fix is confirmed rather than inferred.

Why it bites: TrimGalore passes FASTQ quality strings through verbatim —
`cli.phred_offset()` feeds only `BamWriter`, never the FASTQ writer — so under
`--phred64` the trimmed FASTQ on disk really is Phred+64. The old pin read that
correctly by heuristic; the new pin has to be told.

## The non-obvious part: BAM output must NOT get the flag

`--phred64` is rejected for BAM *input* (`main.rs:343`) but not for BAM *output*.
On the uBAM-output path `BamWriter` subtracts the input offset, so BAM `QUAL`
holds raw Phred by spec — there is no ASCII offset left to declare. Passing
`phred64: true` to a FastQC run over a `.bam` would subtract 64 a second time.

So the mapping is **not** `cli.phred64` uniformly. Classifying all 11 call sites:

| Call site | Output | `phred64` |
|---|---|---|
| `main.rs:1564` | FASTQ SE | `cli.phred64` |
| `main.rs:1881`, `:1887` | FASTQ PE | `cli.phred64` |
| `main.rs:1899` | FASTQ passthrough carrier | `cli.phred64` |
| `main.rs:2353`, `:2511`, `:2630` | **BAM** | `false` |
| `clump_only.rs:420` | FASTQ SE | threaded `phred64` |
| `clump_only.rs:610`, `:611` | FASTQ PE | threaded `phred64` |
| `clump_only.rs:914`, `:1173` | **BAM** | `false` |

Specialty modes (`--hardtrim5/3`, `--clock`, `--implicon`) need no change —
`--fastqc` is already refused there under #421, so they reach no call site.

Rather than trust eleven call sites to each pick the right constant, `run()`
zeroes the flag for `.bam` itself. That puts the invariant on the thing it
constrains and makes every call site pass the same expression.

---

## Edits

### 1. `src/fastqc.rs` — `run()` gains the parameter and the BAM guard

```diff
 /// `cores`       — threading budget for `fastqc-rust`'s internal rayon
 ///                 pool. Always at least 1.
+/// `phred64`     — input quality encoding declared by `--phred64`. Ignored
+///                 for BAM output, which stores raw Phred.
 /// Threads `run` will give FastQC: `--fastqc_args -t N` overrides `cores`.
```

```diff
 pub fn run(
     output_path: &Path,
     fastqc_args: Option<&str>,
     output_dir: Option<&Path>,
     cores: usize,
+    phred64: bool,
 ) -> Result<()> {
+    // BAM QUAL stores raw Phred, so the ASCII offset never applies there.
+    let phred64 = phred64
+        && !output_path
+            .extension()
+            .is_some_and(|e| e.eq_ignore_ascii_case("bam"));
+
     let mut config = FastQCConfig {
         output_dir: output_dir.map(|p| p.to_path_buf()),
         threads: Some(cores.max(1)),
         // trim-galore prints "Running FastQC on …" itself; suppress
         // fastqc-rust's per-file progress chatter to avoid
         // double-logging.
         quiet: true,
+        phred64,
         ..FastQCConfig::default()
     };
```

### 2. `src/main.rs` — pass `cli.phred64` at all seven sites

Each of `:1564`, `:1881`, `:1887`, `:2353`, `:2511`, `:2630` gains a final
argument:

```diff
         fastqc::run(
             &output_path,
             cli.fastqc_args.as_deref(),
             output_dir,
             cli.cores,
+            cli.phred64,
         )?;
```

and `:1899`:

```diff
-            fastqc::run(pt_out, cli.fastqc_args.as_deref(), output_dir, cli.cores)?;
+            fastqc::run(
+                pt_out,
+                cli.fastqc_args.as_deref(),
+                output_dir,
+                cli.cores,
+                cli.phred64,
+            )?;
```

The three BAM sites pass `cli.phred64` like the rest; `run()` discards it. Their
existing "relies on output_path ending in `.bam`" comments already state the
precondition the guard keys on.

### 3. `src/clump_only.rs` — thread the flag through

Four public functions gain `phred64: bool`, placed immediately after
`fastqc_args` so the FastQC arguments stay one visual group:

- `clump_only_single` (`:285`)
- `clump_only_paired` (`:432`)
- `clump_only_single_to_bam` (`:790`)
- `clump_only_paired_to_bam_one_pair` (`:951`)

```diff
     fastqc_requested: bool,
     fastqc_args: Option<&str>,
+    phred64: bool,
     no_report_file: bool,
 ) -> Result<ClumpOnlyStats> {
```

and the five `fastqc::run` calls inside them gain `phred64`. The five call sites
in `main.rs` (`:689`, `:725`, `:769`, `:832`, `:873`) pass `cli.phred64`.

**Risk, stated plainly:** this makes a fourth adjacent `bool` in signatures that
already carry `gzip_output`, `fastqc_requested` and `no_report_file`. A
transposition between `phred64` and `no_report_file` compiles silently. The test
below catches the FASTQ-path case; if you would rather not add the fourth bool,
the alternative is grouping the FastQC arguments into a small struct, which is a
wider diff than a release fix wants. Flagging it rather than deciding it.

### 4. `Cargo.toml` — `zip` as a dev-dependency

Needed to read `fastqc_data.txt` out of the report. Already in the graph via
`fastqc-rust`, so no new tree — same argument as the existing `bstr` entry.

```diff
 # Already in the dep graph (transitive via noodles); promoted to a direct
 # dev-dep so test code can construct `bstr::BStr` for noodles-aux-Value
 # variants that take `&'a BStr` (e.g. `Value::Hex`).
 bstr = "1"
+# Already in the dep graph (transitive via fastqc-rust); promoted to a direct
+# dev-dep so test code can read fastqc_data.txt out of the report archive.
+zip = "9"
```

### 5. `tests/integration_fastqc_phred64.rs` — new file

Every existing FastQC test asserts only that `*_fastqc.zip` exists; none reads
inside it. That is the gap this defect went through, so the test opens the report.

```rust
//! `--phred64` must reach the bundled FastQC config.
//!
//! FastQC 0.13.0 keys its encoding guess off the observed quality characters, so
//! a Phred+64 file read as Phred+33 reports every mean 31 too high. TrimGalore
//! writes FASTQ quality verbatim, so `--phred64` output is Phred+64 on disk and
//! the engine has to be told. BAM output stores raw Phred and must not be.

use std::io::Read;
use std::process::Command;

use tempfile::tempdir;

fn report_field(zip_path: &std::path::Path, stem: &str, measure: &str) -> String {
    let file = std::fs::File::open(zip_path).expect("report archive");
    let mut archive = zip::ZipArchive::new(file).expect("valid zip");
    let mut data = String::new();
    archive
        .by_name(&format!("{stem}_fastqc/fastqc_data.txt"))
        .expect("fastqc_data.txt in archive")
        .read_to_string(&mut data)
        .expect("readable");
    data.lines()
        .find(|l| l.starts_with(measure))
        .unwrap_or_else(|| panic!("no {measure} line in:\n{data}"))
        .split('\t')
        .nth(1)
        .expect("tab-separated value")
        .to_string()
}

#[test]
fn phred64_reaches_the_fastqc_config_on_fastq_output() {
    let dir = tempdir().unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_trim_galore"))
        .args([
            "--phred64",
            "--fastqc",
            "-q",
            "0",
            "--length",
            "0",
            "--output_dir",
            dir.path().to_str().unwrap(),
            "test_files/phred64_test.fastq",
        ])
        .status()
        .expect("trim_galore runs");
    assert!(status.success());

    let zip = dir.path().join("phred64_test_trimmed_fastqc.zip");
    assert_eq!(
        report_field(&zip, "phred64_test_trimmed", "Encoding"),
        "Illumina 1.5",
        "a Phred+64 file read as Phred+33 inflates every mean by 31"
    );
}
```

A companion case asserting the `.bam` path still reports `Sanger / Illumina 1.9`
under `--phred64 --output-format ubam` would pin the other half of the guard.
Drafting that needs a check I have not run — whether `--phred64` with FASTQ input
and uBAM output is reachable in a single invocation — so it is left out rather
than written blind.

---

## Proposed commit message

```
fix(fastqc): declare --phred64 to the bundled FastQC engine

fastqc-rust 1.1.1 added a phred64 field to FastQCConfig and FastQC 0.13.0
changed the encoding heuristic. TrimGalore writes FASTQ quality verbatim,
so --phred64 output is Phred+64 on disk, and the engine now reads it as
Phred+33: on test_files/phred64_test.fastq the per-base means come out
71 and 33 where the data is Q40 and Q2. The 0.12.1 heuristic guessed the
encoding correctly, so this path was right before the pin moved.

BAM output is excluded. BamWriter already subtracts the input offset, so
BAM QUAL holds raw Phred and a second subtraction would repeat the bug in
the other direction. run() zeroes the flag for .bam itself rather than
asking eleven call sites to each pick the right constant.

The report archive is now opened in a test. Every FastQC test so far
asserted only that the .zip exists, which is how this reached a release
candidate with 824 passing tests.
```

---

## Follow-on

With this in place, `test_files/README.md:14` (held back from the doc sweep) can
say that TrimGalore declares the encoding, so the heuristic no longer applies to
its own FastQC runs — and the fixture's two distinct quality values remain
necessary for the writer-offset assertion, which is what that note is really for.
