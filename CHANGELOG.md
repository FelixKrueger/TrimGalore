# Trim Galore Changelog


### Unreleased (queued for v2.1.0-beta.6)

#### Performance (since v2.1.0-beta.5) — Buckberry-scale audit (#248)

Profiled, prototyped, and benchmarked by @an-altosian against an
84M-read single-end Bisulfite-Seq fixture (Buckberry et al. 2023,
SRR24827373, 2.1 GiB gzipped, 38% adapter rate) using `hyperfine`
with 10 trials per condition. Two of the three confirmed wins from
that audit landed here; item #4 (Myers' bit-parallel adapter
alignment, ~13% additional reduction) is queued for a separate PR
from @an-altosian.

- **Output gzip compression level lowered 6 → 1.** Compression CPU
  was the dominant wall-time consumer on saturated workers
  (cores=8) — at level 6, gzip output dwarfed actual trimming work.
  Lowering to level 1 (fastest) measured −23.3% wall-clock
  (69.693 ± 0.200 s → 53.471 ± 2.135 s) and −43% user-CPU
  (593.9 s → 338.7 s) at Buckberry scale. Trade: output `.fq.gz`
  files are roughly 75% larger (273 KB → 479 KB on the local
  BS-seq_10K_R1 smoke). Decompressed content is byte-identical to
  the previous output — the CI validation matrix's `gzip -dc |
  md5sum` comparisons against Perl v0.6.11 still pass on every
  flag path. New `pub const OUTPUT_GZIP_LEVEL` in `src/fastq.rs`
  centralises the level so a future `--high-compression` opt-in
  flag is a one-line change for storage-conscious users. (#248 #1)
- **Single buffered write per FASTQ record.** Pre-format the
  4-line record into a `Vec<u8>`, then issue one `write_all`
  instead of four separate `writeln!` calls. Byte-identical output
  (md5 match verified end-to-end on BS-seq_10K_R1). At Buckberry
  scale this measured 9.8% wall-clock reduction
  (69.693 ± 0.200 s → 62.862 ± 0.352 s) by amortising the per-call
  Write-trait overhead. (#248 #2)

#### Bug fixes (since v2.1.0-beta.5) — Perl-parity regressions (contributor-reported)

- **`--basename foo --paired` now produces `foo_val_1.fq.gz` /
  `foo_val_2.fq.gz`, not `foo_R1_val_1.fq.gz` / `foo_R2_val_2.fq.gz`.**
  The `_R{1,2}` segment was interpolated between the basename and the
  `_val_{1,2}` suffix in `io::paired_end_output_names` and
  `io::unpaired_output_names`, breaking the documented Perl v0.6.5+
  filename contract — every nf-core / Snakemake pipeline globbing
  the documented `${basename}_val_*` path silently missed outputs
  under v2.1.0-beta.5. Trimmed read content was unaffected, only the
  filenames differed. Two regression tests added covering the
  basename branches of both functions. Reported by @an-altosian
  during the Phase-1C parity hunt. (#244)

- **Lowercase clip-flag spellings (`--clip_r1`, `--clip_r2`,
  `--three_prime_clip_r1`, `--three_prime_clip_r2`) are now accepted.**
  Perl `trim_galore` accepted both lowercase and uppercase spellings;
  the Rust port matched the uppercase canonical only and rejected the
  lowercase forms with a clap parse error (exit 2). Every Perl-era
  pipeline using the documented lowercase spelling broke under v2.x.
  Added lowercase aliases on all four flags. Reported and diagnosed by
  @an-altosian during a Phase-1B Perl-parity hunt. (#242)

- **Output gzip compression now mirrors input compression by default.**
  Plain `.fastq` input → plain `.fq` output; `.fastq.gz` → `.fq.gz`.
  Restores Perl v0.6.x behaviour. Rust v2.1.0-beta.5 always emitted
  gzipped output regardless of input, breaking pipelines that globbed
  `*.fq` (no `.gz`) for outputs from plain-text inputs. `--dont_gzip`
  still works as the explicit "always plain" override. The first input
  determines the mode for the whole run; mixing plain and gzipped
  inputs in one invocation isn't a supported configuration. Reported
  by @an-altosian via #245 (item A).

- **`--retain_unpaired` now routes both mates independently to their
  unpaired files when both mates fail the discard `--length` cutoff,
  if each is individually long enough for the per-side `--length_{1,2}`
  threshold.** Rust v2.1.0-beta.5 had an extra `!r{1,2}_short` clause in
  `filters::filter_paired_end` that gated unpaired rescue on the read
  itself passing the discard cutoff, which silently dropped reads when
  both mates failed `--length` together but were individually long
  enough. Matches Perl `master:trim_galore:2325-2343`. Slight caveat:
  Perl's behaviour diverges from its own user-guide wording ("rescue
  the surviving mate"); we match the implementation, not the docs, to
  preserve byte-identity for the documented byte-identity flag paths.
  Regression test added. Reported by @an-altosian via #245 (item C).

#### Behavioural notes (v2.x intentional widenings, since v2.1.0-beta.5)

- **`--clock` and `--implicon` now imply `--paired`** — passing either
  flag without `--paired` is no longer rejected. Both modes are
  inherently paired-end specialty modes; requiring users to pass
  `--paired` redundantly was noise. Pipelines using the explicit Perl
  form (`--clock --paired` / `--implicon --paired`) continue to work
  unchanged. Consistent with the multi-pair widening pattern documented
  for these specialty modes in beta.4. Surfaced by @an-altosian via
  #245 (item D).

#### Bug fixes (since v2.1.0-beta.5) — Perl-parity regressions (contributor-reported, continued)

- **`--max_n` fraction-mode now logs the Perl-style notice on entry**
  ("`--max_n will be interpreted as a fraction of the read length
  (0.5)`"). Investigating @an-altosian's #243 confirmed the dispatch
  and filter chain behave correctly: values in `(0.0, 1.0)` already
  build `MaxNFilter::Fraction` and `n_count/length > threshold`
  filtering matches Perl v0.6.8+ behaviour byte-for-byte. The
  reproducer-fixture's max-N fraction (3/153 ≈ 0.02) just doesn't
  exceed the 0.5 threshold, so neither implementation filters
  anything — *not* a bug. Adding the same warning Perl prints
  (`master:trim_galore:3328`) makes the selected mode visible at
  runtime so users don't have to derive it from output statistics.
  (#243)

#### Bug fixes (since v2.1.0-beta.5)

- **`-o/--output_dir DIR` no longer hangs when `DIR` doesn't exist.** The
  parallel paired-end path opened output files via raw `File::create`,
  which fails immediately on a missing parent — but by then reader and
  worker threads were already spawned, the `?` exit dropped the receiver
  channel, and the process deadlocked at near-zero CPU (workers stuck
  producing into a queue with no consumers). Reported via beta.5 user
  feedback (24h wall / 4s CPU on a SLURM cluster). Fix: hoist
  `create_dir_all` into `main()` immediately after CLI parse, covering
  every downstream code path (parallel, single-threaded, paired,
  single-end, every specialty mode) in one place. Restores Perl v0.6.x
  behaviour ("If an output directory which was specified with -o
  output_directory did not exist, it will be created for you", v0.6.0
  changelog).

#### Tests (since v2.1.0-beta.5) — coverage gaps closed (#246)

Five new unit tests landed across `report.rs`, `demux.rs`, `fastq.rs`,
and `adapter.rs`. Closes 5 of the 6 `§5.x` items from the test-coverage
audit. Total test count: 171 → 177. Items still open from #246:
parallel/serial stat-tracking parity (§5.2), comprehensive
`parallel.rs` coverage, optional upstreaming of @an-altosian's
proptest harness — tracked as separate followups.

- **§5.4 PE param-summary `removed-end:` regression guard**
  (`report.rs::tests`). Beta.3 fixed a stray `-end` suffix in the
  paired-end parameter-summary line (`...sequence pair gets
  removed-end: 20 bp` → `...removed: 20 bp`); MultiQC parsers grep
  for the literal `removed:` form. Test renders the PE header and
  asserts both `!contains("removed-end")` and the canonical phrase —
  any reintroduction of the typo class fails the assertion.

- **§5.5 Demux CRLF samplesheet handling** (`demux.rs::tests`).
  Windows-authored barcode sheets use `\r\n` line endings; without
  the `trim_end_matches('\r')` strip in `read_barcode_file`, the
  trailing `\r` would pollute the parsed barcode and fail the
  ACGTN-only validator with a confusing "barcode must contain only
  A, C, T, G, N" error. Test writes a CRLF samplesheet and asserts
  no stray `\r` survives on any parsed entry.

- **§5.6 Demux short-read NoCode routing** (`demux.rs::tests`).
  When a trimmed read is shorter than the barcode length,
  `demultiplex` (`src/demux.rs:178-187`) routes it to the NoCode
  bucket instead of slicing past the read end. End-to-end test:
  fixture with one 5 bp read, one non-matching 16 bp read, and one
  matching 16 bp read; asserts both NoCode-bound reads land in
  `*_NoCode.fq` (with cleared seq+qual for the too-short one) and
  the matching read lands in the per-sample bucket. Together with
  §5.5 closes the `demux.rs` zero-tests module gap.

- **§5.1 Multi-member gzip reader round-trip** (`fastq.rs::tests`).
  The parallel writer (`--cores N`) emits each worker's chunk as
  its own gzip member; the concatenated stream is a valid RFC 1952
  multi-member gzip file. Test crafts a 2-member buffer with
  `GzEncoder::finish` twice, concatenates, and asserts
  `FastqReader` (using `MultiGzDecoder`) yields records from BOTH
  members in order. Originally fixed in `9dcf519` (pre-beta.1) but
  never had a unit-level regression test.

- **§5.3 Adapter auto-detect `MAX_SCAN_READS` boundary**
  (`adapter.rs::tests`). The 1M-read scan cap was unverified at the
  unit level; generating a >1M-read fixture per test run is too
  slow. Refactored: extracted `autodetect_adapter_with_max_scan`
  (crate-private) so tests can exercise the same `increment-then-
  break` control flow at scale 7. Two paired tests: cap-bounded
  (input has 100 records, max_scan=7, asserts `reads_scanned == 7`
  and matches < 100) and cap-unbounded (input has 50 records,
  max_scan=1M, asserts full-file scan). Public `autodetect_adapter`
  API unchanged.

#### Tests (since v2.1.0-beta.5) — parallel/serial stats parity (#246 §5.2)

- **`parallel::run_single_end_parallel` and `trimmer::run_single_end`
  must yield field-identical `TrimStats` on the same input.** First
  unit test in `src/parallel.rs` (closes the zero-tests module gap
  noted in #246). Beta.0/1 had per-field stat drift between the two
  paths (commits 82d1e34, 3996fc5 fixed `total_bp_after_trim` /
  `rrbs_r2_clipped_5prime`); this test locks the invariant down at
  the unit level. `TrimStats` gained a `PartialEq` derive so a
  single `assert_eq!` covers every field — any future field added
  to the struct is automatically covered without test edits.
  Closes #246 §5.2.

#### Infrastructure (contributor-facing, since v2.1.0-beta.5) — CI hardening (#247)

The five remaining deferred items from the original CI audit landed
in this round (items 2, 3, 4, 7 — item 8 cargo-nextest deferred
pending a focused per-test-isolation audit since some existing tests
share `std::env::temp_dir().join(...)` paths):

- **#247 item 2 — macOS matrix on `rust-tests`.** Job now runs on
  both `ubuntu-latest` and `macos-latest` (Apple Silicon hosted
  runner). Catches Apple-Silicon-specific regressions before
  release-tag time. `fail-fast: false` so an OS-specific failure on
  one entry doesn't kill the other.
- **#247 item 3 — release-profile test step.** New `Run tests
  (release)` step alongside the existing debug-profile run. Catches
  LTO + `codegen-units=1` interactions that the default debug build
  doesn't see. Cheap because the next step (`cargo build --release`)
  was already populating the same target dir.
- **#247 item 4 — line/branch coverage reporting via
  `cargo-llvm-cov`.** New `coverage` job emits an LCOV file as a
  CI artifact (downloadable from the Actions run UI) plus a text
  summary in the job log. No third-party uploader integration —
  Codecov / Coveralls is a separate decision.
- **#247 item 7 — `justfile` for local CI parity.** New top-level
  `justfile` with targets `fmt`, `clippy`, `test`, `test-release`,
  `ci`, `reproduce`, `validate-paired-end`, `logos`, `docs`. Run
  `just` (or `just ci`) to execute the same portable checks CI
  runs on every push. Pairs cleanly with the contributor docs
  flow.



Three CI improvements landed from @an-altosian's audit (#247). Touches
only `.github/workflows/ci.yml`; no runtime change.

- **Validation outputs uploaded on failure.** When any md5 oracle
  step fails, the `/tmp/tg`, `/tmp/op*`, and `/tmp/*.log` paths that
  triggered the mismatch were previously lost when the runner cleaned
  up. New `if: failure()` step uploads them as a 7-day artifact under
  `validation-outputs-<run_id>-<attempt>`. Especially useful for
  reviewing perf PRs that intentionally change output bytes (e.g. a
  default gzip-level change) — the new artefacts can be diffed
  against the Perl baseline directly.
- **Perl Trim Galore source fetched from local `master` instead of
  `raw.githubusercontent`.** The Perl v0.6.x release line lives at
  `master:trim_galore` in this repo; replacing the curl with
  `git fetch --depth=1 origin master && git show
  origin/master:trim_galore` gives byte-identical content with zero
  external network dependency, eliminating a class of CI flake.
- **Cutadapt bioconda revision pinned (`cutadapt=5.2=*_0`).** The
  validation matrix uses Cutadapt's output as the Perl-side oracle,
  so an unannounced bioconda revision bump (5.2-1, etc.) could
  silently shift the md5 baseline. Pin to the first build of 5.2 so
  any rev bump becomes a visible CI failure rather than invisible
  drift.

#### Bug fixes (since v2.1.0-beta.5) — surfaced by the nf-core pre-GA validation review

- **RRBS samples: `Total written (filtered)` cutadapt-section line now matches
  Perl v0.6.x byte-for-byte.** Beta.4's `eedbc66` MultiQC-parity fix introduced
  `total_bp_after_trim` (incremented after the full per-read trimming pipeline)
  but didn't account for the `--rrbs` 2 bp 3' truncation. The reported value
  drifted from v0.6.x by `RRBS-trimmed-reads × 2 bp`. `TrimResult` now carries
  a `bp_after_cutadapt` field captured immediately after quality + adapter
  trimming and before RRBS / poly-A/G / N-trim / clipping, and that's what
  drives `total_bp_after_trim` in both single-end and paired-end pipelines.
  Trimmed FASTQ output is unchanged — this fix only affects the reported
  count. (#232)
- **`RUN STATISTICS` filter-removed lines are now always emitted, even when
  the count is 0.** The `if stats.too_short > 0` / `too_long` / `too_many_n`
  guards (and the PE counterpart `pairs_removed_n` line) caused MultiQC's
  canonical fallback parser to break on samples that pass 100% of reads
  through length / max-N filters — the parser greps for the exact line and
  treats absence as a parse failure. v0.6.x always emits the line. Now we do
  too, with zero-protected percentage display. (#233)

#### Documentation (since v2.1.0-beta.5)

- **Migration notes: trimming-report behaviour changes vs Perl v0.6.x.** The
  pre-GA review surfaced four intentional report-text changes that were not
  filed because they match the v2.x reference report attached to MultiQC #3529.
  Documented for users / parsers expecting the v0.6.x shape:
  - **RRBS quality-trim line shape** changed from `Sequences were truncated to
    a varying degree because of deteriorating qualities …: N (P%)` (counts
    *reads*, v0.6.x RRBS only) to the cutadapt-style `Quality-trimmed: N bp
    (P%)` (counts *bp*, v2.x both modes). Non-RRBS mode is unchanged in both
    implementations. v2.x is more consistent across modes, but a regex tuned
    to one shape won't match the other. (#234)
  - **Adapter family-name annotation** dropped — v2.x emits the bare
    sequence (e.g. `'AGATCGGAAGAGC'`) where v0.6.x emitted family names
    (Illumina TruSeq, Nextera, smallRNA). The family is still tracked
    internally but not rendered in the report.
  - **"Bases preceding removed adapters" histogram** omitted from the
    `=== Adapter N ===` block.
  - **Per-adapter "Minimum overlap" line** not repeated under each adapter
    block (the same datum is in the parameter summary at the top).
  - **Length-distribution `max.err` column** uses the modern Cutadapt formula
    `floor(L × error_rate)`. v0.6.x's display capped this at 1 in many
    positions. The `count` column is unchanged byte-for-byte.


### Version 2.1.0-beta.5 (Release on 27 Apr 2026)

#### Bug fixes (since v2.1.0-beta.4)
- **Bundled FastQC: percentage precision in `>>Overrepresented sequences`.**
  Bumps `fastqc-rust` dep from v1.0.0 to v1.0.1, which restores Java
  FastQC 0.12.1 byte-identity in that section. v1.0.0 rounded the
  percentage column to 2 decimals (`7.16`) instead of emitting Java's
  full `Double.toString()` precision (`7.160449112640348`). Detected
  during the `nf-core/rnaseq#1789` integration matrix on NF 25.04.3,
  where the older pinned MultiQC preserves the literal percentage
  string when aggregating into `fastqc_trimmed_overrepresented_sequences_plot.txt`,
  so the truncated value cascaded into a downstream MD5 mismatch. Fix
  filed and merged upstream as
  [ewels/FastQC-Rust#2](https://github.com/ewels/FastQC-Rust/pull/2).
  Detected sequences, counts, and source classification were always
  correct — this was a cosmetic formatting deviation, not a scientific
  one.


### Version 2.1.0-beta.4 (Release on 26 Apr 2026)

#### New features (since v2.1.0-beta.3)
- **Bundled FastQC.** `--fastqc` now uses the
  [fastqc-rust](https://crates.io/crates/fastqc-rust) library directly
  instead of shelling out to an external `fastqc` binary. Removes Java
  and the FastQC tarball as runtime dependencies — completes the
  "single static binary, zero external runtime deps" story for v2.x.
  Output files (`*_fastqc.html`, `*_fastqc.zip`) are FastQC 0.12.1-
  compatible (the same version we previously bundled in the Docker
  image), so MultiQC parsers and downstream pipelines see identical
  structure. The Docker image is correspondingly slimmer (no
  `default-jre-headless`, no `perl`, no FastQC tarball; saves
  approximately 350 MB at the runtime layer). (#226)
  - `--fastqc_args` continues to accept the common subset (`--nogroup`,
    `--expgroup`, `--quiet`, `--svg`, `--nano`, `--nofilter`,
    `--casava`, `-t`/`--threads`, `-o`/`--outdir`); other flags emit a
    warning and are ignored — forward-compat with future fastqc-rust
    additions.
  - `--help` text for `--fastqc` and `--fastqc_args` refreshed to
    describe the bundled integration and enumerate the translated flag
    set; `docs/SUMMARY.md` architecture-shift paragraph and parity
    table updated accordingly. (#227)

#### Bug fixes (since v2.1.0-beta.3)
- `--clock` and `--implicon` now accept multi-pair input (an even
  number of files as consecutive R1/R2 pairs), restoring v0.6.x
  semantics that the v2.x rewrite had narrowed to "exactly 2 input
  files". Same widening as the `--paired` fix in beta.2 — the two
  specialty run-and-exit modes had their own validation that wasn't
  updated at the time. Each pair gets a per-pair header
  (`=== Clock pair N of M ===` / `=== IMPLICON pair N of M ===`) and
  the same output-collision pre-flight (case-insensitive on full path)
  that `--paired` runs. (#224)

#### Infrastructure (contributor-facing)
- `rust-version` bumped from 1.85 → 1.88 (required by fastqc-rust).
- `.gitattributes` added so the GitHub repo language bar reflects the
  actual Rust content rather than HTML in `docs/`. (#225, contributed
  by @ewels)


### Version 2.1.0-beta.3 (Release on 24 Apr 2026)

#### New features (since v2.1.0-beta.2)
- **BGI/DNBSEQ added to adapter auto-detection.** Users running
  BGI/MGI/DNBSEQ data no longer need to pass `--bgiseq` explicitly; the
  32 bp BGI adapter is now probed alongside Illumina, Nextera, and
  smallRNA on the first 1 M reads. Stranded Illumina stays
  explicit-only (`--stranded_illumina`) because its sequence differs
  from Nextera by a single A-tail base and would produce constant
  ambiguous ties if probed. Tie-break semantics unchanged — the
  zero-count fallback is still Illumina.
- **Repeatable `-a` / `-a2` multi-adapter syntax.** `-a SEQ1 -a SEQ2`
  now works directly — no need for the v0.6.x embedded-string
  (`-a " SEQ -a SEQ"`) or FASTA file workaround. Embedded-string and
  `file:adapters.fa` forms still work and can be mixed with repeated
  flags in a single invocation (e.g. `-a " SEQ1 -a SEQ2" -a SEQ3`).
- **Perl-style `A{N}` single-base expansion** for `-a` / `-a2` — e.g.
  `-a A{10}` expands to `-a AAAAAAAAAA`, matching Perl v0.6.x syntax.
  Only applied to the single-adapter case (not to multi-adapter or
  FASTA entries), mirroring Perl behaviour.
- **Perl-era `-r1` / `-r2` / `-a2` short-flag forms are accepted.**
  Clap's single-character short-flag rule meant `-r1 40` previously
  parsed as `-r=1` with `40` becoming a stray positional (producing
  a confusing "odd count" error). A small pre-parse hook now
  transparently rewrites the exact tokens `-r1`, `-r2`, `-a2` (and
  their `=VALUE` variants) to the existing `--r1`, `--r2`, `--a2`
  long-alias forms so Perl-era scripts keep working.

#### Bug fixes (since v2.1.0-beta.2)
- `--trim-n` is now suppressed under `--rrbs`, restoring Perl v0.6.x
  byte-identical output for users who combine both flags. Perl's RRBS
  code path omitted `$trim_n` from its Cutadapt invocations; beta.2
  was applying N-trimming unconditionally, which narrowly violated the
  byte-identity guarantee for that specific flag combination.
- The paired-end parameter-summary line in the text trimming report
  previously emitted a stray `-end` suffix
  (`...before a sequence pair gets removed-end: 20 bp`). Now correctly
  emits `...before a sequence pair gets removed: 20 bp`. Single-end
  output (`...length single-end: 20 bp`) is unchanged.

#### Documentation
- **Flag-by-flag help-text polish (#221).** 25 docstring edits across
  `src/cli.rs`. Most notable: `--paired` no longer claims "exactly 2
  input files" (the multi-pair fix from beta.2 made this stale);
  `--rrbs` help now warns against Tecan Ovation kit incompatibility;
  `--small_rna` surfaces its length auto-lowering side-effect;
  `--bgiseq` notes it is also probed by auto-detection.
- **Positioning reframed from "byte-identical" to "faithful rewrite
  with useful additions" (#222).** The original framing no longer held
  given new capabilities (poly-G handling, generic poly-A trimmer,
  per-pair adapter detection, the above BGI auto-detect, etc.).
  Updated in README, SUMMARY, User Guide, CHANGELOG, and the
  `--help` preamble. Benchmarks' "byte-identical across core counts"
  claims (which are about `--cores` determinism, not Perl equivalence)
  are retained where accurate.
- **User guide refreshed (#223).** 534 → 325 lines. Dropped a ~220-line
  duplicate of `--help` that had been drifting out of sync on every
  polish cut; replaced with a curated "Flag reference" section on
  cross-flag interactions, RRBS-specific guidance (Tecan, MseI), and
  adapter-specification recap. Added IMPLICON coverage (missing from
  the original guide). Modernised the intro/framing, removed the
  floating Babraham logo (project now maintained solo outside
  Babraham), dropped the hand-maintained "Last update" line, and
  renamed the `Version 0.6.11` section heading to `Introduction`.
  Historical attribution preserved as a footer credit.

#### Infrastructure
- Test count grew from **106 → 147** across the beta.2→beta.3 window.
  New coverage: multi-pair validation branches, specialty modes
  (`--clock`, `--implicon` end-to-end), adapter brace expansion,
  four-probe auto-detection set, Perl-era flag rewriting,
  `--trim-n`/`--rrbs` interaction, `parse_adapter_specs` mixed-form
  path.


### Version 2.1.0-beta.2 (Release on 20 Apr 2026)

#### New features (since v2.1.0-beta.1)
- `--version` now prints build provenance on a second line: `<git-hash> — <target>
  — built <ISO-8601 UTC timestamp>`. The short form `-V` remains unchanged (one
  line, matches the original terse format). Useful for bug reports — users can
  paste `trim_galore --version` to pinpoint the exact build.
- Builds are now **reproducible**: setting `SOURCE_DATE_EPOCH` to a fixed Unix
  timestamp (Debian reproducible-builds spec) produces a bit-identical binary
  across runs. Unset, the build stamps the current wall-clock time as before.
  Malformed values hard-fail the build rather than silently falling back.

#### Infrastructure (contributor-facing, no runtime effect)
- New CI gates on every PR: `cargo fmt --check` + `cargo clippy -D warnings`
  (lint), a dedicated reproducibility job that builds the release binary twice
  under a fixed `SOURCE_DATE_EPOCH` and asserts bit-identity, and a weekly
  `rustsec/audit-check` for dependency advisories.
- Dependabot enabled for cargo + github-actions ecosystems (weekly, Monday,
  limit 5, routed to `@FelixKrueger`).

#### Bug fixes (since v2.1.0-beta.1)
- `--paired` now accepts any even number of input files and processes them as
  consecutive R1/R2 pairs, restoring v0.6.x Perl behaviour. Beta 1 rejected
  more than 2 files with "Paired-end mode requires exactly 2 input files".
  Common shell-glob invocations like `trim_galore --paired *fastq.gz` now
  work again. Adapter auto-detection and poly-G scanning run **per pair**
  (intentional deviation from Perl v0.6.x, which detected once on the first
  input file). This is safer for shell-glob invocations that mix library
  types or 2-colour/4-colour chemistries across samples, at negligible cost
  (header-only peek, dominated by FASTQ I/O). The paired-end loop is now
  symmetrical with the single-end loop, which has always detected per file.
- Paired-end validation catches when R1 and R2 are the exact same filename
  (byte-equal path comparison; does not follow symlinks or canonicalise —
  matches v0.6.x behaviour).
- Paired-end invocations pre-flight-check for output-path collisions across
  pairs and abort before writing rather than silently overwriting. Comparison
  is case-insensitive (ASCII) so filenames differing only in letter-case are
  caught on APFS/NTFS (macOS/Windows default) as well as ext4 (Linux). Fixes
  issue #216.
  - On opt-in case-sensitive APFS/ZFS volumes, filenames differing only in
    letter-case are legitimately distinct; the pre-flight will still reject
    them. Use distinct base names or `--basename` per sample to work around
    this.
  - On partial failure at pair K, pairs 1..K&#8722;1 retain their complete
    outputs on disk and pair K may have a partial output file; pairs
    K+1..N are not attempted. Inspect and re-run only the failing pair —
    matches v0.6.x Perl behaviour (no rollback across pairs).


### Version 2.1.0 (Beta, Release on 18 Apr 2026)

**Major release — Rust rewrite (Oxidized Edition).** Faithful Rust rewrite of Trim Galore, delivered as a single binary with zero external dependencies and designed as a drop-in replacement for v0.6.x workflows. Outputs match v0.6.x for the core feature set; new capabilities beyond the Perl version include poly-G auto-detection and trimming, a generic poly-A trimmer, per-pair adapter auto-detection, cleaner multi-adapter invocation, a JSON MultiQC-native report, and other extensions. Built from `src/main.rs` (Cargo crate at repo root); the historical Perl script will be preserved at `legacy/trim_galore` once v2.1.0 GA ships (retained in the `0.6.11` tag during the beta window).

**Note on v2.0.0:** v2.0.0 was a pre-release cut inadvertently published to crates.io on 13 Apr 2026. It will be yanked when v2.1.0 GA ships. Users should install v2.1.0 or later.

#### Features (since v2.0.0)
- Multi-adapter support for both R1 and R2 via repeated `-a`/`-a2` flags and `file:adapters.fa` (c36b7fe)
- `--discard-untrimmed` flag (b0db3db)
- Multi-file single-end input (b0db3db)
- JSON trimming report for MultiQC native parsing (efedb95)
- MultiQC-compatible Cutadapt section in text trimming reports (121b821)

#### Bug fixes (since v2.0.0)
- Parallel-path `total_bp_after_trim` and `r2_clipped_5prime` stats now tracked correctly (82d1e34, 3996fc5)
- Cutadapt-section report values match v0.6.x MultiQC-parsed values (eedbc66)
- `--fastqc_args` accepts hyphenated values like `--nogroup` (def0344)
- Multi-member gzip FASTQ files decode correctly (9dcf519)
- Adapter auto-detection scans exactly 1M reads (9129650)

### Version 0.6.11 (Release on 24 Feb 2026)

- Added option `--rename` to write clipped bases to the read ID. Works in all modes with options `--clip_(r1/r2)` and `--three_prime_clip_(r1/r2)`, as well as `--hardtrim5` and `--hardtrim3`. Requested in [this issue](https://github.com/FelixKrueger/TrimGalore/issues/166).

- Added option `--bgiseq` to trim BGISEQ/DNBSEQ/MGISEQ adapters instead of the default auto-detection. Uses `AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA` for Read 1 (BGI/MGI forward), and
 `AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG` for Read 2 (BGI/MGI reverse). Requested in [issue#196](https://github.com/FelixKrueger/TrimGalore/issues/196)

- Added option  `--demux <barcode_file>` to demultiplex files from a 3'-end barcode after trimming is completed. Requested in [#199](https://github.com/FelixKrueger/TrimGalore/issues/199)

- Added option `--cutadapt_args "<ARGS>"` to pass extra arguments to Cutadapt, enabling use of advanced Cutadapt options without modifying Trim Galore.

- Changed `--clock` Epigenetic Clock processing behaviour for the 5' end of Read 2.

- Fixed `--demux` handling of CR (carriage return) characters in barcode files; fixed barcode length issue with NoCode; added barcode description to demux summary output.

- Fixed RRBS-specific trimming being silently bypassed when `--nextseq` and `--rrbs` are used together ([#210](https://github.com/FelixKrueger/TrimGalore/issues/210)).

### Version 0.6.10 (Release on 02 Feb 2023)

- Fixed a missing default value of `gzip` as the default decompression path (see [here](https://github.com/FelixKrueger/TrimGalore/commit/a3c6a64ae71657f1a282e01134293e424177a7d5)).

### Version 0.6.9 (Release on 29 Jan 2023)

- Fixed a declaration bug for `maxn_fraction` which had crept in during merging of different branches (see [here](https://github.com/FelixKrueger/TrimGalore/commit/cf9a9d97b723d3829dd902f1229d9c9b7cff8ba0)).

### Version 0.6.8 (Release on 28 Jan 2023)

- Added new option `--stranded_illumina` to allow trimming of the adapter sequence `ACTGTCTCTTATA` (whick looks like the Nextera sequence but with an additional A from A-tailing). See also here: https://github.com/FelixKrueger/TrimGalore/issues/127.

- Trim Galore will now preferentially use `igzip` for decompression, if installed. [More info here](https://github.com/FelixKrueger/TrimGalore/pull/149)

- finally dropped the option `--trim1` entirely. It wasn't useful beyond Bowtie 1 paire-end mode, and hence people should cease using it

- the option `--max_n COUNT` now interprets value between 0 and 1 as fraction of the read length (see [here](https://github.com/FelixKrueger/TrimGalore/issues/137))

- enabled the option `--max_length` also for paired-end trimming (of small RNAs)


### Version 0.6.7 (Release on 23 Jul 2021)

- just to get a DOI via Zenodo

### Version 0.6.6 (Release on 04 Sep 2020)

- Changed the way in which we test for the version of Cutadapt, more here: https://github.com/FelixKrueger/TrimGalore/issues/85

- Allowed specifying of multiple adapters for special cases. Works either via the command line, e.g.: `-a  " AGCTCCCG -a TTTCATTATAT -a TTTATTCGGATTTAT"` or via a FastA file, like so: `-a "file:multiple_adapters.fa"`  More info here: https://github.com/FelixKrueger/TrimGalore/issues/86.

- Added new special trimming mode for UMIs for the IMPLICON method ([`--implicon`](https://github.com/FelixKrueger/TrimGalore/issues/90)). In this mode, an 8bp UMI (unique molecular identifier) sequence is transferred from the start of Read 2 to the readID of both reads to allow UMI-aware deduplication (e.g. with `deduplicate_bismark --barcode` or [UmiBam](https://github.com/FelixKrueger/Umi-Grinder). Following this, Trim Galore will exit.

### Version 0.6.5 (Release on 19 Nov 2019)

- Added checks for whitespace(s) within input filenames, or a potential output folder name (supplied with `-o`). `[FATAL ERROR]` messages will advise users to use `_` instead.

- In a `--paired --basename BASE` scenario, the output files will now be called `BASE_val_1.fq.gz BASE_val_2.fq.gz` as described in the documentation (we previously also added `_R1` and `_R2`). This had to be addressed twice (0f631e5f979281fd4f18faef39818399a068a4b3 and 9ad019635a8a7f1aebb56f309889a7841a0ae42e) as the first approach was generating the Read 1 twice.

- removed a superflous warning statement for directional RRBS mode


### Version 0.6.4 (Release on 01 Aug 2019)

- Changed the adapter auto-detection procedure so that inconclusive detection always defaults to `--illumina`, unless none of the top 2, equal contaminants was 'Illumina', in which case it now defaults to `--nextera`. A warning message about this is now printed to the screen as well as to the trimming report.

- In addition to that, added the option `--consider_already_trimmed INT`. If no specific adapter exceeds this limit during the adapter auto-detection, the file is considered 'already adapter-trimmed' and will not be adapter trimmed again. Quality trimming is carried out as usual (technically, the adapter sequence is set to `-a X`). This option was added so that pipelines that are being fed either already-trimmed or untrimmed data will do the right thing in both cases.

- Changed the trimming mode for paired-end `--rrbs` in conjunction with `--non_directional`: previously, Read 2 was only trimmed for `CGA` or `CAA` at the 5' end, but not trimmed for read-through contamination at the 3' end if no 5' contamination had been removed. This problem had been introduced in v0.4.3, but since non-directional RRBS is not very common it had not been spotted so far. 

- File names for single-end trimming are now changed correctly when both `--output_dir` and `--basename` were specified together (was working correctly for PE mode already) 

### Version 0.6.3 (Release on 27 06 2019)

- Also added the number of PolyA trimmed bases to the start of the read in the format `trimmed_bases:A:`

So an example trimmed read would look like this:
```
 @READ-ID:1:1102:22039:36996 1:N:0:CCTAATCC
GCCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAATAAAAACTTTATAAACACCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

@32:A:READ-ID:1:1102:22039:36996_1:N:0:CCTAATCC_PolyA:32
GCCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAATAAAAACTTTATAAACACC

```


### Version 0.6.2 (Release on 08 05 2019)

- Changed the version checking mechanism so that Trim Galore even works in single-core mode if the version of Cutadapt was 7 years old...

- Fixed setting `-j $cores` for Cutadapt versions 2.X or above.


### Version 0.6.1 (Release on 20 Mar 2019)

* Added a check for very old versions Cutadapt, so that single-core trimming still works with Cutadapt versions prior to 1.15. 

- Fixed the way single-core trimming was dealt with in paired-end mode (which was introduced by the above 'fix')

- the option `--basename preferred_name` should now correctly work when specified in conjunction with `--output_dir`

### Version 0.6.0 (Release on 1 Mar 2019)

* Added option `--hardtrim3 INT,` which allows you to hard-clip sequences from their 5' end. This option processes one or more files (plain FastQ or gzip compressed files) and produces hard-trimmed FastQ files ending in .{INT}bp_3prime.fq(.gz). We found this quite useful in a number of scenarios where we wanted to remove biased residues from the start of sequences. Here is an example :

```
before:         CCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAAT
--hardtrim3 20:                                                   TTTTTAAGAAAATGGAAAAT
```

* Added new option `--basename <PREFERRED_NAME>` to use `PREFERRED_NAME` as the basename for output files, instead of deriving the filenames from the input files. Single-end data would be called `PREFERRED_NAME_trimmed.fq(.gz)`, or `PREFERRED_NAME_val_1.fq(.gz)` and `PREFERRED_NAME_val_2.fq(.gz)` for paired-end data. `--basename` only works when 1 file (single-end) or 2 files (paired-end) are specified, but not for longer lists (see #43).

* Added option `--2colour/--nextseq INT` whereby INT selects the quality cutoff that is normally set with `-q`, only that qualities of G bases are ignored. `-q` and `--2colour/--nextseq INT` are mutually exclusive (see [#41](https://github.com/FelixKrueger/TrimGalore/issues/41))

* Added check to see if Read 1 and Read 2 files were given as the very same file.

* If an output directory which was specified with `-o output_directory` did not exist, it will be created for you.

* The option `--max_n INT` now also works in single-end RRBS mode.

* Added multi-threading support with the new option `-j/--cores INT`; many thanks to Frankie James for initiating this. Multi-threading support works effectively if Cutadapt is run with Python 3, and if parallel gzip (`pigz`) is installed:

<img title="Multi-threading benchmark" style="float:right;margin:20px 20 20 600px" id="Multi-threading support" src="docs/Images/pigz_bench.png" >

For Cutadapt to work with multiple cores, it requires Python 3 as well as parallel gzip (pigz) installed on the system. The version of Python used is detected from the shebang line of the Cutadapt executable (either 'cutadapt', or a specified path). If Python 2 is detected, `--cores` is set to 1 and multi-core processing will be disabled. If `pigz` cannot be detected on your system, Trim Galore reverts to using `gzip` compression. Please note however, that `gzip` compression will slow down multi-core processes so much that it is hardly worthwhile, please see: [here](https://github.com/FelixKrueger/TrimGalore/issues/16#issuecomment-458557103) for more info).

Actual core usage: It should be mentioned that the actual number of cores used is a little convoluted. Assuming that Python 3 is used and `pigz` is installed, `--cores 2` would use:

- 2 cores to read the input (probably not at a high usage though)
- 2 cores to write to the output (at moderately high usage)
- 2 cores for Cutadapt itself
- 2 additional cores for Cutadapt (not sure what they are used for)
- 1 core for Trim Galore itself

So this can be up to 9 cores, even though most of them won't be used at 100% for most of the time. Paired-end processing uses twice as many cores for the validation (= writing out) step as Trim Galore reads and writes from and to two files at the same time, respectively.

`--cores 4` would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) +     1 (Trim Galore) = ~15 cores in total.

From the graph above it seems that `--cores 4` could be a sweet spot, anything above appear to have diminishing returns.


### 28-06-18: Version 0.5.0 

* Adapters can now be specified as single bases with a multiplier in squiggly brackets, e.g. -a "A{10}" to trim poly-A tails

* Added option `--hardtrim5 INT` to enable hard-clipping from the 5' end. This option processes one or more files (plain FastQ or gzip compressed files) and produce hard-trimmed FastQ files ending in `.{INT}bp.fq(.gz)`. 

* Added option `--clock` to trim reads in a specific way that is currently used for the Mouse Epigenetic Clock (see here: [Multi-tissue DNA methylation age predictor in mouse, Stubbs et al., Genome Biology, 2017 18:68](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1203-5)). Following the trimming, Trim Galore exits.

In it's current implementation, the dual-UMI RRBS reads come in the following format:

```
Read 1  5' UUUUUUUU CAGTA FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF TACTG UUUUUUUU 3'
Read 2  3' UUUUUUUU GTCAT FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF ATGAC UUUUUUUU 5'
```

Where UUUUUUUU is a random 8-mer unique molecular identifier (UMI), CAGTA is a constant region,
and FFFFFFF... is the actual RRBS-Fragment to be sequenced. The UMIs for Read 1 (R1) and
Read 2 (R2), as well as the fixed sequences (F1 or F2), are written into the read ID and
removed from the actual sequence. Here is an example:

```
R1: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 1:N:0: CGATGTTT
    ATCTAGTTCAGTACGGTGTTTTCGAATTAGAAAAATATGTATAGAGGAAATAGATATAAAGGCGTATTCGTTATTG
R2: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 3:N:0: CGATGTTT
    CAATTTTGCAGTACAAAAATAATACCTCCTCTATTTATCCAAAATCACAAAAAACCACCCACTTAACTTTCCCTAA

R1: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 1:N:0: CGATGTTT:R1:ATCTAGTT:R2:CAATTTTG:F1:CAGT:F2:CAGT
                 CGGTGTTTTCGAATTAGAAAAATATGTATAGAGGAAATAGATATAAAGGCGTATTCGTTATTG
R2: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 3:N:0: CGATGTTT:R1:ATCTAGTT:R2:CAATTTTG:F1:CAGT:F2:CAGT
                 CAAAAATAATACCTCCTCTATTTATCCAAAATCACAAAAAACCACCCACTTAACTTTCCCTAA
```
Following clock trimming, the resulting files (.clock_UMI.R1.fq(.gz) and .clock_UMI.R2.fq(.gz))
should be adapter- and quality trimmed with Trim Galore as usual. In addition, reads need to be trimmed
by 15bp from their 3' end to get rid of potential UMI and fixed sequences. The command is:

`trim_galore --paired --three_prime_clip_R1 15 --three_prime_clip_R2 15 *.clock_UMI.R1.fq.gz *.clock_UMI.R2.fq.gz`

Following this, reads should be aligned with Bismark and deduplicated with UmiBam in `--dual_index` mode (see here: https://github.com/FelixKrueger/Umi-Grinder). UmiBam recognises the UMIs within this pattern: R1:(**ATCTAGTT**):R2:(**CAATTTTG**): as (UMI R1=**ATCTAGTT**) and (UMI R2=**CAATTTTG**).


### 13-11-17: Version 0.4.5 released

* Trim Galore now dies during the validation step when it encounters paired-end files that are not equal in length


### 28-03-17: Version 0.4.4 released

* Reinstated functionality of option `--rrbs` for single-end RRBS reads which had gone amiss in the previous release. What happened in detail was that RRBS trimming was de facto skipped if there was only a single file specified.

* Updated User Guide and Readme documents, added Installation instruction and Travis functionality - thanks Phil!


### 25-01-17: Version 0.4.3 released

* Changed the option `--rrbs` for paired-end libraries from removing 2 additional base pairs from the 3' end of both reads to trim 2 bp from the 3' end only for Read 1 and set `--clip_r2 2` for Read 2 instead. This is because Read 2 does not technically need 3' trimming since the end of Read 2 is not affected by the artificial methylation states introduced by the [end-repair] fill-in reaction. Instead, the first couple of positions of Read 2 suffer from the same fill-in problems as standard [paired-end libraries](https://sequencing.qcfail.com/articles/library-end-repair-reaction-introduces-methylation-biases-in-paired-end-pe-bisulfite-seq-applications/). Also see [this issue](https://github.com/FelixKrueger/TrimGalore/issues/3).

* Added a closing statement for the REPORT filehandle since it occasionally swallowed the last line...

* Setting `--length` now takes priority over the smallRNA adapter (which would set the length cutoff to 18 bp).


### 07-09-16: Version 0.4.2 released
* Replaced `zcat` with `gunzip -c` so that older versions of Mac OSX do not append a .Z to the end of the file and subsequently fail because the file is not present. Dah...
*	Added option `--max_n COUNT` to remove all reads (or read pairs) exceeding this limit of tolerated Ns. In a paired-end setting it is sufficient if one read exceeds this limit. Reads (or read pairs) are removed altogether and are not further trimmed or written to the unpaired output
* Enabled option `--trim-n` to remove Ns from both end of the reads. Does currently not work for RRBS-mode
* Added new option `--max_length INT` which removes reads that are longer than INT bp after trimming. This is only advised for smallRNA sequencing to remove non-small RNA sequences

### 12-11-15: Version 0.4.1 released: Essential update for smallRNA libraries!
*	Changed the Illumina small RNA sequence used for auto-detection to `TGGAATTCTCGG` (from formerly `ATGGAATTCTCG`). The reason for this is that smallRNA libraries have ssRNA adapters ligated to their -OH end, a signature of dicer cleavage, so there is no A-tailing involved. Thanks to Q. Gouil for bringing this to our attention
*	Changed the length cutoff for sequences to 16bp (down from 20bp) for smallRNA libraries before sequences get removed entirely. This is because some 20-23bp long smallRNAs species that had removed T, TG, or TGG etc. might just about pass the 20bp cutoff
*	Added a small description to the `--help` message for users of the NuGEN Ovation RRBS kit as to *NOT* use the `--rrbs` option (see `--help`)

### 06-05-15: Version 0.4.0 released
*	Unless instructed otherwise Trim Galore will now attempt to auto-detect the adapter which had been used for library construction (choosing from the Illumina universal, Nextera transposase and Illumina small RNA adapters). For this the first 1 million sequences of the first file specified are analysed. If no adapter can be detected within the first 1 million sequences Trim Galore defaults to --illumina. The auto-detection behaviour can be overruled by specifying an adapter sequence or using `--illumina`, `--nextera` or `--small_rna`
*	Added the new options `--illumina`, `--nextera` and `--small_rna` to use different default sequences for trimming (instead of `-a`): Illumina: `AGATCGGAAGAGC`; Small RNA: `TGGAATTCTCGG`; Nextera: `CTGTCTCTTATA`
*	Added a sanity check to the start of a Trim Galore run to see if the (first) FastQ file in question does contain information at all or appears to be in SOLiD colorspace format, and bails if either is true. Trim Galore does not support colorspace trimming, but users wishing to do this are kindly referred to using Cutadapt as a standalone program
*	Added a new option `--path_to_cutadapt /path/to/cudapt`. Unless this option is specified it is assumed that Cutadapt is in the PATH (equivalent to `--path_to_cutadapt cutadapt`). Also added a test to see if Cutadapt seems to be working before the actual trimming is launched
*	Fixed an open command for a certain type of RRBS processing (was open() instead of open3())

### 16-07-14: Version 0.3.7 released
*	Applied small change that makes paired-end mode work again (it was accidentally broken by changing @ARGV for @filenames when looping through the filenames...)

### 11-07-14: Version 0.3.6 released
*	Added the new options `--three_prime_clip_r1` and `--three_prime_clip_r2` to clip any number of bases from the 3' end after adapter/quality trimming has completed
* Added a check to see if Cutadapt exits fine. Else, Trim Galore will bail a well
*	The option `--stringency` needs to be spelled out now since using `-s` was ambiguous because of `--suppress_warn`

### late 2013: Version 0.3.5 released
*	Added the Trim Galore version number to the summary report

### 19-09-13: Version 0.3.4 released
*	Added single-end or paired-end mode to the summary report
*	In paired-end mode, the Read 1 summary report will no longer state that no sequence have been discarded due to trimming. This will be stated in the trimming report of Read 2 once the validation step has been completed

### 10-09-13: Version 0.3.3 released
*	Fixed a bug what was accidentally introduced which would add an additional empty line in single-end trimming mode

### 03-09-13: Version 0.3.2 released
*	Specifying `--clip_R1` or `--clip_R2` will no longer attempt to clip sequences that have been adapter- or quality-trimmed below the clipping threshold
*	Specifying an output directory with `--rrbs` mode should now correctly create temporary files

### 15-07-13: Version 0.3.1 released
*	The default length cutoff is now set at an earlier timepoint to avoid a clash in paired-end mode when `--retain_unpaired` and individual read lengths for read 1 and read 2 had been defined

### 15-07-13: Version 0.3.0 released
*	Added the options `--clip_R1` and `--clip_R2` to trim off a fixed amount of bases at from the 5' end of reads. This can be useful if the quality is unusually low at the start, or whenever there is an undesired bias at the start of reads. An example for this could be PBAT-Seq in general, or the start of read 2 for every bisulfite-Seq paired-end library where end repair procedure introduces unmethylated cytosines. For more information on this see the M-bias section of the Bismark User Guide.

### 10-04-13: Version 0.2.8 released
* Trim Galore will now compress output files with GZIP on the fly instead of compressing the trimmed file once trimming has completed. In the interest of time temporary files are not being compressed
*	Added a small sanity check to exit if no files were supplied for trimming. Thanks to P. for 'bringing this to my attention'

### 01-03-13: Version 0.2.7 released
*	Added a new option `--dont_gzip` that will force the output files not to be gzip compressed. This overrides both the `--gzip` option or a .gz line ending of the input file(s)

### 07-02-13: Version 0.2.6 released
*	Fixes some bugs which would not gzip or run FastQC correctly when the option `-o` had been specified
*	When `--fastqc` is specified in paired-end mode the intermediate files '_trimmed.fq' are no longer analysed (only files '_val_1' and '_val_2')

### 19-10-12: Version 0.2.5 released
*	Added option `-o/--output_directory` to redirect all output (including temporary files) to another folder (required for implementation into Galaxy)
*	Added option `--no_report_file` to suppress a report file
*	Added option `--suppress_warn` to suppress any output to STDOUT or STDERR

### 02-10-12: Version 0.2.4 released
*	Removed the shorthand `-l` from the description as it might conflict with the paired-end options `-r1/--length1` or `-r2/--length2`. Please use `--length` instead
*	Changed the reporting to show the true Phred score quality cutoff
*	Corrected a typo in stringency...

### 31-07-12: Version 0.2.3 released
* Added an option `-e ERROR RATE` that allows one to specify the maximum error rate for trimming manually (the default is 0.1)

### 09-05-12: Version 0.2.2 released
*	Added an option `-a2/--adapter2` so that one can specify individual adapter sequences for the two reads of paired-end files; hereby the sequence provided as `-a/--adapter` is used to trim read 1, and the sequence provided as `-a2/--adapter2` is used to trim read 2. This option requires `--paired` to be specified as well

### 20-04-12: Version 0.2.1 released
*	Trim Galore! now has an option `--paired` which has the same functionality as the validate_paired_ends script we offered previously. This option discards read pairs if one (or both) reads of a read pair became shorter than a given length cutoff
*	Reads of a read-pair that are longer than a given threshold but for which the partner read has become too short can optionally be written out to single-end files. This ensures that the information of a read pair is not lost entirely if only one read is of good quality
*	Paired-end reads may be truncated by a further 1 bp from their 3' end to avoid problems with invalid alignments with Bowtie 1 (which regards alignments that contain each other as invalid...)
*	The output may be gzip compressed (this happens automatically if the input files were gzipped (i.e. end in .gz))
*	The documentation was extended substantially. We also added some recommendations for RRBS libraries for MseI digested material (recognition motif TTAA)

### 21-03-12: Version 0.1.4 released
* Phred33 (Sanger) encoding is now the default quality scheme
*	Fixed a bug for Phred64 encoding that would occur if several files were specified at once

### 14-03-12: Version 0.1.3 released
*	Initial stand-alone release; all basic functions working
*	Added the option `--fastqc_args` to pass extra options to FastQC

