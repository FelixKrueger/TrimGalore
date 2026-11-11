# Code review — Reviewer A — `1b16bed` + `8d5567a`

**Target:** `/Users/fkrueger/Github/TrimGalore`, branch `chore/fastqc-rust-1.1.1`, commits
`1b16bed` (docs sweep) and `8d5567a` (`--phred64` wiring), diffed against `d5f41c5`.
Reviewed in place; no worktree created, no branch switched.

**Verdict:** the wiring is correct at every call site. Two gaps worth acting on before release:
a new hard-failure path with an unhelpful diagnostic (§1), and a test that does not constrain
most of what the commit changed (§2, §3). The doc sweep is accurate but incomplete (§4).

---

## Build verification

All run with `DEVELOPER_DIR=/Library/Developer/CommandLineTools` and
`CARGO_TARGET_DIR=<scratchpad>/target_revA`.

| Check | Result |
|---|---|
| `cargo clippy --all-targets --release -- -D warnings`, **cold cache** | exit 0, 0 warnings. 46 crates compiled, 30.3 s — a real compile, not a cache hit |
| `cargo test --release` | **826 passed, 0 failed, 0 ignored.** Both new tests pass |
| `cargo fmt --all -- --check` | fails — but only on **untracked** `examples/fastqc_only.rs` (§6). CI unaffected |
| `--version` / `-V` provenance | long form keeps `8d5567a — macos/aarch64 — built …` (em-dash); `-V` omits it. Intact |

---

## Q1 — Call-site audit: all correct, but there are **12**, not 11

I classified each site from its own `output_path` construction rather than from the commit.
`src/main.rs:910` returns `run_ubam_output(...)` before `run_paired`/`run_single_file` are
reachable (`src/main.rs:1045`, `:1094`), so the FASTQ and BAM drivers are cleanly disjoint.

**FASTQ output — 7 calls, all pass `phred64`:**

| Site | Output path origin | Arg |
|---|---|---|
| `src/main.rs:1566` `run_single_file` | FASTQ trimmer output | `cli.phred64` ✓ |
| `src/main.rs:1884` `run_paired` R1 | `outputs.val_1` (`src/main.rs:1933`) | `cli.phred64` ✓ |
| `src/main.rs:1891` `run_paired` R2 | `outputs.val_2` (`src/main.rs:1934`) | `cli.phred64` ✓ |
| `src/main.rs:1904` `--passthrough` carrier | `passthrough_output` (`src/main.rs:1619`) | `cli.phred64` ✓ |
| `src/clump_only.rs:421` `clump_only_single` | `naming::clumped_output_name` (`:323`) | `phred64` ✓ |
| `src/clump_only.rs:612` `clump_only_paired` R1 | `:472` | `phred64` ✓ |
| `src/clump_only.rs:613` `clump_only_paired` R2 | `:472` | `phred64` ✓ |

**BAM output — 5 calls, all pass `false`:**

| Site | Output path origin | Arg |
|---|---|---|
| `src/main.rs:2364` `run_ubam_output_single` | BAM | `false` ✓ |
| `src/main.rs:2523` `run_ubam_output_paired_two_files` | BAM | `false` ✓ |
| `src/main.rs:2643` `run_ubam_output_paired_single_file` | BAM | `false` ✓ |
| `src/clump_only.rs:917` `clump_only_single_to_bam` | `naming::clumped_bam_output_name` (`:823`) | `false` ✓ |
| `src/clump_only.rs:1177` `clump_only_paired_to_bam_one_pair` | `naming::clumped_paired_bam_output_name` (`:1003`, `:1049`) | `false` ✓ |

No site was missed and none is misclassified. `src/fastqc.rs:62` is the **only** `FastQCConfig`
construction in production code (`examples/fastqc_only.rs` is the only other, and is untracked),
so there is no second entry point that could bypass the parameter.

The brief said 11 sites. There are 12 — most likely `run_paired`'s R1/R2 counted as one, or the
`--passthrough` call was not counted. Flagging it because a miscount is exactly how a site gets
skipped; in this case nothing was.

The premise holds independently: `cli.phred_offset()` (`src/cli.rs:1223`) reaches only
`quality.rs` comparisons and `BamWriter` — `grep` over `src/` shows no quality-string rewrite on
the FASTQ writer — and `docs/src/content/docs/guide/quality.md:33` states the same contract
("FASTQ-to-FASTQ runs preserve it, so Phred+64 in gives Phred+64 out").

## Q2 — Reachability: every wired site is live; no dead code

`--phred64` is rejected only for BAM **input** (`src/main.rs:343`). `Cli::validate` never rejects
it otherwise — `grep -n phred64 src/cli.rs` returns only the clap declaration (`:155`) and
`phred_offset()` (`:1223`). So:

- All 7 FASTQ sites are reachable with `--phred64`, including `--clump_only` (lossless reorder
  preserves the quality bytes verbatim, so declaring Phred+64 is right) and `--passthrough`.
- The 5 BAM sites are reachable with `--phred64` for **FASTQ input → uBAM output**
  (`phred64_ubam_out_se_stores_true_phred` et al. exercise exactly that), so hardcoding `false`
  there is a live decision, not dead code. It is the correct one — `known_phred_encoding()`
  returns `Some(SANGER)` for BAM and short-circuits the `.or()` at
  `fastqc-rust-1.1.1/src/runner.rs:396-398`.

---

## Findings

### 1. [High] `--phred64 --fastqc` now aborts the run after writing all output, with a diagnostic that never names `--phred64`

**Verified empirically, before and after.** Built the branch binary, ran it on a 2-read Phred+33
fixture with `--phred64 --fastqc -q 0 --length 0`:

```
Running FastQC on …/p33_trimmed.fq
Failed to process p33_trimmed.fq: Phred64 encoding is incompatible with having ASCII char '#' with value 35) in the file
Error: FastQC analysis failed (exit code 1)
EXIT=1
```

Output directory afterwards: `p33_trimmed.fq`, `p33.fastq_trimming_report.txt`,
`p33.fastq_trimming_report.json` — **and no `_fastqc.html` / `_fastqc.zip`**.

I then forced `phred64: false` at `src/fastqc.rs:65` to reproduce the pre-commit state and reran
the identical command: **exit 0**, report produced. So `8d5567a` introduces this failure path.

The mechanism is upstream and deliberate: `phred::detect(lowest_char, true)` returns `Err` when
any quality char is below ASCII 64 (`fastqc-rust-1.1.1/src/utils/phred.rs:66-72`), and
`BasicStats::write_text_report` converts that to an `io::Error` with the comment *"Java fails the
whole file when the encoding can't be resolved"* (`src/modules/basic_stats.rs:394-397`).

**The new behaviour is correct** — a hard failure beats silently inflating every mean by 31, and
it matches Java FastQC. Three things around it are not:

1. **The diagnostic never mentions `--phred64`.** The user sees upstream's string (complete with
   an unbalanced `)` — upstream's typo, not ours) and has to work backwards to their own flag.
2. **Trim Galore already has the right message, scoped to the wrong path.** The `NOTE:` at
   `src/main.rs:506` explains this exact mistake in detail, but it sits inside
   `if matches!(cli.output_format, OutputFormat::UBam)` (`src/main.rs:487`) — so it fires on the
   path that *cannot* fail this way and stays silent on the path that now can. I confirmed it does
   not appear in the failing run's output (`grep -c "NOTE: input is Phred+64"` → 0).
3. **Exit 1 with trimmed output and both reports on disk.** A pipeline gating on exit status now
   treats a successful trim as a failed step because an optional QC add-on rejected an encoding
   declaration. This repo cares about exactly this shape (`tests/integration_partial_output.rs`,
   `plans/08132026_partial-output-on-refusal/`).

Realistic triggers: user error (the common case), and genuine Solexa / Illumina 1.0 data, whose
quality alphabet legitimately starts at ASCII 59.

**Recommendation (low-risk, no behaviour change):** hoist the `--phred64` NOTE out of the
`UBam`-only block at `src/main.rs:487` so it fires whenever `--phred64` is given, and add a
sentence to the `:::caution` block at `docs/src/content/docs/guide/quality.md:39-41`, which
currently says only that the NOTE appears "when `--phred64` is combined with uBAM output". If you
want the trim to survive the QC rejection, that is a larger call and belongs in its own issue —
I am not recommending it here.

### 2. [Medium] The new test covers one of seven FASTQ call sites; a PE regression survives the whole suite

`tests/integration_fastqc_phred64.rs` exercises single-end FASTQ output only. I mutated **both**
`run_paired` sites (`src/main.rs:1889`, `:1896`) from `cli.phred64` to `false` and ran the full
suite:

```
TOTAL passed=826 failed=0
```

The mutant survives everything. The same holds by inspection for `src/main.rs:1904`
(`--passthrough`) and `src/clump_only.rs:421/612/613` — no test combines `--clump_only` or
`--passthrough` with `--phred64 --fastqc` and then reads inside the zip.

PE is the common case for this tool, and §1 means a PE regression is no longer silent (it would
now hard-fail on genuinely Phred+64 PE data only if the data dipped below ASCII 64 — otherwise it
is silently 31 too high again, exactly the bug being fixed).

**Recommendation:** add a third test using `test_files/phred64_test_R{1,2}.fastq` (already tracked
— same 24×`h` / 10×`B` shape, confirmed) asserting `Encoding` and the per-base means on both
`_val_1` and `_val_2` reports. That one test closes three of the four uncovered FASTQ sites'
failure modes by sharing the `cli.phred64` plumbing; a `--clump_only` case would close the rest.

### 3. [Medium] `ubam_output_is_not_told_about_phred64` does not test what its name claims

I mutated all three uBAM sites (`src/main.rs:2364`, `:2523`, `:2643`) from `false` to
`cli.phred64` and ran `integration_fastqc_phred64`, `integration_ubam_out` and
`integration_clump_only_ubam`:

```
test phred64_reaches_the_fastqc_config_on_fastq_output ... ok
test ubam_output_is_not_told_about_phred64 ... ok
test result: ok. 2 passed; 0 failed
test result: ok. 17 passed; 0 failed
test result: ok. 20 passed; 0 failed
```

The mutant survives. It has to: `known_phred_encoding()` returns `Some(SANGER)` for BAM and
short-circuits `.or()`, so the report reads `Sanger / Illumina 1.9` whether we pass `true` or
`false`. The test asserts the *upstream precedence contract*, not the argument.

The in-test comment is honest about this ("this pins the contract the BAM call sites rely on").
The **name** is not: it reads as a guard on the call sites, and a future maintainer trusting it
will believe the BAM arguments are pinned when they are not. Strictly, the test fails only on the
conjunction (upstream flips precedence **and** we pass `true`) — it cannot detect either alone.

**Recommendation:** rename to something like `bam_format_encoding_wins_over_the_phred64_flag` and
drop one line into the module doc saying the BAM sites' `false` is a local statement of intent
that no test constrains. No code change needed.

Other mutants I considered: hardcoding `phred64: true` in `src/fastqc.rs` **is** caught —
indirectly, by every existing `--fastqc` test on a Phred+33 fixture, since §1's hard failure would
make them exit non-zero (`tests/integration_preflight_tripwire.rs:508`, `:759` and four others).
Worth knowing that this is incidental coverage, not designed.

### 4. [Medium] The doc sweep missed `test_files/README.md`, which now describes behaviour 0.13.0 removed

`test_files/README.md:14` (tracked) says of the `phred64_test.fastq` fixture:

> a uniform-quality fixture … would additionally trip FastQC's encoding heuristic (upstream
> FastQC 0.12.1 behaviour, reproducible on either output format: it keys off the lowest observed
> quality character, and an all-Q40 record leaves nothing below ASCII 64, so the guess falls to
> Phred+64 and every mean comes out 31 too low)

Under 1.1.1 that no longer happens. `detect(lowest_char, false)` now returns `SANGER`
unconditionally (`fastqc-rust-1.1.1/src/utils/phred.rs:64-66`), and the upstream unit test
comments the change in so many words: `// Formerly misdetected as Illumina 1.5`
(`src/utils/phred.rs:103`). 0.13.0 emits a *warning* via `phred64_suspicion` instead.

This matters more than a stale sentence usually would: it is the note a maintainer reads when
deciding whether the fixture's two-value design is still needed, and it now argues from a
mechanism that is gone. The surviving `0.12.1` references in `CHANGELOG.md:1398,1420` and
`docs/src/content/docs/reference/changelog.md:654,676` are historical release entries and are
correctly left alone.

**Recommendation:** one sentence — the 0.12.1 misdetection is gone in 0.13.0; the two-value
fixture is retained because it exercises both offsets.

### 5. [Medium] Nothing in the CHANGELOG covers any of this

`CHANGELOG.md` has an **empty** `### Unreleased` section; `### Version 2.4.0 (Release on 3
November 2026)` is finalized and `Cargo.toml:3` is still `2.4.0`. So three user-visible changes
have no entry: the FastQC 0.13.0 retarget (report gains Total Bases / Mean Length / Median Length,
drops "Sequences flagged as poor quality"), `--svg` becoming a no-op with `--png` added, and the
`--phred64` fix itself.

This may well be deliberate while the version number is undecided — flagging it so the decision is
explicit rather than forgotten. Per `CLAUDE.md`, one sentence each.

### 6. [Low] `examples/fastqc_only.rs` fails `cargo fmt --all -- --check`

```
Diff in /Users/fkrueger/Github/TrimGalore/examples/fastqc_only.rs:8:
-    trim_galore::fastqc::run(Path::new(&args[1]), None, Some(Path::new(&args[2])), 1, false)?;
+    trim_galore::fastqc::run(
+        Path::new(&args[1]), …
```

The file is **untracked** (`git ls-files examples/` lists only `mk_bam_fixture.rs` and
`mk_clump_fixture.rs`), so CI's `cargo fmt --all -- --check` on a clean checkout passes and the
shipped crate is unaffected. It was updated with the new 5th argument to keep the local
`--all-targets` build compiling, which is why it showed up in my cold-cache clippy run. Either
`cargo fmt` it or leave it — just be aware `cargo fmt --all -- --check` is red in this working
tree for a reason unrelated to the commits.

### 7. [Low] Two adjacent same-typed bools in the `clump_only` signatures

`phred64: bool` was inserted immediately before `no_report_file: bool` in both
`clump_only_single` (`src/clump_only.rs:295-296`) and `clump_only_paired` (`:444-445`). A
transposition at a call site compiles silently.

I checked all 14 callers. The two in `src/main.rs` (`:700`, `:736`) are in signature order. The 12
unit-test callers all pass `false, None, false, false`, and the new `false` was appended at the end
of each — so even if the author had meant to insert it one position earlier, both values are
`false` and nothing shifts. No test anywhere passes `no_report_file: true` through these
functions. Correct today; no type safety for the next caller. Mentioning it only because the
functions already carry `#[allow(clippy::too_many_arguments)]`, which is the signal that a
parameter struct would serve better — not worth doing for this commit.

### 8. [Low] Duplicated archive-reading in the new test helpers

`basic_stat` and `per_base_mean` (`tests/integration_fastqc_phred64.rs:19-37`, `:40-62`) each open
the zip and read `fastqc_data.txt` in full; each test opens the same archive three times. They also
disagree on error detail — `unwrap_or_else(|e| panic!("opening {}: {e}", …))` vs
`.expect("report archive")`. A single `fn report_data(zip, stem) -> String` would halve the file.
Cosmetic.

---

## Q4 — Doc claims: verified against the binary and the crate, all accurate

I did not take the commit message's word for any of these.

- **"FastQC 0.13.0-compatible"** — `fastqc-rust-1.1.1/src/lib.rs:10` is
  `pub const VERSION: &str = "0.13.0"`; its README line 22 reads "Currently tracking upstream Java
  FastQC version `0.13.0`". The built binary emits the claim in `--help`. Consistent with the
  crate's own statement.
- **`--svg` is a no-op; `--png` selects PNG** — verified by running both and inspecting the HTML.
  Default/`--svg`: **8** `data:image/svg+xml` charts, **0** PNG charts (the 21 `data:image/png`
  hits are the tick/warning/error icons and the logo). `--png`: **0** SVG, **29** PNG — the same
  21 icons plus the 8 charts converted. Upstream's own CLI agrees:
  `/// No-op: SVG is now the default. Kept so existing commands still work.`
  (`fastqc-rust-1.1.1/src/main.rs:140-142`) and `/// Embed charts in the HTML report as PNG
  instead of SVG.` (`:136-138`). Neither flag emits the "not translated" warning (count 0 for
  both). The `svg_output` field is genuinely gone from `FastQCConfig` in 1.1.1; `png_output` is
  new.
- **Tripwire comment** (`tests/integration_preflight_tripwire.rs:31-32`), "neither `--png` nor the
  no-op `--svg` writes an image file beside the primary artefacts" — verified: the `--png` run's
  output directory holds only `_fastqc.html`, `_fastqc.zip`, the trimmed `.fq` and the two
  trimming reports. (Both `Images/*.svg` and `Images/*.png` are inside the zip either way, which
  the comment's "beside" wording correctly does not claim otherwise.)

**Soundness of the new test's fixture reasoning** (Q3, mechanical part): `h` = ASCII 104 → Q40,
`B` = ASCII 66 → Q2 under Phred+64; lowest char 66 ≠ 65 so `detect` names it `Illumina 1.5`,
matching the assertion. Read misread as Phred+33 would give `71.0` / `33.0`, so the per-base
assertions carry real signal, not just the `Encoding` row. The 34 bp read is short enough that
FastQC does not group bases, so a row labelled exactly `"25"` exists — and the adapter does not
trim (read ends `…GTAC`, Illumina adapter starts `A`), so base 25 survives. All confirmed by the
test passing.

## Q5 — `zip` dev-dependency: cannot affect the shipped binary, and the tree really is unchanged

The spec at `Cargo.toml:56` is byte-identical to fastqc-rust's own
(`fastqc-rust-1.1.1/Cargo.toml:124-130`: `version = "9"`, `default-features = false`,
`features = ["deflate", "time"]`), so the comment's claim holds. The `Cargo.lock` diff is exactly
one line — `+ "zip",` appended to `trim-galore`'s dependency list — with **no** new `[[package]]`
block and no version change, which is the direct evidence that resolution did not move. `zip`
stays at 9.3.0.

Dev-dependencies are excluded from `cargo build`, so the release binary cannot change regardless;
feature unification could in principle widen the `cargo test` build only, and here it does not
because the feature sets are identical. Reproducibility is untouched (no new build scripts, no
time/path inputs).

Worth recording while looking: `flate2` in the lock now resolves to
`crc32fast + miniz_oxide + zlib-rs` with **no** `libz-sys` (`Cargo.lock:529-537`), so the
C-zlib leak via fastqc-rust's old `default = ["native-zlib"]` is gone. The memory note
`reference_fastqc_rust_zlib.md` is now stale.

---

## What I could not verify

- **Byte-identity of the report against Java FastQC 0.13.0.** No Java FastQC on this machine. The
  "0.13.0-compatible" claim is verified only as far as the crate's own `VERSION` and README
  assert it. `CLAUDE.md:67` is now honest that the CI smoke test checks file presence, not content
  — so nothing in this repo verifies it either, and that is a conscious, documented position.
- **The CI `validation` job.** It installs Perl Trim Galore 0.6.11 via conda from a pinned
  upstream URL; not reproducible here. The FastQC smoke test asserts `fastqc` is absent from
  `$PATH` and that the zip contains `summary.txt`, `fastqc_data.txt`, `Images/`, `Icons/` — all
  four present in the archives I produced, so that assertion should hold.
- **Downstream MultiQC parsing of 0.13.0's new Basic Statistics rows** (Total Bases, Mean Length,
  Median Length) and the dropped "Sequences flagged as poor quality" row. Out of scope for these
  two commits, but it is the kind of thing the retarget changes and no test here covers it.

---

## Priority summary

| # | Priority | Finding |
|---|---|---|
| 1 | **High** | `--phred64 --fastqc` aborts after writing output; diagnostic never names the flag; the explanatory NOTE is scoped to the uBAM-only block (`src/main.rs:487`) |
| 2 | Medium | PE / passthrough / clump_only FASTQ sites untested — PE mutant survives all 826 tests |
| 3 | Medium | `ubam_output_is_not_told_about_phred64` tests the upstream contract, not the call sites; mutant survives |
| 4 | Medium | `test_files/README.md:14` describes 0.12.1 misdetection that 0.13.0 removed |
| 5 | Medium | No CHANGELOG entry for the retarget, the `--svg`/`--png` change, or the fix |
| 6 | Low | Untracked `examples/fastqc_only.rs` is unformatted (CI unaffected) |
| 7 | Low | Two adjacent `bool` params in the `clump_only` signatures (all callers verified correct) |
| 8 | Low | Duplicated zip-reading in the new test's helpers |

No fixes applied — findings 1, 2, 3 and 5 all involve judgement calls that are yours (where the
NOTE belongs, how much test surface to buy, which version these commits ship under). Finding 4 is
a one-sentence edit I would have made directly if it were not in a file the commits deliberately
left alone.
