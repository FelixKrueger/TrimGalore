# Code review — Reviewer B — `1b16bed` + `8d5567a`

**Scope:** `git diff d5f41c5..8d5567a` on `chore/fastqc-rust-1.1.1`, reviewed in place (no worktree).
**Date:** 2026-11-05
**Verdict:** The fix is correct and does what the commit says on the path it targets. One real
behaviour change is undocumented and makes an existing source comment false; the test coverage is
narrower than the commit message implies; two pieces of prose (one of them in the new test file) state
the upstream mechanism backwards.

Established context taken as given: the `=1.0.1` → `=1.1.1` pin move in `d5f41c5`, and the decision to
retarget at Java FastQC 0.13.0.

---

## 1. What actually changed for a user

Built both commits (release, separate `CARGO_TARGET_DIR`s, `d5f41c5` from a `git archive` pin so the
working tree was never switched) and compared real report content out of `fastqc_data.txt`.

| Invocation (fixture) | `d5f41c5` | `8d5567a` |
|---|---|---|
| `--phred64 --fastqc -q 0 --length 0` (p64) | `Sanger / Illumina 1.9`, base1 **71.0**, base25 **33.0** | `Illumina 1.5`, base1 **40.0**, base25 **2.0** |
| `--phred64 --fastqc` **default `-q 20`** (p64) | `Sanger / Illumina 1.9`, base1 **71.0** | `Illumina 1.5`, base1 **40.0** |
| `--paired --phred64 --fastqc` (p64 ×2) | both mates 71.0 / 33.0 | both mates 40.0 / 2.0 |
| `--clump_only --phred64 --fastqc` (p64) | `Sanger / Illumina 1.9`, 71.0 / 33.0 | `Illumina 1.5`, 40.0 / 2.0 |
| `--output-format ubam --phred64 --fastqc` (p64) | `Sanger / Illumina 1.9`, 40.0 / 2.0 | **identical** |
| `--clump_only --output-format ubam --phred64 --fastqc` | — | `Sanger / Illumina 1.9`, 40.0 / 2.0 |
| gzipped input, `--cores 4`, `--phred64 --fastqc` | — | `Illumina 1.5`, 40.0 / 2.0 |
| `--fastqc_args "-t 4"` / `"-o OTHERDIR"` + `--phred64` | — | declaration survives both; `Illumina 1.5`, 40.0 |

So: the FASTQ paths are fixed, the uBAM paths are bit-for-bit unchanged, and the mates are told the
same thing. The diagnosis in the commit message is confirmed at source —
`fastqc-rust-1.1.1/src/utils/phred.rs:56-65` returns `PhredEncoding::SANGER` unconditionally when
`phred64` is false, and `src/runner.rs:396-398` resolves
`seq_file.known_phred_encoding().or(config.phred64.then_some(PHRED64))`, with
`src/sequence/bam.rs:359-360` returning `Some(SANGER)` for BAM.

**Internal consistency (review question 3):** the fix *removes* a contradiction rather than creating
one. The trimming report prints `Quality encoding type selected: ASCII+64` on both the SE and PE FASTQ
paths; at `d5f41c5` the FastQC report next to it said `Sanger / Illumina 1.9`. They now agree. No path
I probed tells R1 and R2 different things.

---

## 2. Findings

### HIGH — none.

### MEDIUM-1 — `--phred64 --fastqc` now aborts the run on data that is not really Phred+64, and `src/main.rs:332-336` says it cannot

`fastqc-rust` 1.1.1 treats a declared Phred+64 as a constraint, not a hint:
`phred::detect(lowest_char, true)` returns `Err` when the lowest observed quality character is below
ASCII 64 (`utils/phred.rs:67-72`), and `BasicStats::write_text_report` turns that into an
`io::Error` (`modules/basic_stats.rs:395-397`), which becomes a non-zero exit from `runner::run` and
therefore an `Err` out of `fastqc::run` (`src/fastqc.rs:77-78`).

Verified with a Phred+33 fixture (24×`I` then 10×`,`, lowest char 44):

```
d5f41c5  --phred64 --fastqc -q 0 --length 0   rc=0, report written (Sanger, 40.0/11.0)
8d5567a  --phred64 --fastqc -q 0 --length 0   rc=1
         Failed to process p33_low_trimmed.fq: Phred64 encoding is incompatible
         with having ASCII char ',' with value 44) in the file
         Error: FastQC analysis failed (exit code 1)
```

The trimmed FASTQ and both trimming reports are already on disk when this fires; only the FastQC
artefacts are missing. On the trimming path the default `-q 20` usually masks it (Phred+33 data read as
Phred+64 is trimmed to nothing, the output is empty, and `detect(NO_QUALITY_SEEN, true)` succeeds — I
confirmed a 0-sequence report with `Encoding Illumina 1.5` and `rc=0`).

**`--clump_only` has no such escape hatch.** It does no quality arithmetic, so the low characters
always survive to the output, and the failure is unconditional:

```
d5f41c5  --clump_only --phred64 --fastqc  (p33 fixture)   rc=0, report written
8d5567a  --clump_only --phred64 --fastqc  (p33 fixture)   rc=1, no FastQC artefacts
```

This makes the comment at `src/main.rs:332-336` false:

> `--hardtrim5/3` and `--clump_only` currently accept the flag as an inert no-op, so this does remove
> working invocations

`--hardtrim5/3` still do (`--fastqc` is refused there). `--clump_only` no longer does once `--fastqc`
is requested. The comment is the stated justification for rejecting `--phred64` uniformly on BAM input,
so it is worth keeping accurate.

I am not arguing the new behaviour is wrong — failing loudly on a mislabelled encoding is defensible,
and the upstream message names the offending character. But it is a new exit-1 path reached after the
primary output is complete, on a flag combination the codebase currently documents as inert, with no
CHANGELOG line and no `--help` note. Pick one of:

1. Amend `src/main.rs:333` to `"--hardtrim5/3 accept the flag as an inert no-op (and --clump_only does
   too, unless --fastqc is requested)"`, plus a CHANGELOG entry; or
2. Downgrade an encoding-resolution failure to a warning so a report problem does not fail a run whose
   trimmed output is correct — `src/fastqc.rs:77` is the single place to do it, and it would also cover
   the pre-existing "any FastQC failure fails the whole run" behaviour.

### MEDIUM-2 — the new test covers 1 of the 6 call sites that received a real value; `ubam_output_is_not_told_about_phred64` cannot fail

Mutation-tested against a pinned `git archive 8d5567a` copy with its own target dir, full suite each
time (baseline: **826 passed, 0 failed**):

| Mutation | Suite result |
|---|---|
| **M1** — all three `main.rs` uBAM sites `false` → `cli.phred64` | **826 passed, 0 failed** |
| **M2b** — PE FASTQ ×2 + passthrough + both `clump_only` FASTQ sites `cli.phred64` → `false` | **826 passed, 0 failed** |
| **M2a** — SE FASTQ site (`main.rs:1571`) `cli.phred64` → `false` | `phred64_reaches_the_fastqc_config_on_fastq_output` FAILED (`left: "Sanger / Illumina 1.9"`) |

Two things follow.

First, M2a confirms the commit message's mutation claim for the direction it checked. But M1 shows the
opposite direction is unconstrained: the five BAM-site `false` arguments
(`main.rs:2369`, `2528`, `2648`; `clump_only.rs:917`, `1177`) can all be flipped to `cli.phred64`
without a single test noticing. That is inherent — `.or()` makes BAM's own `Some(SANGER)` win either
way — so `ubam_output_is_not_told_about_phred64` is a characterization test of the dependency, not a
test of the arguments its name refers to. Worth keeping as an upstream tripwire, but the name and the
in-test comment ("this pins the contract the BAM call sites rely on",
`tests/integration_fastqc_phred64.rs:101-102`) oversell it. Say plainly that the argument is inert and
that the test pins `fastqc-rust`'s precedence.

Second, and more useful: M2b shows the **paired-end FASTQ path and both `--clump_only` FASTQ paths are
untested**. `clump_only.rs:295` and `:444` gained a parameter that none of the twelve unit-test call
sites exercises — every one passes `fastqc_requested: false`, so the new argument is dead in those
tests. Given that `--clump_only` is the path where this change has the largest user-visible effect
(MEDIUM-1), one PE case and one `--clump_only` case in
`tests/integration_fastqc_phred64.rs` would be cheap: both reuse `phred64_test.fastq` and assert the
same two numbers.

Also note the ordering hazard this created: in `src/clump_only.rs` the new `phred64` parameter was
inserted *before* `no_report_file` in the signature, while all twelve test call sites appended a new
`false` at the end. That happens to be harmless because both arguments are `false` at every site, but
the pre-existing literal silently changed meaning. A struct or builder would be the real answer for an
11-parameter function with four bools; short of that, named-argument comments at the call sites.

### MEDIUM-3 — the new test file states the upstream change backwards

`tests/integration_fastqc_phred64.rs:3-4`:

> FastQC 0.13.0 keys its encoding guess off the observed quality characters, so a Phred+64 file read as
> Phred+33 reports every mean 31 too high.

0.13.0 is precisely the release that **stopped** keying the guess off the observed characters. In
`fastqc-rust-1.1.1/src/utils/phred.rs:56-65`, `detect(lowest, false)` returns `SANGER` unconditionally;
the observed characters are used only for the feasibility check and to name Illumina 1.3 vs 1.5 once
Phred+64 has been declared. Upstream's own test comment says so: `// Formerly misdetected as Illumina
1.5` (`phred.rs:106-110`). The consequence stated in the second half is correct; the mechanism is
inverted, and the inversion is the whole reason this commit exists. Suggested replacement for the first
clause: *"FastQC 0.13.0 no longer guesses the encoding — it reads Phred+33 unless told otherwise — so a
Phred+64 file reports every mean 31 too high."*

The commit message carries a smaller version of the same problem: *"`phred_offset()` feeds only
`BamWriter`"* is not right. `cli.phred_offset()` also feeds `trimmer::TrimConfig.phred_offset`
(`src/main.rs:1280`, i.e. the quality-trimming arithmetic) and `report::TrimConfig.phred_encoding`
(`src/main.rs:1529`, `2328`, `2710`). The conclusion is still sound — I confirmed `FastqWriter` writes
`self.qual` verbatim (`src/fastq.rs:93`, no offset anywhere in the file) — but a future reader could
take the parenthetical to mean quality trimming ignores the offset, which it does not.

### LOW-1 — the doc-sweep commit missed `test_files/README.md:14`, which now describes an engine the repo no longer bundles

`1b16bed`'s message says ten places still asserted 0.12.1. There is an eleventh, and it matters because
it is the stated rationale for the very fixture the new test uses:

> …would additionally trip FastQC's encoding heuristic (upstream FastQC 0.12.1 behaviour, reproducible
> on either output format: it keys off the lowest observed quality character, and an all-Q40 record
> leaves nothing below ASCII 64, so the guess falls to Phred+64 and **every mean comes out 31 too low**)

Under the bundled 1.1.1 engine a uniform-Q40 record no longer trips anything, and the error direction
has inverted — a Phred+64 file now reads 31 too *high*, which is exactly what `8d5567a` fixes. Either
retarget the parenthetical at 0.13.0 or mark it explicitly as history.

(I checked the rest of the sweep: `CHANGELOG.md:1398/1420` and the synced
`docs/src/content/docs/reference/changelog.md:654/676` correctly keep 0.12.1 as a statement about what
v2.1.0 shipped, and `wiki/usage/install.md:58-61` is version-agnostic. No other stale reference.)

### LOW-2 — neither commit adds a CHANGELOG entry

`### Unreleased` is empty at `CHANGELOG.md:4`. Between them these three commits retarget the bundled
report at a new Java FastQC release, change `--fastqc_args`'s flag vocabulary (`--svg` → no-op, `--png`
added), change report content under `--phred64`, and add the exit-1 path in MEDIUM-1. The pin bump is
established context, but shipping it with no changelog line is the kind of thing the project's own
one-sentence-per-entry rule exists to catch. Three sentences would do it.

### LOW-3 — the `phred64` parameter doc landed on the wrong function

`src/fastqc.rs:23-39` is a parameter block written for `run`, but it sits immediately above
`pub fn resolved_threads` (line 40) and is therefore that function's rustdoc; `pub fn run` (line 55) has
no doc comment at all. The misattribution predates this commit, but `8d5567a` added lines 36-37 to it,
so `resolved_threads` now documents a `phred64` parameter it does not take. Move lines 23-37 down to
line 54 and leave lines 38-39 (the `resolved_threads` sentences) where they are.

### LOW-4 — `--fastqc_args "--phred64"` warns that the bundled engine does not understand a flag it does understand

`apply_fastqc_args` (`src/fastqc.rs:95-125`) has no `--phred64` arm, so the token falls to the
unknown-flag branch:

```
Warning: --fastqc_args flag '--phred64' is not translated to the bundled FastQC; ignored
```

`fastqc-rust`'s own CLI does have `--phred64` (`fastqc-rust-1.1.1/src/main.rs:134`), and now that
trim_galore wires the field, this is the one curated-subset gap where the right answer exists one flag
away. Either add `"--phred64" => config.phred64 = true` next to the `--png` arm the same sweep added
(idempotent with trim_galore's own flag), or leave it and accept that the message is literally true.

### LOW-5 — `basic_stat()`/`per_base_mean()` duplicate the archive-opening block

`tests/integration_fastqc_phred64.rs:20-31` and `:36-47` repeat the same five statements to get
`fastqc_data.txt` out of the zip, and the two differ gratuitously in their panic messages (`"opening
{path}: {e}"` vs `"report archive"`). One `fn report_data(zip, stem) -> String` helper, two thin parsers
on top.

---

## 3. Things I checked that are fine

- **Clippy, cold cache:** `cargo clippy --all-targets --release -- -D warnings` on a fresh target dir —
  exit 0, no warnings (27.6 s, full dependency graph compiled, so not a cache hit).
- **Tests:** `cargo test --release` — 826 passed, 0 failed across 21 test binaries, including both new
  cases.
- **`--help` output from the built binary** (not the source): `--fastqc` says `0.13.0-compatible`,
  `--fastqc_args` lists `--png` and notes `--svg` as a no-op. Accurate.
- **The `--svg` no-op claim is true:** with `--png` the HTML carries 29 PNG data-URIs and 0 SVG; without
  it, 21 PNG + 8 SVG. The 8 that flip are the graphs, so "SVG graphs are the default" is right.
- **`integration_preflight_tripwire.rs:31` is right too:** neither `--png` nor `--svg` writes any file
  beside `*_fastqc.html` / `*_fastqc.zip` — verified by listing the output dir for both.
- **The `zip` dev-dep claim:** `zip = { version = "9", default-features = false, features = ["deflate",
  "time"] }` matches `fastqc-rust`'s own spec character for character
  (`fastqc-rust-1.1.1/Cargo.toml:123-129`), and `Cargo.lock` gained one line in trim-galore's own
  dependency list, no new `[[package]]`.
- **All 11 `fastqc::run` call sites were updated**, and each gets the value its output format implies.
  (The commit message says "the three uBAM call sites"; there are five — the two in `clump_only.rs` have
  their own comment, so this is wording, not a gap.)
- **`cli.phred64` vs `cli.phred_offset()`:** equivalent here (`cli.rs:1222-1224`), and the direct bool
  is the better choice at a `bool` parameter.
- **No path where `cli.phred64` is true but the FASTQ output is not Phred+64**, other than the
  theoretical `--passthrough` case where a user hands a Phred+33 carrier file alongside Phred+64
  R1/R2. uBAM input is rejected at `main.rs:343`, the specialty modes reject `--fastqc`, and everything
  else writes quality verbatim.

## 4. Not verified

- **Behaviour on real Illumina 1.5 data.** Everything above uses the 4-read, 34 bp
  `phred64_test.fastq` and a hand-built Phred+33 counterpart. The per-base rows are ungrouped at that
  length; on a read long enough to trigger `BaseGroup` the test's `per_base_mean(zip, stem, "25")` would
  panic with `no row for base 25` rather than assert anything. Loud, so acceptable, but untested.
- **Report byte-identity against Java FastQC 0.13.0.** Out of scope for this review and, as
  `1b16bed` correctly now says in `CLAUDE.md`, not something CI checks either — the smoke test asserts
  four filenames are present in the `.zip`.
- **The CI `validation` job.** I did not install Perl Trim Galore 0.6.11 or run the md5 matrix; none of
  the paths it covers (SE, PE, hardtrim5, clock, demux) touch `fastqc::run`.
- **Whether `--phred64` reaching the `--clumpify` (not `--clump_only`) FastQC path is correct** beyond
  the code reading — it goes through `run_single_file`, the same site M2a proves is covered, so I did
  not probe it separately.
