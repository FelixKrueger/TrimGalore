# Code review B — PR #382, `fix(io): strip .bgz/.bgzf so output and report share a stem (#381)`

- **Reviewer:** B (independent; fresh context, own worktree at `pr382-b` = `ae55301`)
- **Base:** `dev` (`c6528e5`)
- **Files touched:** `src/io.rs`, `tests/integration_gzip_non_gz_extension.rs`, `CHANGELOG.md`
- **No repository file was modified by this review.** All suggested edits are quoted below for the maintainer to apply.

## Verdict

**MERGE.** The rewrite is behaviour-preserving for every filename that does not end in `.bgz`/`.bgzf` — verified by construction and by brute force, not by sampling — and the `.bgz`/`.bgzf` changes are exactly the ones #381 asks for. Nothing in the `validation` byte-identity matrix moves. `cargo test` (383 lib tests + 12 integration suites), `cargo clippy --all-targets --release -- -D warnings` and `cargo fmt --all -- --check` are all clean.

Nothing I found rises to blocking. There is **one factual overstatement in the CHANGELOG/PR description** (M1) worth fixing before merge, because a maintainer would otherwise rely on it, plus a short list of Low nits.

## The crux: old-vs-new equivalence

I extracted both implementations verbatim into a standalone harness (old = the function at `src/io.rs:435` on `dev`, via `git show dev:src/io.rs`; new = `src/io.rs:451-477` on the PR) and diffed their outputs. Full table:

| input | old | new | same |
|---|---|---|---|
| `sample.fastq.gz` | `sample` | `sample` | yes |
| `sample.fq.gz` | `sample` | `sample` | yes |
| `sample.fastq` | `sample` | `sample` | yes |
| `sample.fq` | `sample` | `sample` | yes |
| `sample.gz` | `sample` | `sample` | yes |
| `sample.txt.gz` | `sample` | `sample` | yes |
| `sample.bam` | `sample` | `sample` | yes |
| `sample` | `sample` | `sample` | yes |
| `sample.fastq.bz2` | `sample.fastq` | `sample.fastq` | yes |
| `archive.fq.tar.gz` | `archive.fq` | `archive.fq` | yes |
| `lib.2_R1.fastq.bgz` | `lib.2_R1.fastq` | `lib.2_R1` | **NO** — intended (#381) |
| `SAMPLE.FASTQ.GZ` | `SAMPLE.FASTQ` | `SAMPLE.FASTQ` | yes |
| `sample.FQ.gz` | `sample` | `sample` | yes |
| `.fq` | `` (empty) | `` (empty) | yes |
| `.fastq.gz` | `` (empty) | `` (empty) | yes |
| `<empty name>` | `` (empty) | `` (empty) | yes |
| `sample.fastq.bgz` | `sample.fastq` | `sample` | **NO** — intended (#381) |
| `sample.fq.bgz` | `sample.fq` | `sample` | **NO** — intended (#381) |
| `sample.fastq.bgzf` | `sample.fastq` | `sample` | **NO** — intended (#381) |
| `sample.fq.bgzf` | `sample.fq` | `sample` | **NO** — intended (#381) |
| `sample_R1.fq.bgz` | `sample_R1.fq` | `sample_R1` | **NO** — intended (#381) |
| `sample.bgz` | `sample` | `sample` | yes |
| `sample.bgzf` | `sample` | `sample` | yes |
| `sample.txt.bgz` | `sample.txt` | `sample` | **NO** — side effect, benign (see M1) |
| `sample.bam.bgz` | `sample.bam` | `sample` | **NO** — side effect, benign (see M1) |
| `sample.bgz.gz` | `sample` | `sample` | yes |
| `sample.fastq.bgz.bgzf` | `sample.fastq.bgz` | `sample.fastq` | **NO** — degenerate double suffix |
| `sample.fq.tar.bgz` | `sample.fq.tar` | `sample.fq` | **NO** — matches existing `.tar.gz` handling |
| `sample.fastqq.gz` | `sample` | `sample` | yes |
| `sample.gz.fastq.gz` | `sample.gz` | `sample.gz` | yes |
| `sample.fq.fastq.gz` | `sample.fq` | `sample.fq` | yes |
| `s.a.m.p.l.e.fastq.bgz` | `s.a.m.p.l.e.fastq` | `s.a.m.p.l.e` | **NO** — intended (dotty stem, #381) |
| `sample.bgzip` | `sample` | `sample` | yes |
| `sample.fastq.bgzip` | `sample.fastq` | `sample.fastq` | yes |
| `sample.fq.zst` | `sample.fq` | `sample.fq` | yes |
| `/data/dir.with.dots/sample.fastq.bgz` | `sample.fastq` | `sample` | **NO** — intended; the dotty *directory* is not consumed |
| `sample.fastq.BGZ` | `sample.fastq` | `sample.fastq` | yes — case-sensitive, see M4 |
| `..fq.bgz` | `..fq` | `.` | **NO** — pathological |
| `sample.` | `sample` | `sample` | yes |
| `sample.fq.` | `sample.fq` | `sample.fq` | yes |
| `.bgz` | `.bgz` | `` (empty) | **NO** — pathological dotfile |
| `.fastq.bgz` | `.fastq` | `` (empty) | **NO** — pathological dotfile; `.fastq.gz` already gave `""` on `dev` |
| `sample.fastq.gz.bgz` | `sample.fastq.gz` | `sample.fastq` | **NO** — degenerate double suffix |

A 43-row table is still sampling, so I then brute-forced it: 6 stems (`sample`, `s.a`, `x`, `sample_R1`, `""`, `.`) × up to three extension tokens drawn from `{.fastq, .fq, .gz, .bgz, .bgzf, .txt, .bam, .tar, .FQ, .GZ, .BGZ, .bz2, ""}` = **13,182 filenames**:

```
cases checked: 13182
differences on names ending .bgz/.bgzf (expected): 2118
differences on names NOT ending .bgz/.bgzf (must be 0): 0
```

Uppercase `.BGZ`-terminated names were deliberately classified as *unexpected*, and none appeared — confirming the new strip is case-sensitive exactly as the old one was, so `SAMPLE.FASTQ.GZ` still yields `SAMPLE.FASTQ` and Perl-parity on casing is untouched.

**Conclusion: no non-`.bgz` name changes.** The change is confined to the spelling #381 is about.

Two structural facts, because they are what makes the rewrite safe rather than merely untested:

- The two lists are genuinely disjoint, as the comment at `src/io.rs:458-459` claims: `.bgz` does not end in `.gz` (its last three characters are `bgz`), and `.bgzf` ends in neither. So `find_map` order is immaterial.
- The old ordered list's "order of specificity" was never load-bearing either — `sample.fastq.gz` does not end in `.fq.gz` (the six-character tail is `stq.gz`). Both versions are order-independent, which is why splitting the list could not reshuffle any `.gz` case.

## Verification of the PR's own claims

| Claim | Verdict | Evidence |
|---|---|---|
| "Output *filenames* only" | **Almost** — one behavioural exception | see M2 |
| "Nothing in `test_files/` … is named `.bgz`" | **True** | `ls test_files \| grep -i bgz` → empty |
| "nothing in the `validation` matrix is named `.bgz`" | **True** | `grep -n bgz .github/workflows/*.yml` → no hits; the validation job's inputs are `illumina_10K.fastq.gz`, `BS-seq_10K_R{1,2}.fastq.gz`, `clock_10K_R{1,2}.fastq.gz`, `demux_test.fastq.gz`, `ubam_test.bam` (`.github/workflows/ci.yml:383,394,709,718,729,943`) — every one of those stems is an unchanged row above |
| "the Perl-0.6.11 byte-identity comparison does not move" | **True** | follows from the row above; the hardtrim5 oracle even asserts the literal name `illumina_10K.30bp_5prime.fq.gz`, whose stem is table row 1 |
| "`.bam` inputs … keep the stems they had" | **True** | table row 7; `test_clumped_bam_output_name_from_bam_input` (`src/io.rs:807`) still green |
| "bare `.bgz`/`.bgzf` already worked via `file_stem`" | **True** | table rows `sample.bgz`, `sample.bgzf` |
| "Names with no inner `.fastq`/`.fq` … keep the stems they had" | **False as written** | `sample.txt.bgz`, `sample.bam.bgz` — see M1 |
| `cargo test` green, clippy clean, fmt clean | **True** | `cargo test` exit 0, 13 suites, 383 lib tests; clippy exit 0; fmt exit 0. Exit codes captured **unpiped** — a `\| tail` would have masked a failure. |

## Blast radius

`strip_fastq_extensions` has 20 call sites, all of them name construction; none feeds a decision, comparison, or parse.

- FASTQ trim outputs — `single_end_output_name` (`src/io.rs:90`), `paired_end_output_names` (`src/io.rs:179-180`), `unpaired_output_names` (`src/io.rs:210-211`), `passthrough_output_name` (`src/io.rs:245`).
- `--clump_only` — `clumped_output_name` (`src/io.rs:283`), `clumped_paired_output_names` (`src/io.rs:315-316`).
- uBAM outputs — `single_end_bam_output_name` (`src/io.rs:116`), `paired_bam_output_name` (`src/io.rs:143`), `clumped_bam_output_name` (`src/io.rs:343`), `clumped_paired_bam_output_name` (`src/io.rs:369`). A `.bam` input reaches all four and its stem is unchanged (table row 7), so `--output-format ubam` names do not move.
- Specialty modes — `hardtrim_output_name` (`src/specialty.rs:430`), `hardtrim_bam_output_name` (`src/specialty.rs:452`), `clock_output_name` (`src/specialty.rs:466`), `implicon_output_name` (`src/specialty.rs:482`).

Report namers (`report_name`, `json_report_name`, `clumping_report_name`, `src/io.rs:400-430`) do **not** use it — they use the whole input filename, which is why the divergence existed and why the report side is deliberately unchanged.

The specialty modes are the piece the PR description does not mention. Their names change for `.bgz` input too, in the same direction. Confirmed against the built binary:

```
$ trim_galore --hardtrim5 30 -o out4 sample.fastq.bgz
$ ls out4
sample.30bp_5prime.fq          # dev would have produced sample.fastq.30bp_5prime.fq
```

That is the desirable direction and consistent with the fix; it is simply unmentioned and untested. `--implicon`'s extra `_R1`/`_R2` strip (`src/specialty.rs:485-498`) operates on the already-stripped stem, so it is unaffected in kind.

## Issues

### M1 — Medium (documentation accuracy). "keep the stems they had" is not true for `<name>.<ext>.bgz`

`CHANGELOG.md` (new entry) and the PR description both say:

> Names with no inner `.fastq`/`.fq` (`sample.bgz`), `.bam` inputs and unrelated extensions keep the stems they had.

Counterexamples from the table: `sample.txt.bgz` was `sample.txt`, is now `sample`; `sample.bam.bgz` was `sample.bam`, is now `sample`. Cause: the gzip-family strip runs first, so the `file_stem()` fallback then sees `sample.txt` and drops a second component.

This is *benign* — it is exactly what `sample.txt.gz` and `sample.bam.gz` have always done, so it makes `.bgz` consistent with `.gz`, which is the point of the PR. But the sentence as written is the kind of claim a maintainer relies on instead of re-deriving the table. The accurate boundary is "names with no gzip-family suffix".

Suggested edit (`CHANGELOG.md`, replacing that one sentence in the new `#### Bug fixes` entry):

```markdown
  `.bgz`/`.bgzf` is now treated exactly as `.gz` always has been, so
  `sample.txt.bgz` yields `sample` just as `sample.txt.gz` does. Names with no
  gzip-family suffix — `.bam` inputs, a bare stem, unrelated extensions — keep
  the stems they had.
```

`src/io.rs:449-450` carries a milder version of the same imprecision ("Anything else … falls through to `Path::file_stem`, which drops the final component only" — true of the fallback, but after a `.bgz` strip two components are gone in total). If L6 is taken, that paragraph is rewritten anyway.

### M2 — Medium/Low. The one place the change is not "filenames only"

A `--paired` invocation that succeeded on `dev` now hard-errors, because the newly-shared stem trips the output-collision pre-flight (`src/main.rs:706-724`):

```
$ trim_galore --paired -o out6 p_R1.fastq.bgz p_R2.fastq.bgz p_R1.fq.bgz p_R2.fq.bgz
Error: Output path collision (case-insensitive, for APFS/NTFS safety):
out6/p_R1_val_1.fq and out6/p_R1_val_1.fq would be written to the same file.
exit=1
```

On `dev` the two pairs produced `p_R1.fastq_val_1.fq` and `p_R1.fq_val_1.fq` and the run completed. Contrived (it needs the same sample present under both `.fastq.bgz` and `.fq.bgz`), loud rather than silent, and the `.gz` spelling has always behaved this way — I would not block on it. But "Output *filenames* only" is strictly not true. One clause covers it if you want the CHANGELOG airtight: *"a multi-input run whose inputs differ only in `.fastq` vs `.fq` before the `.bgz` suffix now shares one output name, so paired mode's collision pre-flight fires where it previously did not — the same behaviour the `.gz` spelling has always had."*

### M3 — Low. Pre-existing gap that one more spelling now reaches: no SE-trim collision pre-flight

The same collision on the **single-end trim** path is silent:

```
$ trim_galore -o out1 sample.fastq.bgz sample.fq.bgz
exit=0
out1/sample_trimmed.fq                        # written twice; the second input wins
out1/sample.fastq.bgz_trimming_report.{txt,json}
out1/sample.fq.bgz_trimming_report.{txt,json}
```

I checked whether this is the PR's doing. It is not — the exact `.gz` analogue already behaves identically on this same build:

```
$ trim_galore -o out5 s2.fastq.gz s2.fq.gz
exit=0 → one s2_trimmed.fq.gz, two reports, no error
```

So the hazard class predates #382; the PR extends it to the `.bgz` spelling. Worth a **separate issue**, because the asymmetry looks unintentional: `--clump_only` SE has a case-folded pre-flight (`src/main.rs:481-499`), SE uBAM has one (`src/main.rs:609-615`), `--paired` has one (`src/main.rs:706`), specialty multi-pair has one (`src/main.rs:2410`) — the plain single-end FASTQ trim loop is the only dispatch path without one. That is issue #216's guard missing from one branch, not a #382 regression. **Do not gate this PR on it.**

### M4 — Low. The fix is spelling-sensitive: `.BGZ` still reproduces #381

```
$ trim_galore --dont_gzip -o out3 sample.fastq.BGZ
$ ls out3
sample.fastq_trimmed.fq                 # the #381 symptom, verbatim
sample.fastq.BGZ_trimming_report.txt
```

Since #374 the reader is content-based, so an uppercase-named bgzip file reads fine and lands in exactly the state #381 describes. I am **not** recommending a case-insensitive match: `.gz` matching is case-sensitive too (`SAMPLE.FASTQ.GZ` → `SAMPLE.FASTQ`, old and new alike), that matches Perl 0.6.11, and folding the whole family would silently rename outputs for anyone with uppercase inputs — its own discussion, like the compression-mirroring limitation. What I would do is pin it, so the next person who "obviously" adds `to_ascii_lowercase()` meets the `.GZ` consequence first. Assert in L5.

### L5 — Low (test sufficiency). Four table rows worth pinning

The three new tests are genuinely capable of failing:

- `test_strip_fastq_extensions_bgz` (`src/io.rs:505`) — on `dev` the first assert gets `sample.fastq` for `sample.fastq.bgz`; my harness reproduces exactly the `left: "sample.fastq", right: "sample"` the author reports.
- `test_strip_fastq_extensions_leaves_other_names_alone` (`src/io.rs:534`) — a real regression net for the four rows it covers.
- `bgz_output_and_report_agree_on_the_stem` (`tests/integration_gzip_non_gz_extension.rs:305`) — on `dev`, `dir.join("sample_trimmed.fq").exists()` is false because the file is `sample.fastq_trimmed.fq`; the assert message dumps the directory listing. Present and green in my run.

Left unpinned, in descending order of how likely a future rewrite is to move them unnoticed:

1. `sample.txt.bgz` / `sample.bam.bgz` — the *only* newly-changed rows that are not FASTQ names, i.e. exactly the ones M1 mis-describes.
2. Uppercase `.BGZ` / `.GZ` — currently undocumented in code (M4).
3. `archive.fq.tar.gz` → `archive.fq` and `sample.fastq.bz2` → `sample.fastq` — the two "surprising but longstanding" rows a future `strip_suffix` shuffle could move.

Suggested addition to `test_strip_fastq_extensions_leaves_other_names_alone` (`src/io.rs:534-539`):

```rust
        // `.bgz`/`.bgzf` is treated exactly as `.gz`: a second component goes
        // with it, same as `sample.txt.gz`.
        assert_eq!(strip_fastq_extensions(Path::new("sample.txt.bgz")), "sample");
        // Case-sensitive, as `.gz` has always been — `SAMPLE.FASTQ.GZ` keeping
        // `SAMPLE.FASTQ` is Perl-0.6.11 behaviour, so `.BGZ` matches it.
        assert_eq!(
            strip_fastq_extensions(Path::new("sample.fastq.BGZ")),
            "sample.fastq"
        );
        // Longstanding, and surprising enough to pin: only one compression
        // suffix is ever stripped, so a non-gzip compressor leaves the inner
        // extension in place.
        assert_eq!(
            strip_fastq_extensions(Path::new("sample.fastq.bz2")),
            "sample.fastq"
        );
```

Also unpinned, and cheaper to leave so: specialty-mode naming on `.bgz`. One assert against `hardtrim_output_name(Path::new("sample.fastq.bgz"), 30, "5prime", None, false)` == `sample.30bp_5prime.fq` in `src/specialty.rs`'s test module would cover the family, since all four share the helper. Optional.

### L6 — Low (convention). The 19-line doc comment does not fully earn its length

`src/io.rs:432-450` is 19 doc lines for a 27-line function. Against `CLAUDE.md`'s "state the fact, not the evidence … detail that genuinely matters goes in the commit message, PR description, or a linked issue": the first paragraph (what it handles, why the groups are split) and the last (the fallback) are load-bearing. The middle paragraph (`src/io.rs:440-447`) retells #381's narrative, which already exists in four other places in this same diff — the issue, the CHANGELOG entry, the unit-test doc (`src/io.rs:498-503`), and the integration-test doc (`tests/integration_gzip_non_gz_extension.rs:291-303`). Five copies of one story is where a comment stops paying rent.

House style in `io.rs` is generous with doc comments (`is_gzipped` at `src/io.rs:12-28`, `ensure_output_dir` at `src/io.rs:48-64`), so this is not out of place and I would not block on it. Suggested trim:

```rust
/// Strip FASTQ extensions from a filename, returning just the base stem.
///
/// Handles `.fastq` / `.fq`, each optionally followed by `.gz`, `.bgz` or
/// `.bgzf`. The two groups are stripped in sequence rather than enumerated as
/// combined suffixes, so a new compression name costs one list entry (#381).
///
/// Anything else — a `.bam` input, an unrelated extension — falls through to
/// `Path::file_stem`.
```

The two-line inline comment at `src/io.rs:458-459` *does* earn its length: "these are disjoint, so order does not matter" is a non-obvious fact a reader would otherwise have to re-derive.

### L7 — Low (nit). "longest first" implies a constraint that does not exist

`src/io.rs:465`: `// Then the FASTQ extension itself, longest first.` — `.fastq` does not end in `.fq`, so this list is order-independent too, which the comment six lines above goes out of its way to establish for the other list. As written it invites the next reader to believe the order is load-bearing. Suggested edit:

```rust
    // Then the FASTQ extension itself.
```

## Test-rewrite direction (review item 4) — correct, nothing weakened

The five assertions moved from #374's `sample.fq_trimmed` / `pair_R{1,2}.fq_val_{1,2}` to `sample_trimmed` / `pair_R{1,2}_val_{1,2}`. That is the direction #381 mandates. They cannot silently pass at the new stem: `read_output_raw` (`tests/integration_gzip_non_gz_extension.rs:95-108`) tries `{stem}.fq` then `{stem}.fq.gz` and **panics** with `no output at … or …` if neither exists, and `read_output` additionally rejects an empty file ("every read was filtered, so this test would pass without checking anything"). The `.contains("@pe1_0")`-style content assertions were left untouched. The module doc (`tests/integration_gzip_non_gz_extension.rs:79-82`) moved from "known cosmetic wart" to a pointer at the new test — the right bookkeeping.

One small imprecision in the new test's doc (`tests/integration_gzip_non_gz_extension.rs:302-303`): "the assertion is the same for both: this is a contrast test, not a hard-coded expectation" — the expected stem *is* hard-coded (`"sample"` on all three rows). The parity is real (all three names must produce the same stem), so the test is sound; the sentence oversells it. Not worth a round trip on its own.

## CHANGELOG placement (review item 5) — correct

The entry lands in `### Unreleased` → `#### Bug fixes` (`CHANGELOG.md:33`), immediately after the #374 entry it amends, and it correctly demotes #374's "Two known limitations" to "One known limitation remains". Right section for a fix; `#### Changes` begins after it.

`docs/src/content/docs/reference/changelog.md` is a hand-synced copy whose `Unreleased` section has no `#### Bug fixes` block at all — it does not carry #374's entry either, so leaving it alone matches the precedent set by the PR this one amends. No action.

## Informational

- **`.bgzf` inclusion** answers the issue's open question at the cost of one array element. `bgzip` itself defaults to `.gz`; `.bgz`/`.bgzf` come from the surrounding tooling. Reasonable; maintainer's call.
- **The MultiQC consequence is asserted, not measured**, in both #381 ("I have reproduced the filenames but not run MultiQC over them") and the PR. Nothing in the change or its tests depends on it — the durable justification is internal consistency with the `.gz` spelling, which the PR does establish and pin.
- **The documented asymmetry still bites**: `io::is_gzipped` (`src/io.rs:29`) remains filename-based, so a `.bgz` input still yields plain `sample_trimmed.fq` (confirmed). A pipeline globbing `*_trimmed.fq.gz` still misses it. Correctly out of scope and still recorded in the CHANGELOG.
- **Empty-stem dotfiles**: a file literally named `.fastq.bgz` now yields an empty stem (output `_trimmed.fq`). `dev` already did this for `.fastq.gz` → `""`, so the hazard is pre-existing and one more spelling now reaches it. Pathological; no action recommended.
- **Spotted in passing, unrelated to this PR**: the collision pre-flight messages at `src/main.rs:721` and `src/main.rs:1815` suggest `` `--output-dir` ``, but the flag is declared `long = "output_dir"` (`src/cli.rs:152`). Wrong flag name in a user-facing error.

## Fixes recommended (none applied — report-only review)

| # | Priority | Where | Action |
|---|---|---|---|
| M1 | Medium | `CHANGELOG.md`, new entry | Replace "Names with no inner `.fastq`/`.fq` … keep the stems they had" with the gzip-family-suffix wording above. The claim is currently false for `sample.txt.bgz` / `sample.bam.bgz`. |
| M2 | Medium/Low | `CHANGELOG.md`, new entry | Optional clause noting the paired-mode collision pre-flight can now fire where it previously did not. |
| M3 | Low | separate issue | Single-end FASTQ trim loop is the one dispatch path with no output-collision pre-flight (compare `src/main.rs:481`, `:609`, `:706`, `:2410`). Pre-existing; file separately, do not gate #382. |
| M4 | Low | `src/io.rs` tests | Pin case-sensitivity rather than change it — assert in L5. |
| L5 | Low | `src/io.rs:534-539` | Add the three asserts quoted above (`.txt.bgz`, `.fastq.BGZ`, `.fastq.bz2`); optionally one specialty-mode assert in `src/specialty.rs`. |
| L6 | Low | `src/io.rs:432-450` | Trim the doc comment's middle paragraph; replacement quoted above. |
| L7 | Low | `src/io.rs:465` | Drop "longest first". |

## Checks run

| Check | Result |
|---|---|
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0 |
| `cargo test` (crate root, exit code captured unpiped) | exit 0 — 13 suites, 383 lib tests, 0 failures; `test_strip_fastq_extensions_bgz`, `test_strip_fastq_extensions_leaves_other_names_alone`, `bgz_output_and_report_agree_on_the_stem`, `single_end_bgz_input_is_decompressed` all `ok` |
| old-vs-new harness, 43 hand-picked names | differences confined to `.bgz`/`.bgzf` |
| old-vs-new brute force, 13,182 names | 0 differences outside `.bgz`/`.bgzf` |
| binary repro: #381 scenario, uppercase `.BGZ`, SE stem collision, `.gz` control, PE pre-flight, `--hardtrim5` | as reported above |
| `test_files/` and `.github/workflows/` scanned for `.bgz` | no hits — byte-identity claim holds |
