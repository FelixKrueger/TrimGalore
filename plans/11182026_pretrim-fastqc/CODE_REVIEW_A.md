# CODE REVIEW A — `--fastqc_raw`

**Target:** `git diff origin/dev..feat/fastqc-raw` (base `61ae72c`, head `cf3bc84`)
**Reviewer:** A (of two independent reviewers)
**Date:** 2026-11-18
**Verification:** release build, `cargo clippy --all-targets --release -- -D warnings` (forced recheck,
exit 0, zero diagnostics), `cargo fmt --all -- --check` clean, `--lib` 557 pass,
`integration_fastqc_raw` 7 pass, `integration_preflight_tripwire` 57 pass, plus 14 hand-run probes
against the release binary.

> **Note on the working tree.** While I was verifying, someone applied uncommitted edits to
> `src/fastqc.rs`, `src/main.rs`, `tests/integration_fastqc_raw.rs` and
> `tests/integration_preflight_tripwire.rs` — three of them overlap with fixes I had queued (see
> §"Fixes applied"). I made **no source edits**, to avoid conflicting with concurrent work in a
> shared tree. Findings below were re-verified against the tree **as it stands with those edits in
> place** unless stated otherwise.

---

## Summary

The core design is sound and I verified its load-bearing mechanism directly rather than taking it on
trust: a uBAM copied to `masquerade.fastq.gz` is QC'd as BAM (10 sequences, 65 bp) because
`SeqFormat::Bam` is declared, and an extensionless gzipped FASTQ is still decompressed (10 000 × 150 bp)
because fastqc-rust sniffs gzip by magic rather than by extension. The `_raw` alias does put the
intended string in the `Filename` field. All 13 pre-flight sites route through the wrapper, exactly one
`preflight_output_collisions` call survives in `main.rs`, every dispatch terminal in `main()` reaches a
wrapper, and `#[track_caller]` really does chain through it — `every_site_is_reached_by_some_case`
passes, which is the empirical proof, not the author's say-so.

The defect that matters is not in any of that. **The staged alias and the raw reports are new writers
whose paths never enter the pre-flight's candidate list, and their directory is chosen by a different
rule than every other output's.** On an ordinary multi-sample command line with no `--output_dir`, that
silently discards all but one raw report, exit code 0. This is the same shape as #383 / #388 / #391 /
#409 — the family `integration_preflight_tripwire` exists to catch — and the tripwire does not cover the
new flag at all.

**Verdict: changes needed before this gates the nf-core bump.** H1 is a silent-wrong-output bug on a
plausible command line. H2 is the missing guard that would have caught it. H3 makes the documented
specialty-mode support unusable the moment `--fastqc_args` is added, with an error naming a flag the
user never passed.

---

## Issues by area

### Logic

**H1 (High) — raw reports silently overwrite each other when `--output_dir` is absent.**

`run_pretrim_fastqc` pins the alias/report directory:

```rust
let dir = output_dir.unwrap_or_else(|| Path::new("."));
```

Every other output namer in the codebase, when `output_dir` is `None`, names relative to the *input's
own parent* (`naming::pair_output_dir`, `src/io.rs:121-125`, and the SE/PE trim namers). So the raw
reports are the only artifact whose directory is flat. Two inputs sharing a basename across
directories produce distinct trimmed outputs — the pre-flight passes, correctly — and one raw report.

Reproduced, single-end:

```
$ ls a b
a/s.fastq.gz   (illumina_10K,   10 000 reads, 150 bp)
b/s.fastq.gz   (smallRNA_100K, 100 000 reads,  51 bp)
$ trim_galore --fastqc_raw a/s.fastq.gz b/s.fastq.gz ; echo $?
0
$ ls
a/  b/  s_raw_fastqc.html  s_raw_fastqc.zip
$ unzip -p s_raw_fastqc.zip '*/fastqc_data.txt' | grep 'Total Sequences'
Total Sequences	100000          # ← pair b only; pair a's report is gone
$ ls a
s_trimmed.fq.gz  s.fastq.gz  s.fastq.gz_trimming_report.txt  s.fastq.gz_trimming_report.json
```

and paired-end, same result (`--paired a/s_R1 a/s_R2 b/s_R1 b/s_R2`, exit 0, one `s_R1_raw_fastqc.zip`
reporting 100 000 × 51 bp).

Two separate problems compound here:

1. **Wrong directory.** The trimmed output and the trimming report land in `a/` and `b/`; the raw
   reports land in the CWD. Even with no collision, a MultiQC run scoped to the trimmed output's
   directory sees only the trimmed set — which defeats the stated goal of "one invocation produces
   both report sets". Specialty modes are *not* affected: `--hardtrim5` writes its primary to the CWD
   too (verified: input `a/s.fastq.gz` → `./s.20bp_5prime.fq.gz` + `./s_raw_fastqc.zip`), so they are
   co-located there and a basename collision is caught by the pre-flight.
2. **Outside the pre-flight.** Neither the alias path nor the report paths are in `planned`, so the
   mechanism that exists precisely to refuse this cannot see it. The irony is structural: the author
   built `preflight_then_raw_qc` so the raw pass could not run *before* the refusal, but the paths the
   raw pass writes are not among the things the refusal reasons about.

The narrow fix is to derive the alias directory the way the mode's primary output is derived
(input-parent for the trim paths, CWD for the specialty paths) — `fastqc::run` then follows
automatically, since it is already handed the same `output_dir` Option and falls back to the analysed
file's parent. The durable fix, and the one this repo's history argues for, is to put the alias and the
two report paths into the `planned` list so the pre-flight owns them. I'd do both; either alone leaves
one of the two problems standing. **Recommend, not fixed** — it is a naming-convention decision across
five modes.

**H3 (High) — `--fastqc_args` turns the advertised specialty-mode support into a refusal, and the error
names a flag the user did not pass.**

```
$ trim_galore --hardtrim5 20 --fastqc_raw --fastqc_args "--nogroup" a.fastq.gz ; echo $?
Error: --fastqc is not supported with --hardtrim5, which writes no QC report; run FastQC
separately on the hardtrimmed output, or drop --fastqc
1
```

`--fastqc` appears nowhere on that command line. `fastqc_requested()` is
`self.fastqc || self.fastqc_args.is_some()` (`src/cli.rs:617-619`) and the #421 refusals key on it.
`run_pretrim_fastqc` *does* pass `cli.fastqc_args` through to `fastqc::run`, so the flag genuinely
configures the raw pass — there is simply no way to configure it in the four modes the `--help` text
advertises as supported ("Accepted with the specialty modes that reject --fastqc"), the flags doc
advertises, and the wiki table advertises.

The new unit test `fastqc_raw_accepted_with_every_specialty_mode` asserts exactly the documented claim
across all four modes — and never combines `--fastqc_raw` with `--fastqc_args`, so the suite endorses
the claim without testing the input that breaks it. **Recommend** (the fix is a judgement call: either
narrow the #421 refusals to "`--fastqc` explicitly set, or `--fastqc_args` without `--fastqc_raw`", or
retract the documented support).

**M1 (Medium) — the `_raw` alias can collide with a real input, and the diagnostic blames the user's
own file.**

```
$ trim_galore --fastqc_raw s_raw.fastq.gz s.fastq.gz ; echo $?
Error: Cannot stage the raw-FastQC input alias './s_raw.fastq.gz': a file is already there.
Remove it, or rename the input.
1
$ ls
s_raw_raw_fastqc.html  s_raw_raw_fastqc.zip  s_raw.fastq.gz  s.fastq.gz
```

The occupied path is input #2. Two consequences: the "a refused run creates nothing" invariant that
`preflight_then_raw_qc` was built to protect is broken (input #1's report survives a failed run), and
the remediation instructs the user to delete an input. The reverse argument order is refused cleanly
before any report, so the behaviour is order-dependent. Folding the alias paths into the pre-flight
(H1's second fix) resolves this too, with the pre-flight's own much better diagnostics.

**M4 (Medium) — `--passthrough` inputs get no raw report.**

```
$ trim_galore --fastqc_raw --paired --passthrough p_I1.fastq.gz p_R1.fastq.gz p_R2.fastq.gz
$ ls | grep fastqc
p_R1_raw_fastqc.html  p_R1_raw_fastqc.zip  p_R2_raw_fastqc.html  p_R2_raw_fastqc.zip
```

`run_pretrim_fastqc` iterates `&cli.input` only. Both neighbouring code paths disagree: the `--phred64`
pre-flight iterates `cli.input.iter().chain(cli.passthrough.iter())` (`src/main.rs:500-506`), and
`--fastqc` does QC the passthrough *output* (`src/main.rs:2137`, with a comment acknowledging the
report will look bad). Nothing documents the omission. If it is deliberate — a 16 bp uniform barcode
read gives a meaningless report — say so in the `--fastqc_raw` help text; if not, extend the loop.

### Errors

**M2 (Medium) — no unwind-safe cleanup for the staged symlink.** `let _ = remove_file(&alias); result?;`
covers `Ok` and `Err`, but not a panic. `fastqc-rust` is a third-party library running rayon workers
over attacker-shaped input; a panic there unwinds past the `remove_file` and strands the alias, which
then blocks *every* later `--fastqc_raw` run in that directory with "a file is already there". A drop
guard (`struct AliasGuard(PathBuf)` with `Drop` calling `remove_file`) makes the cleanup
unwind-correct and is strictly simpler than the current two-step. **Inferred** — I did not find an
input that panics fastqc-rust 1.1.1, so I could not reproduce it. SIGINT/SIGKILL leaves the alias
stranded regardless; that is inherent and the occupancy check handles the next run gracefully (verified
after the working-tree `symlink_metadata` fix).

**M3 (Medium) — TOCTOU between the occupancy check and `symlink()`.** Two `--fastqc_raw` runs over the
same input in the same directory both pass `symlink_metadata`, one wins `symlink()` and the loser dies
with a bare `File exists (os error 17)`; worse, the winner's `remove_file` can unlink the alias before
the loser's FastQC opens it. Concurrent runs in one directory already race on the trimmed outputs, so
this is not a new class — but the alias is the one path where TG creates *and deletes* a file whose
name it does not own. An `O_EXCL`-style create (`symlink` already is one; drop the `exists` pre-check
and map `EEXIST` to the friendly message) collapses the window and the duplicated logic. **Inferred.**

### Structure

**M5 (Medium) — `raw_alias_path` is untested and its fallback branch is mis-commented.** `main.rs`
carries no `#[cfg(test)] mod tests` (verified by grep), so all three new functions are reachable only
through the binary, and no test exercises `raw_alias_path` directly. Its `else` branch is commented
"Not a FASTQ-family name — a content-detected uBAM may be called anything", but
`strip_fastq_extensions` already falls back to `file_stem()` (`src/io.rs:837-841`), so a uBAM named
`sample.bam` takes the *if* branch (verified: alias `sample_raw.bam`). The `else` is reached only by a
name with no dot at all, or a leading-dot name (`.hidden` → `_raw.hidden`). Moving the three functions
into `src/lib.rs` would make the table of cases unit-testable; by hand I confirmed no-extension
(`reads` → `reads_raw`, report correct), `.fq.bgz`, multi-dot, and uBAM-with-FASTQ-name all behave.

The byte-slice `&name[stem.len()..]` is safe: both strings come from the same `to_string_lossy()` of
the same `file_name()`, and `strip_suffix_ignore_ascii_case`/`file_stem` both return prefixes, so the
index is always a char boundary — including for multibyte and non-UTF-8 names. I checked this rather
than assuming it.

### Efficiency

**L5 (Low) — the raw pass is strictly serial across inputs.** `for input in &cli.input` runs one
complete FastQC per file before trimming starts. FastQC is internally threaded, so throughput is
fine, but at `--cores 1` a four-file run does four serial full-file passes ahead of any output. The
help text's "Costs one extra pass over each input, run before trimming starts" already says this
honestly. No change recommended.

Deviation §13.4 (re-detect the format rather than thread `input_formats` through six drivers) is a good
trade and I confirmed the cost claim: `detect_input_format` decompresses one block, against a full-file
QC pass.

---

## Tests

**H2 (High) — `--fastqc_raw` is absent from the pre-flight tripwire.**

`fastqc_active` (`tests/integration_preflight_tripwire.rs:233-236`) matches by exact equality:

```rust
.any(|a| *a == "--fastqc" || *a == "--fastqc_args")
```

so `--fastqc_raw` does not enable the `_fastqc.{html,zip}` allowance, and **no tripwire case passes
`--fastqc_raw`** — not one of the 13 shapes in `every_site_is_reached_by_some_case`, not one of the
~40 `Tripwire::new` cases. The feature adds two new files per input, plus a transient symlink, to the
output directory on every one of the 13 arms, and none of it is under the tripwire whose module doc
says: *"Every instance was a writer whose path never entered the pre-flight's candidate list, so the
pre-flight silently stopped covering that output."*

The exact-match is fortunate: because `is_fastqc_artifact` *does* match `_raw_fastqc.zip` while
`fastqc_active` does not light up, adding a `--fastqc_raw` case today goes **red** rather than silently
passing — which forces the decision rather than hiding it. That makes the fix small (add
`|| *a == "--fastqc_raw"` to `fastqc_active`, add at least an SE and a PE case) and high-value: a
`--fastqc_raw` case on two same-basename inputs is exactly how H1 would have been caught.

**The seven new tests do fail on the mistakes they exist to catch.** I checked each against the
earlier-round concern about validations that cannot fail:

- `single_end_raw_report_names_the_sample_raw` — if FastQC were pointed at the input directly, the
  `zip.exists()` assertion fires first (report would be `illumina_10K_fastqc.zip`); the `Filename` and
  150 bp assertions then pin the two separate claims. Sound.
- `a_collision_refused_run_writes_no_raw_report` — I traced the inverted ordering by hand: QC-before-
  pre-flight would write `same_name_raw_fastqc.zip` (second input's alias succeeds, since the first
  was removed) and the `!exists()` assertion fires. Sound, and it is the only direct test of the
  wrapper's ordering guarantee.
- `a_bad_second_input_is_refused_before_any_report_is_written` — removing the per-input sanity loop
  writes R1's report while the paired driver still fails the run, so the assertion fires. Sound.
- `mid_file_truncation_surfaces_from_the_qc_pass` asserts a *weakness* (R1's report does get written).
  That is unusual but deliberate and well-explained, and its message tells a future reader what going
  red means. Iteration-log entry #3 is honest about having corrected an unachievable requirement here;
  I re-confirmed `sanity_check_any` passes on `truncated.fq.gz`.
- **L3 (Low)** `the_staged_alias_does_not_survive_the_run` uses `Path::exists()`, which follows
  symlinks and so reads `false` for a *dangling* alias. The target exists in this test, so it does
  catch a leak today; `symlink_metadata().is_err()` keeps catching it if the target ever moves.

**Not covered by any test** (beyond H2): uBAM input with `--fastqc_raw` — the `SeqFormat` parameter's
entire reason for existing. I verified it by hand (both `sample.bam` and the content-detected
`masquerade.fastq.gz`), but nothing in the suite pins it, so a future `SeqFormat`/`None` transposition
at the 12 threaded call sites would go unnoticed. Also untested: the four specialty modes actually
*producing* raw reports (the CLI tests only assert `validate().is_ok()`; I confirmed `--hardtrim5` and
`--clock` do produce them), the alias-occupied bail, and `raw_alias_path` itself.

**L2 (Low) — the inventory assertions are brittle in a way worth knowing.** `src.matches(...).count() - 1
// minus the definition` silently absorbs one stray textual occurrence of `preflight_then_raw_qc(` in a
`main.rs` comment, and a comment mentioning `preflight_output_collisions(` breaks the second
assertion. For a deliberately brittle tripwire that is acceptable; the working-tree widening from
`naming::preflight_output_collisions(` to `preflight_output_collisions(` closed a real hole (an
`use`-imported direct call would have passed both assertions), and I confirmed the count is still
exactly 1.

---

## Fixes applied

**None by me.** Three of the fixes I had queued were already present in the working tree when I went to
apply them, from concurrent work by another reviewer:

| Finding I had queued | State in the working tree |
|---|---|
| `alias.exists()` misses a dangling symlink → raw `EEXIST` instead of the friendly bail | Fixed: `std::fs::symlink_metadata(&alias).is_ok()` (`src/main.rs:130`) |
| The 18-line `run`-parameter doc block was absorbed onto the new `pub enum SeqFormat` | Fixed: split back onto `run` and `resolved_threads`, with `sequence_format` documented |
| `"a stranded symlink would trip the next run's output-collision pre-flight"` is false | Fixed in both `src/main.rs` and the test message |

I verified all three rather than re-deriving them: tests still 7 + 57 pass, clippy clean on a forced
recheck, and the dangling-symlink probe now produces the intended message
(`Cannot stage the raw-FastQC input alias './a_raw.fastq.gz': a file is already there`) where before the
fix it produced `Failed to create the raw-FastQC input alias … File exists (os error 17)`.

I made no edits of my own to avoid conflicting with uncommitted work in a shared tree. Everything above
is a recommendation.

---

## Recommendations, prioritized

**Critical** — none. No finding corrupts sequence data or the trimmed output.

**High**
1. **H1** — give the alias and the two report paths a directory that matches the mode's primary output,
   *and* put them in the pre-flight `planned` list. Until then, `--fastqc_raw` silently discards raw
   reports on multi-sample runs without `--output_dir`, exit 0.
2. **H2** — add `--fastqc_raw` to `fastqc_active` and add at least one SE and one PE tripwire case,
   including a same-basename pair. This is the guard that makes H1 un-reintroducible.
3. **H3** — resolve the `--fastqc_args` × specialty-mode refusal, or retract the documented support.
   The current error names `--fastqc` on a command line that does not contain it.

**Medium**
4. **M1** — covered by H1's pre-flight fix; verify the `s_raw.fastq.gz` + `s.fastq.gz` case afterwards.
5. **M2** — replace the two-step cleanup with a drop guard so a panic in fastqc-rust cannot strand the
   alias.
6. **M3** — drop the `symlink_metadata` pre-check and map `EEXIST` from `symlink()` to the friendly
   message; one syscall, no window, no duplicated logic.
7. **M4** — decide and document whether `--passthrough` inputs are in scope for the raw pass.
8. **M5** — move `raw_alias_path` / `run_pretrim_fastqc` / `preflight_then_raw_qc` into `src/lib.rs` so
   the alias-naming table can be unit-tested; fix the `else`-branch comment, which describes a case
   that reaches the other branch.

**Low**
9. **L1** — `.bgz`/`.bgzf` input yields `a_raw.fq.bgz_fastqc.zip`. Content is correct and the MultiQC
   sample name is still distinct; cosmetic only.
10. **L3** — `symlink_metadata().is_err()` in `the_staged_alias_does_not_survive_the_run`.
11. **L4** — add the PR link to the CHANGELOG entry when one exists (author acknowledged; every other
    entry in the file carries one).
12. Add a uBAM `--fastqc_raw` test and a specialty-mode `--fastqc_raw` test, so the `SeqFormat`
    mechanism and the specialty wiring are pinned by something other than this review.
13. PLAN §13 now drifts from the tree on two details (the inventory assertion string, and the
    "output-collision pre-flight" rationale for removing the alias). Worth a line when the PR is opened.

---

## Verified vs. inferred

**Verified at the source or by running it:**
13 wrapper call sites and exactly one surviving `naming::preflight_output_collisions` call; every
`return`ing dispatch arm in `main()` reaches a wrapper; all four `run_clump_only` branches and the
`--paired` FASTQ clump branch (via `run_specialty_paired`) are covered; `#[track_caller]` chains
through the wrapper (`every_site_is_reached_by_some_case` passes); `SeqFormat::Bam` overrides
extension sniffing (masquerade probe); fastqc-rust detects gzip by magic, not extension (extensionless
probe); all 12 pre-existing `fastqc::run` sites take `None` and none changed behaviour — including the
three `.bam` ones, which now pass `None` and so keep their extension-based dispatch, which I confirmed
still works via `--clump_only --output-format ubam`; `args_set_outdir` agrees with
`apply_fastqc_args` on every tokenisation I could construct (`-o=dir` warns-and-ignores in both,
trailing `-o` is a no-op in both, `-t -o 4` sets outdir in both); H1, H3, M1, M4, L1 reproduced with
the commands shown; clippy/fmt clean on a forced recheck; 7 + 57 + 557 tests pass.

**Inferred, not reproduced:** M2 (no input found that panics fastqc-rust 1.1.1); M3 (concurrency
window reasoned from the code, not raced).

**Could not check:** the "896 tests pass" total — I ran 621 across three binaries, not the whole suite;
Windows portability (`std::os::unix::fs::symlink`, acknowledged in §13 "Not done", no CI job);
MultiQC's actual parsing behaviour — I verified the `Filename` field the design depends on, not
MultiQC itself.

**Wrong in the committed code:** the comment "a stranded symlink would trip the next run's
output-collision pre-flight". A stranded alias is neither a planned output nor an input, so the
pre-flight never sees it; what it trips is `run_pretrim_fastqc`'s own occupancy check. Already
corrected in the working tree.
