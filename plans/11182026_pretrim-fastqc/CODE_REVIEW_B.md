# Code review — `--fastqc_raw` (Reviewer B)

**Target:** `git diff origin/dev..feat/fastqc-raw` (`61ae72c` → `cf3bc84`)
**Date:** 2026-11-18
**Reviewer:** B of two independent reviewers (no shared state with A)

## Summary

The mechanism is sound and the hard part is right. The staged-symlink trick does what the plan
claims — I verified the `Filename` field, the `sequence_format` override, and the `#[track_caller]`
attribution chain against the dependency source and against running binaries, not against the plan.
The 13-site wrapper is a real structural improvement over 13 `if cli.fastqc_raw` blocks, and
`raw_alias_path` survived every name I threw at it.

Three defects are worth blocking on, and the first two share one root cause:

1. **A raw report is silently overwritten** when `--output_dir` is unset and two inputs share a
   basename across directories. Reproduced. The surviving report's `Filename` field is identical to
   the one it replaced, so it cannot be attributed. PLAN §3.8 states this case is "now caught by
   §3.6's pre-flight" — it is not.
2. **Raw reports land in the current working directory**, not beside the input, when `--output_dir`
   is unset — unlike the trimmed output, the trimming report and the post-trim FastQC report, which
   all go to `input.parent()`. Reproduced. Three user-facing strings say otherwise, including the
   new `-o` refusal message.
3. **`--fastqc_raw --fastqc_args X` is refused on all four specialty modes**, with an error naming
   `--fastqc`, which the user did not pass. Reproduced. This falsifies the `--help` text's
   "Accepted with the specialty modes that reject --fastqc" and PLAN §3.9.

Beyond those: the new output surface is entirely outside the exact-set pre-flight tripwire (the test
that exists because this defect family shipped five times), and one of the seven new integration
assertions cannot fail on the mistake it names.

I applied four low-risk fixes (two guard/test strengthenings, two comment corrections). `cargo fmt`,
`cargo clippy --all-targets --release -- -D warnings` and the full suite are clean after them.

---

## What I verified vs. what I inferred

### Verified by execution

| Claim | Result |
| --- | --- |
| "896 tests pass (883 before, +13 new)" | **True.** `cargo test` → 896 passed, exit 0. |
| `cargo fmt --all -- --check` / `clippy -D warnings` clean | **True**, re-confirmed after my edits. |
| `#[track_caller]` chains through the wrapper to the dispatch arm | **True.** `TRIM_GALORE_DUMP_PLANNED_OUTPUTS=1` reports `src/main.rs:1312` for SE and `src/main.rs:1227` for PE — the arms, not the wrapper at `:183`. |
| The alias is removed on the error path | **True.** `truncated.fq.gz` run leaves no `*_raw.fq.gz`. |
| An input is never overwritten by the alias | **True.** Inputs `a.fastq.gz` + `a_raw.fastq.gz` in one CWD → refusal, input md5 unchanged, nothing written. `symlink(2)` is the real atomic guard; the `exists()` check only buys the better message. |
| V6 — specialty modes do write raw reports | **True** for `--clock` (checked); code path is shared with `--implicon`/`--hardtrim`. |
| uBAM named `reads.dat` is analysed as BAM | **True.** Report produced; format is declared, not sniffed. |
| Input with no extension at all | **Works.** `noext` → `noext_raw_fastqc.{html,zip}`. |

### Verified by reading source (dependency or repo)

| Claim | Where | Result |
| --- | --- | --- |
| `sequence_format` overrides extension dispatch | `fastqc-rust-1.1.1/src/sequence/bam.rs:402-431` | **True.** `"fastq"` and `"bam"` are both accepted; an unrecognised string is an `InvalidInput` error, so the enum cannot produce one. |
| `args_set_outdir` mirrors `apply_fastqc_args` | `src/fastqc.rs` vs `fastqc.rs:158-166` | **True.** I enumerated the interleavings (`-t -o /tmp`, `-o -o`, `--outdir -o`, trailing `-o`, `-o` as a value) — the two agree on "is outdir set" in all of them. |
| `raw_alias_path` cannot panic on a char boundary | `src/main.rs:86-99`, `io.rs:816-842` | **True.** `stem` is a prefix of `name` on both branches (suffix-stripping is byte-exact on ASCII suffixes; `file_stem` of the gz-stripped prefix is itself a prefix), so `&name[stem.len()..]` is always at a boundary. |
| 13 pre-flight sites, each individually exercised | `main.rs`, tripwire | **True.** 13 call sites; exactly one `preflight_output_collisions(` in `main.rs`. I traced every dispatch arm (`main`, `run_clump_only`, `run_ubam_output`, `run_specialty_paired`, `run_paired_trim`, `run_single_end_trim`) — none returns without passing through one. |

### Wrong claims found

| Claim | Where | Reality |
| --- | --- | --- |
| "a stranded symlink would trip the next run's output-collision pre-flight" | `main.rs:155-156` comment, repeated in `integration_fastqc_raw.rs:91` | **False.** `preflight_output_collisions` does no filesystem probing at all (`io.rs:180-250`) — it compares planned paths against each other and against the input list. A stranded alias trips the function's *own* staging check, and only if its target still resolves. The real consequence is that the alias is a FASTQ-named file a later glob picks up as a sample. **Fixed.** |
| "The extension is kept so a gzipped FASTQ is still recognised as one" | `main.rs:85` | **False as stated.** fastqc-rust detects compression from the magic bytes (`sequence/fastq.rs:244-256`), not the extension — and `sequence_format` is always `Some(..)` on this path, so the extension affects nothing but the report's name. Harmless, but the rationale does not hold. |
| "Two inputs sharing a basename collide — now caught by §3.6's pre-flight" | PLAN §3.8 | **False** when `--output_dir` is unset. See Critical-1. |
| "the raw reports are written beside the trimmed output" | `cli.rs` refusal message | **False** when `--output_dir` is unset. See High-1. |
| "`--fastqc_raw` is deliberately not subject to [the #421 refusals]" | PLAN §3.9, and the `--fastqc_raw` help text | **False** the moment `--fastqc_args` is present. See High-2. |

### Could not verify

- **MultiQC's `fn_clean_exts` does not strip `_raw`.** No MultiQC available offline. The design rests
  on it (Assumption 7 covers the `Filename`-field half, which I did verify in fastqc-rust). `_raw` is
  not a suffix I recognise from MultiQC's default list, and `_trimmed`/`trimmed` are, so the raw and
  trimmed names should survive cleaning as distinct — but this is inference, and it is the one
  assumption the whole feature rests on. Worth one real MultiQC run before the nf-core bump.
- **Windows portability.** `std::os::unix::fs::symlink` will not compile for a Windows target
  (author already records this under "Not done"). No CI job covers one, so I could not test it.
- **Behaviour if the process is killed mid-run.** Could not time a signal precisely; reasoned from
  the code that the alias is stranded (no `Drop` guard) — see Medium-3.

---

## Issues by area

### Logic

**Critical-1 — a raw report is silently overwritten across directories (no `--output_dir`).**

`run_pretrim_fastqc` (`src/main.rs:119-120`) derives one directory for the whole run:

```rust
let dir = output_dir.unwrap_or_else(|| Path::new("."));
```

Every alias therefore lands in the CWD, and so does every report. Two inputs with the same basename
in different directories produce the same alias name — sequentially, so the `symlink_metadata` guard
never sees a conflict — and the same report name, which the second run overwrites.

Reproduced:

```
$ ls d1 d2
d1/a.fastq.gz   (illumina_10K, 150 bp)
d2/a.fastq.gz   (BS-seq_10K,    65 bp)
$ cd cwd && trim_galore --fastqc_raw ../d1/a.fastq.gz ../d2/a.fastq.gz
$ ls cwd
a_raw_fastqc.html  a_raw_fastqc.zip        # one report, no warning, exit 0
$ unzip -p a_raw_fastqc.zip a_raw_fastqc/fastqc_data.txt | head -4
Filename        a_raw.fastq.gz            # identical to what d1's would have said
Sequence length 65                        # d1's report is gone
```

The trimmed outputs are both written correctly (`d1/a_trimmed.fq.gz`, `d2/a_trimmed.fq.gz`), so the
pre-flight has nothing to refuse — the two raw reports are the only colliding paths, and they are not
candidates. The surviving report cannot be attributed to either input, because `Filename` is derived
from the alias name, which is the same for both.

This is the #216/#383/#391/#409 shape: a writer whose path never entered the candidate list. It is
narrower than #409 (no read data is destroyed), but a QC report silently attributed to the wrong
sample is exactly what the pre-flight family exists to prevent, and `--output_dir`-less multi-input
runs are ordinary CLI usage.

With `--output_dir` set the hole closes transitively: any two inputs whose raw stems collide also
collide on their trimmed output in that one directory, so the pre-flight refuses first. I checked
`.fastq.gz`/`.fq.gz`, `.FASTQ.GZ`/`.fastq.gz`, `--basename`, and the hardtrim namers — all refuse.

**Recommended fix** (also fixes High-1, one line, and it is the fallback every sibling namer already
uses — `io.rs:302-305`):

```rust
// in the per-input loop, replacing the hoisted `dir`
let dir = output_dir.unwrap_or_else(|| input.parent().unwrap_or(Path::new(".")));
let alias = raw_alias_path(input, dir);
```

fastqc-rust then falls back to the alias's parent for the report (`runner.rs:299-308`), so the raw
reports land beside the input, next to the trimmed output and the post-trim report. The alternative —
adding the raw report paths to all 13 candidate lists — is a much larger change for the same effect,
and would also have to teach the tripwire about them.

This contradicts PLAN §3.8's "`--output_dir` unset → reports and the transient symlink in the working
directory", which is why I am recommending rather than applying it: the author chose the current
behaviour deliberately, and the fix changes a documented output location.

**High-1 — raw reports go to the CWD while everything else goes beside the input.**

Same root cause. Reproduced:

```
$ cd /tmp/x/cwd && trim_galore --fastqc_raw --fastqc /tmp/x/data/sample.fastq.gz
/tmp/x/cwd:   sample_raw_fastqc.html  sample_raw_fastqc.zip
/tmp/x/data:  sample_trimmed.fq.gz  sample_trimmed_fastqc.{html,zip}
              sample.fastq.gz_trimming_report.{txt,json}
```

Three strings disagree with this: the new refusal in `cli.rs` ("the raw reports are written beside the
trimmed output"), `README.md` ("writing `*_raw_fastqc.html`/`.zip` alongside"), and `fastqc.rs`'s
module doc, which states the house convention ("`None` falls back to the parent directory of the
input file, which is where the trimmed output already lives"). A read-only CWD now fails a run whose
output would otherwise have succeeded.

Fix as per Critical-1. If the CWD behaviour is kept deliberately, the refusal message and the README
line both need correcting.

**High-2 — `--fastqc_args` makes `--fastqc_raw` unusable in the specialty modes.**

`fastqc_requested()` is `self.fastqc || self.fastqc_args.is_some()` (`cli.rs:617`), and the #421
refusals key on it, so:

```
$ trim_galore --hardtrim5 20 --fastqc_raw --fastqc_args "--nogroup" x.fastq.gz
Error: --fastqc is not supported with --hardtrim5, which writes no QC report;
       run FastQC separately on the hardtrimmed output, or drop --fastqc
```

The user passed no `--fastqc`. The flag's own help text promises the opposite ("Accepted with the
specialty modes that reject --fastqc"), and there is no way to pass engine args to the raw pass in
those modes at all. PLAN §3.1's last row documents the *trimmed-report* consequence of the implication
and Open-1 leaves it open; neither notes that it makes four of the five modes the flag was widened for
unreachable.

`fastqc_still_refused_in_specialty_modes_alongside_fastqc_raw` passes `--fastqc` explicitly, so it
does not cover this.

Minimal fix, if the implication is to stay: narrow the #421 gate to `self.fastqc || (self.fastqc_args.is_some() && !self.fastqc_raw)`
— i.e. `--fastqc_args` implies `--fastqc` except where the only QC asked for is the raw pass. That
changes `fastqc_requested()`'s meaning at twelve sites, so it needs its own thought; the cheaper
alternative is `--fastqc_raw_args`, as Open-1 already lists. Either way the help text should not
promise something the common combination breaks.

**Medium-1 — a later input's staging failure leaves earlier inputs' reports on disk.**

§3.7's up-front loop sanity-checks every input before analysing any, precisely so that "R1's report is
written and only then does R2 turn out to be unusable" cannot happen. The alias-staging check is not
in that loop — it sits inside the analysis loop — so the same shape still occurs for a staging
conflict rather than a bad record:

```
$ ls .
a.fastq.gz  a_raw.fastq.gz   # a_raw is a real sample, not a leftover
$ trim_galore --fastqc_raw a_raw.fastq.gz a.fastq.gz
# a_raw_raw_fastqc.* written, then:
Error: Cannot stage the raw-FastQC input alias './a_raw.fastq.gz': a file is already there.
```

(With the inputs in the other order nothing is written, which is what I tested first — the ordering
dependence is the tell.) Hoisting the alias computation and the `symlink_metadata` check into the
existing validation loop closes it in ~6 lines and matches the pattern §3.7 established. It cannot
close the concurrent-process race, but that race is not the case this guard is for.

**Medium-2 — `--passthrough` carriers get no raw report, but are still scanned for `--phred64`.**

Open-3 settles the "inputs only" question, and I am not reopening it. The asymmetry worth noting is
that `--fastqc` *does* QC the carrier's trimmed output (`main.rs:2135-2142`), so `--fastqc_raw --fastqc`
on a passthrough run produces a before/after pair for the reads and an after-only report for the
carrier. Meanwhile `main.rs:507` now widens the `--phred64` pre-flight to iterate
`cli.input.chain(cli.passthrough)` on `cli.fastqc_raw`, scanning a file the raw pass will never read.
Harmless (strictly stricter), but the two halves disagree about whether the carrier is in scope.

**Medium-3 — no panic safety on the staged alias.**

Cleanup is a manual `let _ = remove_file(&alias)` between `fastqc::run` and `result?`. That covers
`Err`, which is the case that matters most and which the truncation test exercises. It does not cover
a panic inside fastqc-rust (malformed input is the plausible trigger) or a signal. A `struct
AliasGuard(PathBuf)` with `impl Drop` is ~8 lines and covers the unwind case, converting a stranded
FASTQ-named symlink into nothing. Signals still strand it; see Low-5.

**Low-1 — the run-once `AtomicBool` guards nothing today, and hides a future miss.**

I traced all 13 sites: every dispatch arm returns after its pre-flight, so no process reaches two of
them. The guard is currently unreachable in its "skip" branch. Its failure mode is one-directional and
silent: if a future arm ever pre-flights per-pair or per-chunk with different inputs, the raw pass
covers the first call only, with no diagnostic. If it is kept as belt-and-braces, keying it on the
input list rather than a bare bool would make a second legitimate call work instead of vanish.

(It is also a function-local `static`, so it is per-process, not per-run. That is fine for a binary —
integration tests spawn subprocesses — but it does mean `main` cannot be driven twice in one process.)

**Low-2 — "a rejected run creates nothing" is now narrower than `main.rs:478-482` reads.**

Verified: a failure between the pre-flight and the first write still leaves the raw reports.

```
$ trim_galore --fastqc_raw -a file:/nonexistent/adapters.fa -o out/ x.fastq.gz
Error: ... No such file or directory
$ ls out/
illumina_10K_raw_fastqc.html  illumina_10K_raw_fastqc.zip
```

The reports are complete and correct descriptions of the input, so this is arguably fine — but the
invariant as stated in the comment the plan quotes (§3.6) is now "a run rejected *by the pre-flight* creates
nothing". Worth one word in that comment.

### Errors / robustness

**Low-3 — `exists()` missed a dangling symlink.** `Path::exists()` follows links, so the curated
"a file is already there. Remove it" message did not fire for a dangling alias — precisely the state
a killed previous run leaves once the input moves. It fell through to the generic `symlink(2)` EEXIST
wrapper instead. No safety consequence (the input is never at risk either way), message quality only.
**Fixed.**

**Low-4 — report names leak odd extensions, breaking the documented `*_raw_fastqc.*` contract.**
fastqc-rust's `strip_extensions` (`runner.rs:739-750`) is case-**sensitive** and does not know `.bgz`
or `.bgzf`, whereas TG's `strip_fastq_extensions` is case-insensitive and does. Verified:

```
UP.FASTQ.GZ  →  UP_raw.FASTQ.GZ_fastqc.html     (MultiQC sample: UP_raw.FASTQ.GZ)
reads.dat    →  reads_raw.dat_fastqc.html
```

The raw/trimmed pair stays distinct, so the feature still works; the names are just not what README,
CHANGELOG and the wiki promise. Lower-casing the extension in `raw_alias_path` would fix the uppercase
case; `.bgz` would need a normalisation decision. Since the format is now declared to the engine, the
alias extension is free to be normalised.

**Low-5 — concurrent identical runs into one `--output_dir` now hard-fail at staging.** Verified: two
simultaneous `trim_galore --fastqc_raw -o out/ same.fastq.gz` → one succeeds, the other refuses at the
alias. Such a pair would clobber each other's trimmed output anyway, so refusing is defensible; the
realistic concurrency shape (different inputs, shared `--output_dir`) has distinct aliases and is
unaffected. Noting it because the error text ("Remove it, or rename the input") suggests a stale file
rather than a concurrent run, and a Drop-guard (Medium-3) plus a hoisted check (Medium-1) do not help.

**Low-6 — the alias is staged in the user's CWD when `--output_dir` is unset**, which may be a git
working tree or a read-only directory. Subsumed by Critical-1's fix.

### Structure / style

**Low-7 — a pre-existing doc-comment mis-attachment was relocated and compounded.** On `dev`, `run`'s
parameter block was already attached to `resolved_threads`. Inserting `SeqFormat` between them moved
the whole block onto the new enum, so `rustdoc` for `SeqFormat` now opens with "Run FastQC on a single
trimmed output file" and `run` has no documentation at all. **Fixed** by relocating each block to its
own item (plus the missing `sequence_format` line, and "trimmed output file" → "file", since `run` is
no longer called only on trimmed output).

The module header (`fastqc.rs:10-12`, "Callers in `src/main.rs` invoke it with a trimmed output path")
is stale for the same reason. Left alone — prose, and the author may want to word it.

**Low-8 — redundant `ensure_output_dir`.** `run_pretrim_fastqc` (`main.rs:119`) calls it, but `main` already did
(`main.rs:673`), before every one of the 13 pre-flight sites. §5 step 6c's rationale ("fastqc-rust
does not create its `-o` directory") is about a directory that already exists by then. Harmless and
cheap; flagging only because the comment implies it is load-bearing.

**Low-9 — CHANGELOG entry carries no issue/PR link.** Every neighbouring entry has one; the author
records this under "Not done". Two sentences is in line with the surrounding entries, so I am not
raising the length against CLAUDE.md's one-sentence rule — but the link is a local convention.

**Low-10 — `docs/.../guide/outputs.md` is not updated** (author-recorded). `*_raw_fastqc.*` is a new
user-facing output surface and the brief says an nf-core module bump will consume it; the outputs page
is where a pipeline author looks for the glob.

### Tests

**Medium-4 — the exact-set tripwire has no `--fastqc_raw` case, and would go red if one were added.**
`fastqc_active` (`tests/integration_preflight_tripwire.rs:233`) matches only `--fastqc` and
`--fastqc_args`, so a `--fastqc_raw` case's `*_raw_fastqc.{html,zip}` would be classified as stray
output and fail `compare`. The result is that the one new output surface this branch adds is outside
the audit that exists because the same defect shipped five times — and, as Critical-1 shows, that
surface does have a hole. Adding `--fastqc_raw` to `fastqc_active` and one SE case would have put
Critical-1 within range of the test (a case with two same-basename inputs and no `--output_dir`).

**Medium-5 — one new assertion cannot fail on the mistake it names.**
`raw_and_trimmed_reports_coexist_with_distinct_sample_names`:

```rust
assert_ne!(
    basic_stat(&raw, "illumina_10K_raw", "Filename"),
    basic_stat(&trimmed, "illumina_10K_trimmed", "Filename"),
);
```

commented "Distinct in the field MultiQC reads — the point of the whole feature". Remove the alias
entirely and the two values become `illumina_10K.fastq.gz` and `illumina_10K_trimmed.fq.gz` — still
unequal, so the `assert_ne!` still passes. Only the preceding `raw.exists()` catches the regression,
and that is already covered by test 1. The failure the design prevents is *post-cleaning* collapse
(`illumina_10K` vs `illumina_10K`), which nothing here models. Asserting on the names after stripping
`_trimmed` and the FASTQ/gz extensions — MultiQC's own cleaning, approximated — would make the
assertion load-bearing.

**On the other six:** `single_end_raw_report_names_the_sample_raw`, `the_staged_alias_does_not_survive_the_run`,
`paired_raw_reports_are_named_per_mate`, `a_collision_refused_run_writes_no_raw_report` and
`a_bad_second_input_is_refused_before_any_report_is_written` all fail on the mistake they describe — I
checked each by asking what breaks if the mechanism is removed. The first's `Sequence length 150`
assertion is a genuinely good second lock. `mid_file_truncation_surfaces_from_the_qc_pass` passes
whether the failure comes from the QC pass or from trimming, so its name claims more than it checks;
its assertions (run fails, R1 report exists) are both correct and both informative, so this is a
naming nit, not a dead test.

**On the inventory tripwire:** `src.matches("preflight_then_raw_qc(").count() - 1` is correct (14
occurrences, one of which is the definition) and `every_site_is_reached_by_some_case` really does
still prove per-arm attribution — I confirmed the dump prints `:1227`/`:1312`, not `:183`. The new
"exactly one direct call" assertion was scoped to `naming::preflight_output_collisions(`, which a
future `io::`- or `crate::io::`-qualified call would slip past while adding nothing to the other
count. **Fixed** by dropping the module prefix (still exactly 1).

### Efficiency

Nothing to raise. The extra pass is inherent to the feature and documented; re-detecting the format
per input (deviation 4) is one block-decompress against a full-file QC pass and avoids six signature
changes — the right trade. Sequential execution keeps `resolved_threads` and the clumpify memory model
correct unchanged, which I confirmed: with `--fastqc_raw` alone `fastqc_requested()` is false, so the
model reserves nothing for FastQC — correct, because the raw pass has finished before clumping starts.
The duplicate `sanity_check_any` per input (once here, once in the driver) reads one record. The
double `to_string_lossy` in `raw_alias_path` is one allocation per input.

---

## Fixes applied

All four are in the working tree, unstaged. `cargo fmt --all -- --check` clean,
`cargo clippy --all-targets --release -- -D warnings` clean, `integration_fastqc_raw` (7) and
`integration_preflight_tripwire` (57) green.

1. **`src/main.rs:129-130`** — `alias.exists()` → `std::fs::symlink_metadata(&alias).is_ok()`, so the
   curated staging error also fires for a dangling symlink. Verified before/after against a planted
   dangling link. (Low-3)
2. **`src/main.rs:155-156`** — corrected the cleanup comment, which asserted a pre-flight behaviour
   that does not exist. Now states the actual consequence of a stranded alias, in one line, present
   tense.
3. **`tests/integration_fastqc_raw.rs:91`** — same wrong claim in an assertion message; corrected to
   match.
4. **`tests/integration_preflight_tripwire.rs:1230`** — `"naming::preflight_output_collisions("` →
   `"preflight_output_collisions("`, so a differently-qualified direct call is also caught. Count is
   still exactly 1.
5. **`src/fastqc.rs`** — relocated the hijacked doc block: `run`'s parameter list back onto `run` (with
   the missing `sequence_format` line), the two thread lines back onto `resolved_threads`, leaving
   `SeqFormat`'s own doc intact. Comment-only. (Low-7)

I did **not** touch anything behavioural beyond fix 1 — the three blocking defects all involve a
documented design decision or a choice between remedies.

---

## Recommendations, prioritised

### Critical

1. **Close the silent raw-report overwrite** (Critical-1). Derive the alias directory per input from
   `input.parent()` when `--output_dir` is unset. One line, and it fixes High-1 at the same time.
   Requires updating PLAN §3.8, which currently documents the CWD behaviour and wrongly claims the
   collision case is covered. If the CWD behaviour is to be kept instead, the raw report paths must
   enter the pre-flight candidate lists, or the run must refuse when two inputs share a basename.

### High

2. **Make the raw reports land where every other artefact lands** (High-1), or correct the refusal
   message, the README line and PLAN §3.2/§3.8 to say "the working directory". Shipping the current
   combination means the error message tells users something false about where their files went.
3. **Decide `--fastqc_args` + `--fastqc_raw` in the specialty modes** (High-2). Either narrow the #421
   gate, add `--fastqc_raw_args`, or — at minimum, if it stays — drop "Accepted with the specialty
   modes that reject --fastqc" from the help text and say what actually happens. A refusal naming a
   flag the user did not pass is the worst of the three.

### Medium

4. **Add a `--fastqc_raw` case to the exact-set tripwire** and teach `fastqc_active` the flag
   (Medium-4). The case that would have caught Critical-1 is two same-basename inputs with no
   `--output_dir`.
5. **Hoist the alias-staging check into §3.7's validation loop** (Medium-1), so a staging conflict on
   input 3 cannot leave inputs 1 and 2 reported on.
6. **Make `raw_and_trimmed_reports_coexist_with_distinct_sample_names` able to fail** (Medium-5) by
   asserting on MultiQC-cleaned names rather than raw `Filename` values.
7. **Add a `Drop` guard for the alias** (Medium-3) so a panic in the engine does not strand a
   FASTQ-named symlink in the output directory.
8. **Run one real MultiQC** over a `--fastqc_raw --fastqc` output set before the nf-core bump. The
   `Filename`-reaches-the-report half is verified; the `_raw`-survives-cleaning half is not, and it is
   the half the user-visible promise rests on.

### Low

9. Resolve the carrier asymmetry (Medium-2) — either QC the `--passthrough` carrier too, or stop
   scanning it in the `--phred64` pre-flight on `fastqc_raw` alone.
10. Normalise the alias extension (Low-4) so `.FASTQ.GZ` and `.bgz` inputs still produce
    `*_raw_fastqc.*`.
11. Reconsider or re-key the run-once `AtomicBool` (Low-1) — today it is unreachable; tomorrow it is a
    silent skip.
12. Narrow the "a rejected run creates nothing" comment (Low-2) to "rejected by the pre-flight".
13. Drop the redundant `ensure_output_dir` (Low-8), or note that `main` already guarantees it.
14. Correct `raw_alias_path`'s doc — the extension is kept for the report's name, not for gzip
    recognition, which is magic-byte based.
15. Refresh `fastqc.rs`'s module header ("Callers … invoke it with a trimmed output path") and add
    `--fastqc_raw` to `guide/outputs.md`; add the issue/PR link to the CHANGELOG entry when one exists.
