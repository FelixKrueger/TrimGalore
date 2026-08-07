# Code Review B — #383 output-collision pre-flight

**Reviewer:** B (independent; no coordination with Reviewer A)
**Target:** uncommitted working tree on `fix/383-output-collision-preflight` off `dev` @ `72624c7`
**Plan:** `plans/08062026_se-output-collision-preflight/PLAN.md` (v2, incl. §13)

**No fixes were applied — not even one-liners.** The change is uncommitted in a working tree
shared with a concurrent reviewer, so an edit from me would corrupt the diff we are both reading.
Everything below is a recommendation.

**Two process notes, both material to reading this report:**

1. **The tree moved during the review.** I started against `6 files changed, +503/−99`
   (`io.rs +216`, no `demux.rs`). It is now `7 files changed, +595/−111` (`io.rs +289`,
   `demux.rs +31/−11`), and `PLAN.md` gained a **D6** entry. The added work is
   `demux::demux_base_name` plus a second A2 unit test
   (`distinct_primary_outputs_imply_distinct_secondary_outputs`). Every finding below is
   re-checked against the **current** source; where a finding was overtaken I say so.
   `target/release/trim_galore` (built 2026-08-06 23:34) is **stale** relative to the current
   tree — see the note under Finding 1 for why that does not affect the reproductions.
2. **A different file was already at this path** (64 lines, "In progress", findings B-1/B-2)
   when I came to write. I have preserved it verbatim at
   `…/scratchpad/creview-b/CODE_REVIEW_B.preexisting.md` and reproduced its two findings,
   independently verified, in Appendix A. If another agent is also writing here, that content
   is not lost.

Verification method: constructing invocations the new tests and CI steps do not cover, and
running the binary. Scratch scripts `t1.sh`–`t9.sh` in
`/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/d0ee0171-fc72-476b-bdba-c46c90a5bacd/scratchpad/creview-b/`.
Everything marked *reproduced* was executed.

---

## 1. Summary

The fix is sound on every path I exercised. All eleven call sites pass a candidate list that
matches what the corresponding writer actually opens (A9 holds — checked call-by-call, not
assumed), the four converted hand-rolled sites are behaviour-preserving apart from the two
intended defect fixes, and both new defects are genuinely closed. D1's predicate is right, and
no fifth paired-shaped mode needs excluding. The two dispatch paths that still have no
pre-flight cannot alias their own input — verified, not assumed.

One finding is material: **defect 2 is closed for primary outputs only, and `--demux` still
destroys a named input at exit 0** — reproduced. That also makes a CHANGELOG sentence
inaccurate, which is why it is High rather than a scoped residual.

The rest is diagnostics quality (one error message states something false), one coverage gap on
the widest hardtrim case, and a doc-comment re-parenting introduced by the mid-review D6 change.

Counts: **0 Critical, 1 High, 4 Medium, 9 Low.**

---

## 2. Verified correct (so it need not be re-checked downstream)

**A9 — candidate expression vs. writer, all eleven sites.** Each candidate builder is
character-identical to the namer call inside the writer it guards:

| Site (`main.rs`) | Candidate | Writer's own call |
|---|---|---|
| `:360` `--hardtrim5` | `planned_hardtrim_outputs(…, Five, output_dir, gzip)` | `specialty.rs:45` / `:139` (BAM) |
| `:391` `--hardtrim3` | same, `Three` | `specialty.rs:88` / `:191` |
| `:513` `--clump_only` SE FASTQ | `clumped_output_name(input, output_dir, basename, gzip)` | `clump_only.rs:278` |
| `:548` `--clump_only` PE uBAM N=1 | `clumped_paired_bam_output_name(input[0], None, …)` | `clump_only.rs:922` |
| `:585` `--clump_only` PE uBAM N≥2 | `clumped_paired_bam_output_name(chunk[0], Some(chunk[1]), …)` | `clump_only.rs:969` |
| `:633` `--clump_only` SE uBAM | `clumped_bam_output_name(input, output_dir, basename)` | `clump_only.rs:754` |
| `:734` `--paired` trim FASTQ | `paired_end_output_names` + unpaired + passthrough | `main.rs:1293`, `:1311`, `:1305` |
| `:799` SE trim FASTQ | `single_end_output_name(input, output_dir, basename, gzip)` | `main.rs:1103` |
| `:1825` `--paired` trim uBAM | `paired_bam_output_name(chunk[0], chunk[1], …)` | `run_ubam_output_paired_two_files` |
| `:1869` SE trim uBAM | `single_end_bam_output_name(input, output_dir, basename)` | `main.rs:1893` |
| `:2430` `run_specialty_paired` | caller's `output_names` closure | `specialty.rs:254`/`:367`, `clump_only.rs:416` |

`--rename` does not participate in naming (`specialty.rs:60`, `:104` only mutate read IDs), and
`gzip` is one variable on both sides, so **A4 holds in both argument orders** — *reproduced*
(`t8.sh`): `sample.fastq sample.fq.gz` → rejects on `sample_trimmed.fq`;
`sample.fq.gz sample.fastq` → rejects on `sample_trimmed.fq.gz`.

**The four converted sites, probed individually** (`t2.sh`, `t3.sh`, `t4.sh`, `t6.sh`):

- `--clump_only` SE colliding stems → exit 1, and the message now says `` `--output_dir` ``
  (`grep -rn -- "--output-dir" src/ docs/` is empty, so §2.5's correction landed).
- `--clump_only` SE distinct stems → exit 0, both `_clumped.fq` + both clumping reports.
- `--clock`, two colliding pairs → exit 1 with an empty `-o` dir; one pair → exit 0 with both
  `_UMI_` outputs.
- `--paired` uBAM, two pairs colliding into one `-o` → exit 1, nothing written; two distinct
  pairs → exit 0, both `_val.bam`.
- `--paired --retain_unpaired` where an **unpaired** candidate aliases an input → exit 1 with the
  alias wording, prior file byte-intact — so the `--retain_unpaired` extras survived the
  conversion.
- `--paired` with mixed `./` spelling across two pairs → exit 1: the paired path gained
  defect-1 protection as well as defect-2.
- The alias check is **order-independent** — `trim_galore s_trimmed.fq s.fastq` (alias listed
  first) rejects identically, with the prior file intact.

**`guarded_inputs()` is the right input set.** `--passthrough` is `--paired`-only
(`cli.rs:885`) and single-pair-only (`:890`), so on every other guarded path it contributes
nothing. *Reproduced*: the full `--passthrough` recipe (`BS-seq_10K_R1/R2 + I1`) still exits 0
with all seven expected outputs, i.e. adding the passthrough file to the *input* set produced no
false positive. §3.2's argument holds for the `--demux` **barcode file** — but not for demux
**outputs**, see Finding 1.

**D1's predicate.** `!self.paired && !self.clock && self.implicon.is_none()` is exactly right:
`--paired`, `--clock`, `--implicon` and `--clump_only --paired` (which sets `paired`) are the
only paired-shaped modes, and `--hardtrim5/3` never call `validate_paired_input`
(`cli.rs:617`, `:955`, `:958` are its only call sites) and loop per file, so including them is
correct. *Reproduced*: `--clock R1 R1` → "Read 1 and Read 2 appear to be the same file";
`--implicon=8` with a duplicated pair → "is a duplicate of pair 1"; a duplicated positional
under `--hardtrim5` → "was given more than once".

**The two unguarded dispatch paths cannot alias their input.** `run_paired_ubam_single_file`
(`main.rs:1610`) and `run_ubam_output_paired_single_file` (`:2183`) both derive the stem from
`Path::file_stem()` and append `_val`/`_val_1`, and neither reads `cli.basename`, so output ≠
input for every input name. *Reproduced*: `--paired --output-format ubam x_val.bam` →
`Output: x_val_val.bam`, exit 0. The one namer that does honour `--basename`
(`paired_bam_output_name` → `foo_val.bam`) is reachable only on the two-file path, which is
guarded — and a BAM inside a two-file pair is rejected earlier by
`reject_bam_format_mismatch_in_pair` anyway.

**Test and CI non-vacuity.** All nine acceptance cases (7 integration + 2 unit) drive a guarded
dispatch path and assert success, so an unconditionally-rejecting helper fails every one — the
property the plan wanted does hold. The three new CI steps are each capable of failing: step 1
would see exit 0 plus three files without the SE guard; step 2 greps hint text only the hardtrim
sites can produce; step 3 greps the alias wording *and* md5-pins the named input. The
`validation` job is `ubuntu-latest` only (`ci.yml:315`) and already uses `md5sum` in ~14 places,
so the new `md5sum` introduces no portability problem. `collision_key`'s `Err` arm is unreachable
from the CLI: clap rejects an empty positional with exit 2 before `validate()` runs
(*reproduced*).

---

## 3. Issues

### HIGH

#### H-1 — Defect 2 is closed for primary outputs only; `--demux` still destroys a named input at exit 0

The candidate lists hold primaries plus the `--retain_unpaired` / `--passthrough` extras. Every
other file a run writes — trimming report, JSON report, clumping report, **demux per-barcode
outputs**, FastQC artifacts — is absent, so none is compared against the input set. A2 justifies
that for *output-vs-output* collisions, and that justification is sound (the new
`distinct_primary_outputs_imply_distinct_secondary_outputs` test now demonstrates it properly).
A2 says nothing about output-vs-**input**, which is the defect this change adds. So the second
half of the fix is narrower than §1 ("every prospective output path … checked … against every
input path") and narrower than the CHANGELOG, which claims the alias case "is now refused on
every path with a message of its own".

Reproduction (`t6.sh` scenario Q2 — plain `--demux`, no exotic flags):

```bash
cp test_files/demux_test.fastq.gz test_files/demux_test_samplesheet.txt .
trim_galore --demux demux_test_samplesheet.txt demux_test.fastq.gz     # exit 0
#   → demux_test_trimmed_sample1.fq.gz (1234 reads) among the per-barcode outputs
# a later re-run picks its own earlier output up in the input list:
trim_galore --demux demux_test_samplesheet.txt \
            demux_test.fastq.gz demux_test_trimmed_sample1.fq.gz
```

Observed **with the fix in place: exit 0, no collision message.** Input 2
(`demux_test_trimmed_sample1.fq.gz`, seeded with 10 marker reads for the test) is overwritten by
input 1's demux stage — 0 of the 10 marker reads survive — and is then read in its overwritten
state, so `demux_test_trimmed_sample1_trimmed.fq.gz` holds 1202 reads of the wrong sample. That
is precisely the failure mode the CHANGELOG describes as fixed ("destroyed the second input's
reads and wrote a … file containing the first sample's — a wrong-sample file, which is harder to
notice than a missing one"), with the same reachability story ("re-running over a directory
holding an earlier run's output").

*On the stale binary:* this was run against the 23:34 build, i.e. before D6. It is unaffected —
D6 is a pure extraction (`demux_base_name` recomputes `file_name()` from the same `trimmed_file`
and applies the identical `.gz`-then-`.fq` strips; `demux.rs:145-153`), `main.rs` is byte-identical
to the revision I reproduced on, and no demux path entered any candidate list. The behaviour under
the current source is the same.

The same class is reachable through the report namer, which shows the blind spot is structural
rather than demux-specific (`t7.sh` scenario W): `trim_galore a.fastq a.fastq_trimming_report.txt`,
where the second file is a valid FASTQ → input 1's report write clobbers input 2 (marker reads
gone), *then* the run fails on the corrupted input. Contrived, unlike the demux case, but it
generalises the finding to "every non-primary output".

**Recommendation** — one of:

1. Feed the derived secondaries into the **input-alias comparison only** (they need not join the
   duplicate-output comparison, which A2 covers). `report_name`, `json_report_name` and
   `clumping_report_name` are one line per site. Demux is the harder one: its per-barcode names
   need the barcode file read, which the pre-flight deliberately does not do — so demux
   realistically needs option 2.
2. State the residual the way A8 is stated — name it as an assumption in the plan, pin it with a
   test that documents the limitation, and **soften the CHANGELOG**: "refused on every path" →
   "refused wherever the trimmed/validated output itself is an input", plus one sentence that
   `--demux` per-barcode outputs are not yet covered.

Either way the CHANGELOG sentence must change: as written it claims something the code does not do.

---

### MEDIUM

#### M-1 — The hardtrim hint states something false whenever `--output_dir` is set, in the same sentence that prints the `-o` path

`HARDTRIM_HINT` (`main.rs:56-60`) opens with "Hard-trimmed output is written to the current
working directory". That is true only without `-o`; `hardtrim_output_name` ends
`Some(dir) => dir.join(filename)` (`specialty.rs:450-453`). *Reproduced* (`t9.sh`):

```
trim_galore --hardtrim5 20 -o <w>/out dirA/same.fastq dirB/same.fastq
Error: … <w>/out/same.20bp_5prime.fq and <w>/out/same.20bp_5prime.fq would be written to the
same file. Check that inputs produce distinct output paths (e.g., different source directories
or `--output_dir`). Hard-trimmed output is written to the current working directory, so inputs
sharing a basename collide whatever `--output_dir` is set to — …
```

The message names `out/…` and then denies that `out/` is where output goes. The operative clause
("collide whatever `--output_dir` is set to") is correct and worth keeping; the CWD claim is not.
Both `tests/integration_output_collision.rs:336`/`:361` and the new `ci.yml` hardtrim step use
exactly this `-o` shape, so the contradiction is currently asserted as correct behaviour rather
than caught.

**Recommendation:** reword the first clause to be true in both cases — "Hard-trimmed output is
named from the input's basename only, ignoring the input's directory" — or gate the CWD sentence
on `output_dir.is_none()`.

#### M-2 — The duplicate-output message names the same path twice and never names the inputs

`preflight_output_collisions` reports `{existing} and {p}`, i.e. the two *planned* paths. For the
#383 shape those are byte-identical strings, so the diagnostic is tautological and names neither
input. *Reproduced* on three paths:

```
# --hardtrim5 20 dirA/same.fastq dirB/same.fastq
… same.20bp_5prime.fq and same.20bp_5prime.fq would be written to the same file. …
# -o out sample.fastq sample.fq.gz
… out/sample_trimmed.fq and out/sample_trimmed.fq would be written to the same file. …
# --clump_only --cores 2 sample.fastq sample.fq
… sample_clumped.fq and sample_clumped.fq would be written to the same file. …
```

The identical-string case is the *common* case for the paths this change adds — it is #383 as
filed — whereas the inherited wording was written for #216, where the two paths differ by letter
case. On `--hardtrim5 30 */*.fastq.gz` over a per-sample layout (the pattern §3.4 exists for) the
user gets one repeated bare filename and has to work backwards to find which two of N inputs
produced it. The offending basename *is* recoverable, which is why this is Medium, not High.

**Recommendation:** carry the source input alongside each candidate (`&[(PathBuf, &Path)]`, or a
parallel index vector) and phrase the message around the inputs: "Inputs `dirA/same.fastq` and
`dirB/same.fastq` would both be written to `same.20bp_5prime.fq`." That also makes M-1's hint read
as a consequence rather than a rebuttal of the sentence before it.

#### M-3 — Case-folding now refuses genuinely distinct inputs on case-sensitive filesystems, undocumented

`collision_key` case-folds via `norm_path`, so the SE and hardtrim paths — which had **no** check
before — now reject two inputs whose stems differ only in case. *Reproduced*:
`trim_galore Sample.fastq sAMPLE.fq` → exit 1, "Output path collision (case-insensitive, for
APFS/NTFS safety)". On ext4 those are two files with two distinct outputs, and the run used to
succeed.

A1 records the trade but frames the false-positive surface as "opt-in case-sensitive APFS
volumes". ext4 is case-sensitive **by default** and Linux (containers, nf-core) is the primary
deployment target, so the surface is larger than A1 states — and it grew in this commit from
`--paired` to every path. The CHANGELOG does not mention it, while it does log the much rarer
duplicate-positional refusal under `#### Changes`.

**Recommendation:** no code change — the trade is right. Add one `#### Changes` bullet: inputs
whose output paths differ only in letter case are now refused on all paths, deliberately, for
APFS/NTFS safety, with the remedy. And widen A1's wording beyond APFS.

#### M-4 — No test or CI step covers hardtrim without `--output_dir`, the widest collision class

All four hardtrim rejection tests and the hardtrim CI step pass `-o`, under which the two
candidates collide simply because both are `<-o>/<same stem>` — the same mechanism as SE trim.
The behaviour that motivates the mode-specific hint, §2.2 repro (c) and §2.3's CWD asymmetry is
`--hardtrim5 20 dirA/same.fastq dirB/same.fastq` with **no** `-o`, and nothing exercises it.
§10 V2.5/V2.6 asked for `current_dir(tempdir)` cases; the implementation added `-o` on top, and
§13 does not record the deviation.

It does work — *reproduced* (`t1.sh`): exit 1 with the hint and nothing written into the CWD, and
the matching acceptance case (distinct basenames, no `-o`) exits 0 with correct per-file
attribution. So this is unpinned coverage rather than a defect: a future change that dropped the
`output_dir` argument from `planned_hardtrim_outputs` would still pass every hardtrim test today,
because the bare-filename form collides on the string alone.

**Recommendation:** drop `-o` from one of the two hardtrim FASTQ rejection tests and assert the
CWD afterwards holds only `dirA`/`dirB`.

---

### LOW

1. **`demultiplex`'s doc-comment was re-parented onto `demux_base_name` by D6.** The new function
   was inserted *between* `demultiplex`'s doc block and `demultiplex` itself, so `demux.rs:106-133`
   now reads as one comment — the five-step "For each read: 1. Extract … 5. Write to the matching
   sample's output file" plus "Also writes a summary file with per-barcode counts", followed by the
   new stem/A2 paragraphs, all attached to `demux_base_name` (`demux.rs:122`). `demultiplex`
   (`:134`) is now undocumented and `cargo doc` renders the demultiplexing algorithm as the
   documentation for a filename helper. Move the new function above the doc block, or below
   `demultiplex`. (D6's behaviour claim is otherwise sound — see the stale-binary note in H-1.)
2. **`primary_output_key_is_coarser_than_secondary_keys` still carries a claim it does not
   establish.** Its doc-comment says the shown property means "checking primaries covers the
   secondaries"; it actually exhibits a member of `~primary \ ~secondary`, which shows strictness,
   not containment. The new `distinct_primary_outputs_imply_distinct_secondary_outputs` is the
   test that carries A2, so the older one is now really a #382 stem-equivalence pin — relabel it
   and drop the coverage claim, otherwise a future reader will trust the wrong test. (This was
   Medium before D6 landed; the new test overtakes it.)
3. **The alias message prints the input's spelling as the path that would be written.**
   `io.rs:71-78` reports `input.display()`, not the planned path sharing its key. For a case-only
   or `./` alias those differ textually, so the user is told the run "would write output to
   `s_trimmed.fq`" when it would open `./s_trimmed.fq` or `S_trimmed.fq`. Naming both costs one
   `{}`.
4. **`HashMap<String, &PathBuf>` for `input_keys` keeps only the last of two aliasing inputs.**
   Harmless — the key *set* is unchanged, so whether the check fires is unaffected; only which
   input the message names. A one-line comment on the last-wins behaviour would save a future
   reader the analysis.
5. **The helper's call-site enumeration is incomplete.** `io.rs:53-56` says "every `main.rs`
   dispatch path that writes more than one file", but `run_paired_ubam_single_file` writes two
   FASTQ outputs and does not call it, and `run_ubam_output_paired_single_file` does not either.
   §4.1's stated purpose is that "a fifth omission [is] visible from here"; as written, an auditor
   finds two contradictions and no explanation. Name them as deliberately excluded (one input
   cannot duplicate, and `file_stem` + `_val` cannot alias). Also "all four output formats" should
   read "both output formats" — there are two.
6. **`--clock`/`--implicon` being excluded from the duplicate-input scan costs precision for
   cross-pair repeats.** `--clock A_R1 A_R2 A_R1 B_R2` is neither an R1==R2 nor a duplicate pair,
   so it falls through to the pre-flight and gets the APFS/NTFS wording (*reproduced*: exit 1,
   "Output path collision"). Correctly rejected, just not with the precise message §4.4 argues
   for. Fixable later by ordering the scan after `validate_paired_input` instead of excluding the
   modes.
7. **The duplicate-input message's argument numbers are positions in the input list, not argv.**
   `trim_galore -q 20 a.fq a.fq` reports "arguments 1 and 2" for what the user sees as words 4
   and 5. Say "input files 1 and 2", or drop the numbers.
8. **`tempdir()`'s `canonicalize(&d).unwrap_or(d)`** (`tests/integration_output_collision.rs:35`)
   silently returns to the vacuous-pass state D4 fixed if `canonicalize` ever fails. `.expect()`
   fails loudly instead — the point of D4 is that this exact silence cost an iteration.
9. **New CI directories are missing from the failure-artifact upload list.** `ci.yml:938-946`
   lists `/tmp/rust_case/` (the analogous #216 guard) but not `/tmp/rust_se_coll/`,
   `/tmp/rust_ht_coll/`, `/tmp/rust_alias/`. `/tmp/*.log` does pick up the three logs, so this is
   consistency only. Related: `docs/src/content/docs/modes/hardtrim.md:38` ("Both modes accept
   multiple input files in a single invocation") is now incomplete, since that shape is a hard
   failure when basenames repeat — §11 Open 4, deliberately deferred; line recorded so the
   follow-up is one edit.

---

## 4. Fixes applied

**None**, by instruction — see the header. H-1 and M-1 are the two I would act on before merge;
M-2 is a small signature change; M-3 is a CHANGELOG sentence; M-4 is a test change; the Lows can
ride along or wait.

---

## 5. Priority table

| # | Finding | Priority | Kind |
|---|---|---|---|
| H-1 | Non-primary outputs (esp. `--demux`) can still overwrite an input at exit 0; CHANGELOG overclaims | **High** | correctness + docs |
| M-1 | `HARDTRIM_HINT` asserts CWD output in the same message that prints the `-o` path | **Medium** | diagnostics |
| M-2 | Duplicate message names one path twice, never the inputs | **Medium** | diagnostics |
| M-3 | Case-fold now refuses distinct inputs on case-sensitive filesystems, unlogged | **Medium** | docs |
| M-4 | Hardtrim without `-o` (the widest class) untested | **Medium** | coverage |
| L-1…L-9 | See §3 | Low | polish |

---

## Appendix A — content found at this path before I wrote (preserved)

A 64-line in-progress report was at `CODE_REVIEW_B.md` when I came to write; full text kept at
`…/scratchpad/creview-b/CODE_REVIEW_B.preexisting.md`. Its two findings, both of which I then
verified independently, are carried above:

- its **B-1** (`HARDTRIM_HINT` false under `-o`) → my **M-1**, reproduced independently in `t9.sh`.
- its **B-2** (`demultiplex`'s doc-comment re-parented onto `demux_base_name` by D6) → my
  **L-1**, confirmed from `git diff src/demux.rs`.

Nothing else in that file was lost: it had empty Errors / Efficiency sections and an unfilled
summary.

---

## Appendix B — second independent B pass (author of the pre-existing file)

I am the reviewer whose in-progress file was replaced above. Rather than overwrite this report I
have finished my pass and record only (a) where I independently reached the same conclusion, and
(b) the findings not carried above. **No fixes applied.** My scratch scripts are in
`…/scratchpad/blind-b/` (`p1.sh`–`p5.sh`), run against a `target/release/trim_galore` copied at
mtime 12:57 (newer than every source file, and confirmed current by exercising the new SE guard).

### Independently reproduced, same conclusion

- **H-1.** Reached independently and by a different route: my `p3.sh` scenario H
  (`a.fastq` + `a.fastq_trimming_report.txt` holding valid FASTQ) destroyed the named input via
  the *report* writer, and `p5.sh` produced the exit-0 `--demux` version —
  `--demux bc.txt -q 32 src.fastq`, then `--demux bc.txt src.fastq src_trimmed_sample1.fq` →
  exit 0, `src_trimmed_sample1.fq` md5 changes from `503ad5d9…` (1224 reads) to `6e89e90c…`
  (1234 reads). Two reviewers converging on this from different fixtures is worth weighting.
- **M-1** (was my B-1) and **L-1** (was my B-2): confirmed as stated.
- **A9 across all eleven sites, the four converted sites, `guarded_inputs`, D1's predicate, the
  two unguarded N=1 paths, no fifth hole:** my call-by-call check and probes agree throughout.
  Adding to the above: `--paired --hardtrim5 20 a.fastq a.fastq` correctly yields "Read 1 and
  Read 2 appear to be the same file" (D1's exclusion of `paired` does not lose precision there),
  and `--paired --retain_unpaired` with an unpaired candidate fed back as a fourth input rejects
  with the alias wording and leaves the prior file intact.

### Findings not carried above

**B-3 (Medium) — `distinct_primary_outputs_imply_distinct_secondary_outputs` cannot fail, so it
does not carry A2.** H-1 credits it with demonstrating A2 "properly"; I disagree. Its three
fixture inputs (`d/alpha.fastq.gz`, `d/beta.fastq.gz`, `e/gamma.fq.gz`, `io.rs:1057-1061`) have
**pairwise-distinct basenames**, and all three secondary namers under test — `report_name`,
`json_report_name`, `clumping_report_name` — are injective in the input's filename. So
`namer(inputs[i]) != namer(inputs[j])` is true by construction for every pair the loop reaches,
independent of the pre-flight, the primaries, and the `basename`/`gzip`/`output_dir` matrix. The
two outer loop dimensions add nothing at all: none of the three namers takes `basename` or
`gzip`. The `demux_base_name` sub-assertion is injective in the primary's stem for the same
reason.

The shape that *could* fail is the one the fixture omits: two inputs with the **same basename in
different directories** (`d/same.fastq.gz`, `e/same.fastq.gz`). Without `-o` the primaries differ
and the reports must differ by directory — a real assertion; with `-o` both collapse and the
`a == b` guard makes the case `continue`. Adding that pair is a two-line change and turns a test
that cannot fail into one that can. As it stands this is v1's V6 failure mode in a new costume,
and it is the one the plan's §10 V5 was rewritten specifically to avoid.

**B-4 (Low) — `assert_rejected_cleanly`'s nothing-written check goes vacuous on a missing
directory.** `tests/integration_output_collision.rs:105-111` is
`std::fs::read_dir(dir).map(…).unwrap_or_default()`, so a `read_dir` failure yields an empty
`leftovers` and the assertion passes. Every current caller `create_dir_all`s the `-o` directory
first, so it is sound today; but this is the same silent-fallback shape as L-8's
`canonicalize(&d).unwrap_or(d)`, and D4/D3 are the record of what that costs. `.expect()` here
too.

**B-5 (Low, efficiency) — the new duplicate-input scan is O(n²) in the input count.**
`cli.rs:627-628` does `self.input[..i].iter().position(|other| other == path)` per input, i.e.
`n²/2` `PathBuf` comparisons; the pre-flight it complements is O(n). At a 5 000-file glob that is
~12 M path compares — still well under a second, and it matches the shape of the existing
duplicate-pair check, so this is a note rather than a request. If it is ever touched, a
`HashSet<&PathBuf>` makes it O(n) and the module already carries the imports.

**B-6 (Low) — the `.gz`/`.bgz` regression test writes plain text into gzip-named files.**
`se_trim_rejects_gz_and_bgz_sharing_a_stem` (`:196`) seeds `wide.fastq.gz` and `wide.fastq.bgz`
with uncompressed FASTQ. It does still pin #382's widening (pre-#382 stem-stripping left
`wide.fastq` on the `.bgz` arm, so no collision), but because `is_gzipped` is content-based the
run resolves to `gzip = false` and both candidates end `.fq`, not `.fq.gz` — so the case never
exercises the gzip-suffixed spelling it is named for. Using real gzip/bgzip bytes would cost one
`flate2` write in the fixture helper.

### My counts

0 Critical, 1 High (H-1, converged), 3 Medium (M-1 converged, plus B-3; I did not independently
assess M-2/M-3/M-4, which look right to me on reading), 3 Low unique to this appendix.
