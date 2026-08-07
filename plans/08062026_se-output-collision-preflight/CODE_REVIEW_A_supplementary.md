# Code Review — output-collision pre-flight (#383) — supplementary

> **On the filename.** This review was first written to `CODE_REVIEW_A.md`, which turned out to
> still be owned by another reviewer agent that later wrote its own review there. That file is a
> different agent's **independent** review, not an earlier version of this one — the two are not
> versions of each other and should be read as two separate opinions. This copy is unchanged in
> substance from what I originally wrote; only the filename differs.

**Target:** uncommitted working tree on `fix/383-output-collision-preflight` off `dev` @ `72624c7`
**Files:** `src/io.rs`, `src/main.rs`, `src/cli.rs`, `src/specialty.rs`, `src/demux.rs`,
`.github/workflows/ci.yml`, `CHANGELOG.md`, `tests/integration_output_collision.rs`
**Date:** 2026-08-07

## Caveats on the weight of this review

1. **Not blind.** I audited this same change against its plan earlier in this session, so I
   have prior exposure to the plan's reasoning, its assumption list (A1–A10) and its
   deviation log (D1–D6). A reviewer who had read the design is measurably less likely to
   question the design. Findings below that agree with the plan should be discounted
   accordingly; findings that contradict it are the ones worth the most.
2. **Not dual.** The workflow asks for two independent reviewers in separate context
   windows. This is one, and it is the same agent that produced the coverage audit. A second
   independent reviewer is still owed and should be run before this lands. This is also why the
   filename is "supplementary" rather than "A".
3. **No fixes applied.** The caller requires the diff to stay stable, which overrides the
   skill's "fix directly" clause. Everything below is a recommendation.

## Summary

The change is sound where it matters most. The part I expected to find broken — whether each
of the eleven guards hashes the paths its writer actually writes — is correct at all eleven,
verified expression by expression (C1). `guarded_inputs()` is correct on the one path where
`--passthrough` can be set and provably vacuous on the other ten (C2). D6's `demux_base_name`
extraction is byte-for-byte behaviour-preserving, so the Perl `--demux` byte-identity path is
safe (C3). Both defect fixes work on invocations the test suite does not construct.

**1 High, 4 Medium, 4 Low.** No Critical.

The High finding and two of the Mediums are one root cause: the diff establishes
`collision_key` as this crate's definition of "the same file" but leaves three sibling identity
checks in `cli.rs` on weaker definitions. The sharpest consequence is that
`trim_galore --paired ./a_R1.fq a_R1.fq` exits 0 and writes two byte-identical files labelled
as a validated R1/R2 pair (H1) — pre-existing, not a regression, but this is the commit that
makes it an inconsistency rather than a uniform limitation. Making `collision_key` `pub` and
using it at those three sites closes H1 and M3 together.

The remaining Mediums are independent: `guarded_inputs()` omits the `--demux` barcode file, and
a run whose trim output aliases it destroys that file before any check (M2, reproduced); and
D6 detached `demultiplex`'s doc-comment onto the new helper (M1) — the same orphaning that
Step 2 of this very plan was fixing elsewhere.

Five probes were run from a scratch directory, never the repo root; three reproduced findings
and two confirmed correct behaviour. Nothing was fixed, per instruction.

## Findings

### Clean bills of health (the three things I tried hardest to break)

**C1 — All 11 call sites guard the paths their writers actually write.** This was the
highest-value check available, because a mismatch is invisible to every rejection test: the
guard would hash paths the run never touches, reject nothing, and the acceptance tests would
still pass. I paired each candidate expression against its writer's naming expression:

| Guard | candidate expression | writer's expression | verdict |
|---|---|---|---|
| `main.rs:360` hardtrim5 | `hardtrim_output_name(input, n, Five, output_dir, gzip)` / `hardtrim_bam_output_name(input, n, Five, output_dir)` | `specialty.rs:45` / `:139` — same args | match |
| `:391` hardtrim3 | as above with `Three` | `specialty.rs:88` / `:191` | match |
| `:513` clump_only SE FASTQ | `clumped_output_name(input, output_dir, basename, gzip)` | `clump_only.rs:278`, param `gzip_output` receives `gzip` (4th positional, `clump_only.rs:250-261`) | match |
| `:548` clump_only PE uBAM, 1 file | `clumped_paired_bam_output_name(&cli.input[0], None, output_dir, basename)` | `clump_only.rs:922` `(input, None, output_dir, basename)` | match |
| `:585` clump_only PE uBAM, N pairs | `(&chunk[0], Some(&chunk[1]), output_dir, basename)` | `clump_only.rs:969` `(r1_path, Some(r2_path), …)` | match |
| `:633` clump_only SE uBAM | `clumped_bam_output_name(input, output_dir, basename)` | `clump_only.rs:754` | match |
| `:734` paired FASTQ | `paired_end_output_names` + `unpaired_output_names` + `passthrough_output_name` | `clump_only.rs:416` / `main.rs:1305` | match |
| `:799` SE trim FASTQ | `single_end_output_name(input, output_dir, cli.basename.as_deref(), gzip)` | `main.rs:1103` — character-identical | match |
| `:1825` paired uBAM | `paired_bam_output_name(&chunk[0], &chunk[1], output_dir, cli.basename.as_deref())` | `main.rs:2028` | match |
| `:1869` SE trim uBAM | `single_end_bam_output_name(input, output_dir, cli.basename.as_deref())` | `main.rs:1893` — character-identical | match |
| `:2430` `run_specialty_paired` | caller's `output_names` closure | `specialty.rs:254-255` (clock), `:367-368` (implicon), `clump_only.rs:416` | match |

The one argument that could plausibly have desynchronised the hardtrim pair — `cli.rename` —
does not reach naming at all: it gates a per-read header rewrite inside the record loop
(`specialty.rs:60`, `:104`, `:164`, `:217`), never the output path. The clock/implicon
closures capture the same `output_dir`, `gzip` and `umi_len` bindings their worker closures
do, so they cannot drift.

**C2 — `guarded_inputs()` is harmless where `--passthrough` cannot be set, and correct where
it can.** `cli.rs:885-921` confines `--passthrough` to `--paired` with exactly one pair, and
rejects it alongside `--clock`, `--implicon`, `--hardtrim5/3`, `--demux`, `--clump_only`,
`--clumpify`, `--retain_unpaired` and uBAM output. So on ten of the eleven call sites
`cli.passthrough` is provably `None`, and on the eleventh the paired FASTQ candidate list
(`main.rs:724-731`) adds its output. The `--retain_unpaired` exclusion also means the two
conditional candidate blocks at `:709` and `:724` can never both fire.

**C3 — D6's `demux_base_name` extraction is behaviour-preserving.** `demux.rs:122-135` is the
original statement sequence verbatim: same source (`file_name().unwrap_or_default()
.to_string_lossy().to_string()`), same two `ends_with` strips in the same order, same
non-idempotent single-pass semantics (`x.fq.gz.gz` → `x.fq`, unchanged). `trimmed_name` is
still live at `demux.rs:232` for the progress line, so nothing was orphaned and clippy stays
quiet. The `--demux` byte-identity path in the Perl validation matrix is safe. The one problem
with D6 is structural, not behavioural — see M1.

### Critical

None found. Nothing in the diff loses reads on a path the diff itself introduced, and the
"nothing written" contract holds on all eleven guarded paths (probe 1 and probe 5 below).

### High

**H1 — `--paired` accepts one file as its own mate and silently emits a fake pair
(exit 0, no warning).** Not introduced here, but this diff is the change that decides what
"the same file" means in this crate, and it leaves the R1≠R2 check on the weakest notion.
`cli.rs:534-541` compares `chunk[0] == chunk[1]` — raw `PathBuf` equality — so any two
spellings of one path slip through, and the pre-flight cannot catch it afterwards because
`_val_1` and `_val_2` are genuinely distinct output paths and neither aliases an input.
Measured:

```
$ trim_galore --paired ./a_R1.fq a_R1.fq      # 20-read single-end file
rc=0, no warning
a_R1_val_1.fq  20 reads  md5=984ac651…
a_R1_val_2.fq  20 reads  md5=984ac651…      ← byte-identical to val_1
```

Two byte-identical files presented as a validated R1/R2 pair. An aligner maps that as proper
pairs where every fragment is a self-pair, and the trimming report looks normal. This is the
same failure class as #383 — exit 0 with wrong data — and the plan's own words for it are
"harder to notice than a missing one". The duplicate-*pair* check at `cli.rs:545-557`
(`r1 == pr1 && r2 == pr2`) has the identical weakness.
**Recommendation:** make `collision_key` `pub` and key all three `cli.rs` identity checks on
it. One-line change per site. See M3 — H1 is the sharp end of the same root cause.

### Medium

**M1 — D6 detached `demultiplex`'s doc-comment and hung it on `demux_base_name`.**
`demux.rs:106-115` is `demultiplex`'s doc-comment (the 5-step algorithm and "Also writes a
summary file with per-barcode counts"); `demux_base_name` was inserted at `:116` *inside* it,
so rustdoc now attributes all of it to the new two-line helper, whose own doc-comment reads as
a continuation of it, and `demultiplex` at `:137` has no doc-comment at all. This is precisely
the defect §5 Step 2 set out to repair in `main.rs`, where deleting `preflight_collision_bam`
reunited `resolve_clump_layout` with its orphaned comment. Behaviourally inert; wrong in
`cargo doc` and misleading to read. **Recommendation:** move the new function above
`:106`, or below `demultiplex`.

**M2 — `guarded_inputs()` omits `cli.demux`, and the omission destroys the barcode file
before any check runs.** `main.rs:62-68` builds the guarded set from `cli.input` plus
`cli.passthrough`. `cli.rs` has exactly four path-typed fields — `input`, `output_dir`,
`passthrough`, `demux` — so `--demux`'s barcode file is the one input the run reads that the
pre-flight never sees. The plan reasoned it away by argument (§3.2: output names always carry
a `_trimmed`/`_val_` segment, so they cannot equal a barcode text file "except by deliberate
contrivance"). Measured, the contrivance is destructive:

```
$ trim_galore --dont_gzip --demux sample_trimmed.fq sample.fastq
barcode file md5 before=8a96bd56… after=984ac651…   ← overwritten with trimmed reads
rc=1 (demux then fails to parse the file it just destroyed)
```

`read_barcode_file` is not called until `main.rs:1266`, *after* trimming has written
`output_path`, so the pre-flight is the only thing that could have prevented this — and the
run also breaches the "nothing written" contract, leaving two report files behind. I am
sizing this Medium rather than High because reachability needs a barcode file named exactly
like a trim output, which no sane workflow produces. But the argument was unnecessary:
`guarded_inputs()` exists precisely to be the answer, and adding `cli.demux` to it is one
line that retires the reasoning entirely. `guarded_inputs`'s doc-comment ("Input paths a run
must not overwrite: the positionals plus `--passthrough`") also overstates its completeness.

**M3 — three different notions of "same file" now coexist, and the weakest ones are the
user-facing ones.** After this diff the crate compares paths three ways: raw `PathBuf`
equality (`cli.rs:534`, `:547`, and the new duplicate-input check at `:627`), `norm_path`
(case-folded only — the `--passthrough`-aliases-R1/R2 guard at `cli.rs:931-933`), and
`collision_key` (case-folded *and* absolutised — the pre-flight). The diff upgraded the third
and left the first two, so defect 1's own analysis now applies verbatim to two checks sitting
one layer above the one that was fixed. The visible cost:

```
$ trim_galore ./s.fastq s.fastq
Output path collision (case-insensitive, for APFS/NTFS safety): … would be written to the same file.
```

That is the generic APFS/NTFS message, not the precise "was given more than once" — and §4.4
introduced the `Cli::validate` check specifically so this case would get the precise one,
quoting `cli.rs:504-516`'s doc-comment to that effect. It fires only for byte-identical
spellings, so the plan's stated goal is half-delivered. No data is lost (the pre-flight still
rejects), so this is a diagnostic-quality finding on its own; it is Medium only because it
shares a root cause with H1, where the same weakness does cost correctness.
**Recommendation:** `pub fn collision_key`, then use it at all three `cli.rs` sites. That
closes M3 and H1 together and leaves one definition of file identity in the crate.

**M4 — a seven-line comment above the paired call site documents the key that was just
replaced.** `main.rs:692-698` still reads "Hash key is the full path, case-folded (ASCII
lowercase) … (issue #216)", followed by the "Pragmatic trade-off" paragraph — sitting directly
above `preflight_output_collisions`, which no longer hashes the full path as-is. A reader
debugging a future collision report will take the absolutisation as absent. The same trade-off
paragraph is also now duplicated verbatim in `norm_path`'s doc-comment (`io.rs:40-42`), which
is where it belongs. **Recommendation:** delete `:692-698`; the helper documents itself.

### Low

**L1 — `HashMap<String, &PathBuf>` over inputs does drop an entry, but it cannot cause a
missed rejection.** `inputs.iter().map(…).collect()` keeps the last value for a repeated key,
so when two inputs alias one file the map holds one of them. The map is used only for
membership (unaffected — the key is still present) and to name the offending path in the
message, so the sole consequence is which of the two aliases is printed. Reachable via
`./x.fq x.fq`, and probe 1 confirms the run is still rejected. Worth one word in the message
("one of its inputs") rather than a code change — the current wording already avoids claiming
which argument index it was.

**L2 — `collision_key`'s `Err` fallback is a silent downgrade to the pre-fix key.**
`std::path::absolute(p).unwrap_or_else(|_| p.to_path_buf())` fails on an empty path or a failed
`getcwd`. Degrading to the raw string is the right call — strictly better than skipping the
check — and both triggers are effectively unreachable: an empty positional dies in
`detect_input_format` before dispatch, and a `getcwd` failure is global, so every key degrades
together rather than inconsistently. The inline comment at `io.rs:49` explains what the key
*does*, not that the fallback silently reverts to the behaviour #383 was filed about. Worth
half a sentence there.

**L3 — one violation reported per run.** The loop bails on the first offending planned path,
so an invocation with both a duplicate output and an input alias reports only whichever comes
first in `planned` order (measured: the duplicate). This is §3.3.3 as specified, and
collecting all violations would complicate the two message shapes; noting it only because a
user fixing a glob may need two round trips.

**L4 — the `"R1"`/`"R2"` discriminators stayed `&str` while hardtrim's became an enum.**
`HardtrimEnd` (§4.2) exists because a copy-pasted `"5prime"` in the `--hardtrim3` block would
make the guard hash paths the run never writes — a hole invisible to every rejection test.
`clock_output_name` and `implicon_output_name` take the same kind of discriminator as a bare
`&str`, duplicated between the `output_names` closure (`main.rs:426-427`, `:440-441`) and the
writer (`specialty.rs:254-255`, `:367-368`). The hazard is real but milder: swapping R1/R2 in
the closure yields the same hash *set*, and doubling one value causes a false rejection that
the acceptance tests catch. Left as-is is defensible; it is now the only untyped instance.

### Efficiency

Nothing actionable. `preflight_output_collisions` is O(n+m) with one `String` per path, and
`HashMap::with_capacity(planned.len())` is already there. `std::path::absolute` costs one
`getcwd` per relative path, so n+m syscalls at command-line cardinality — immeasurable, and
paid before adapter auto-detection's ≤1 M-record scan rather than after it, which is a net win
on the rejecting path. Building `input_keys` unconditionally does hash the input list even when
`planned` is empty or length 1; guarding that would trade a branch for a few hundred
nanoseconds and is not worth the line.

## Fixes applied

None, by instruction.

## Recommendations by priority

| # | Priority | Recommendation | Cost | Ship-blocking? |
|---|---|---|---|---|
| H1 | High | `pub fn collision_key`, then key `cli.rs:534` (R1≠R2), `cli.rs:547` (duplicate pair) and `cli.rs:627` (new duplicate input) on it instead of `==` | 4 lines + 1 test | **No** — pre-existing, not a regression. But it is the cheapest correctness win adjacent to this diff, and leaving it makes the commit internally inconsistent. Own issue if not taken here. |
| M2 | Medium | Add `cli.demux` to `guarded_inputs()`; correct its doc-comment | 3 lines | No |
| M3 | Medium | Covered by H1's change; also move `cli.rs:931-933`'s `--passthrough` alias check from `norm_path` to `collision_key` | 1 line | No |
| M1 | Medium | Move `demux_base_name` out of `demultiplex`'s doc-comment | move 20 lines | No |
| M4 | Medium | Delete the stale key description at `main.rs:692-698` | −7 lines | No |
| L2 | Low | One clause at `io.rs:49` noting the `Err` fallback reverts to the pre-#383 key | 1 line | No |
| L1, L3, L4 | Low | Noted for the record; no change recommended | — | No |

**Verdict:** nothing here blocks the commit. H1 is the one I would not leave unfiled — it is a
silent-wrong-output bug in `--paired`, it is one line from being fixed, and the fix belongs with
this change rather than after it, because this is the commit that decides what "the same file"
means. M1 and M4 are tidy-ups I would fold in before pushing, since both are artifacts of this
diff. M2 is a judgement call the author has already made once, consciously; my only argument is
that the one-line fix is cheaper than the paragraph of reasoning defending its absence.

## Probes run

All from `/private/tmp/…/scratchpad/creview-warm/`, using a 20-read slice of
`test_files/illumina_10K.fastq.gz`. Script retained at `creview-warm/probe.sh`.

| Probe | Invocation | Result |
|---|---|---|
| 1 | `trim_galore ./s.fastq s.fastq` | rc=1, duplicate-output wording (not the precise one) → M3; nothing written |
| 2 | `trim_galore --paired ./a_R1.fq a_R1.fq` | **rc=0**, two byte-identical `_val_` files, no warning → H1 |
| 3 | `trim_galore --dont_gzip --demux sample_trimmed.fq sample.fastq` | barcode file **overwritten**, then rc=1 → M2 |
| 4 | duplicate output *and* input alias in one command line | duplicate wording wins (first in `planned` order) → L3 |
| 5 | `trim_galore --dont_gzip --hardtrim5 20 x.fastq x.20bp_5prime.fq` | rc=1, alias wording, aliased file **intact** → correct |

One process note: an early probe attempt wrote a fixture into the repo root because agent-thread
`cd` does not persist between tool calls. It was deleted immediately and `git status` re-verified
against the 7-modified-file baseline before continuing; the diff under review was never touched.
