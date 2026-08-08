# Code Review B — paired report pre-flight (#388)

**Reviewer:** B (independent; Reviewer A reviews the same diff in parallel, no shared state)
**Target:** `fix/388-paired-report-preflight` off `dev` @ `c4599f5`, **uncommitted**
**Diff:** `src/main.rs` (2 blocks, +15), `tests/integration_output_collision.rs` (+151, T1–T6), `CHANGELOG.md` (+17)
**Plan:** `plans/08072026_paired-report-preflight/PLAN.md` (v2 + §8)

**Nothing was fixed.** The tree is shared and uncommitted, so this report is recommend-only.
Every behavioural claim below was run against `target/release/trim_galore` (up to date with
the diff); harness in
`/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/d0ee0171-fc72-476b-bdba-c46c90a5bacd/scratchpad/crev388b/`.

---

## Verdict

**Correct and complete for the two paths it touches.** Candidates match the writers exactly,
the pre-flight sites are the right ones, the tests are discriminating, and no over-rejection
shape survived a sweep. No Critical or High finding on the diff.

One **Medium** is adjacent rather than in-diff: the identical defect — positionally
discriminated primaries plus non-discriminated per-input reports — is live on
`--clump_only --paired`, reproduced below. PLAN §5 A4 records the uBAM-SE residual but not
this one, so it currently reads as covered when it is not.

---

## Focus item 1 — do the candidates match the writers? **Yes, exactly.**

| | candidate site | writer site |
|---|---|---|
| paired FASTQ | `main.rs:765-772` | `run_paired`, `main.rs:1566-1607`; paths at `:1578-1585` |
| paired uBAM out | `main.rs:1867-1873` | `run_ubam_output_paired_two_files`, `main.rs:2164-2202`; paths at `:2178-2185` |

- **Same namers, same arguments.** Both writers call `naming::report_name(input, output_dir)`
  / `json_report_name(input, output_dir)` — the same two-argument form. Notably **neither
  writer applies `cli.basename`** to report paths, and neither does the candidate block:
  no mismatch. (`io.rs:485-514`: both namers key on the *full input filename*, not the stem,
  and ignore `basename` entirely.)
- **Both sides, both formats.** `write_paired_reports` loops `[(stats_r1, &r1), (stats_r2, &r2)]`
  (`main.rs:2382`) writing `txt_path` then `json_path` each — four paths per pair, four
  candidates.
- **Same gate.** Writers gate at `:1566` / `:2164` on `!cli.no_report_file`; candidates use
  the same expression. Confirmed empirically under Focus item 4.
- **No writer left uncovered on these paths.** The paired FASTQ write set is `_val_{1,2}`,
  `_unpaired_{1,2}` under `--retain_unpaired`, the passthrough output (all pre-existing), the
  four reports (now added), and FastQC artifacts (excluded per A3). `--passthrough` folds its
  stats into the **R2** report (`main.rs:2416-2418`) rather than writing a third, so there is
  no fifth path; `--clumpify` on the trim path writes no clumping report
  (`clumping_report_name` is referenced only from `clump_only.rs`). The uBAM paired path
  writes one `_val.bam` plus the same four reports.

## Focus item 2 — single-file interleaved-uBAM paired paths: **correctly out of scope**

`run_paired_ubam_single_file` (`main.rs:1643`, BAM in → FASTQ out) and
`run_ubam_output_paired_single_file` (`:2222`, BAM in → BAM out) do both write reports
(`:1786-1810`, `:2292-2311`) and neither has a pre-flight at all — but both are reached only
under `cli.paired && cli.input.len() == 1` (`:716`, `:1839`), i.e. **one** input. One input
yields one report pair, so no output-vs-output collision is constructible, and
`<bam-filename>_trimming_report.txt` cannot equal the `.bam` input, so the output-vs-input
branch has nothing to catch either. Adding candidates there would be dead code. No finding.

## Focus item 3 — `guarded_inputs` / output-vs-input: **no wrongly-tripped alias branch**

The new candidates are also matched against `guarded_inputs` (`main.rs:64` = inputs +
`--passthrough` + `--demux`; the latter two are single-pair-only / SE-only). A report
candidate can only alias an input if one input is literally named
`<other-input-filename>_trimming_report.txt` — a shape already rejected today, by
`sanity_check_any` on `input[0]` before the pre-flight and on `chunk[1]` at `:795` after it.
So the only change is *which* message `--paired a.fq a.fq_trimming_report.txt` gets (now the
accurate "also one of its inputs" + "drop it from the input list"). Diagnostic improvement,
not a regression.

Over-rejection sweep — A1 holds on every shape I could construct:

```
n1  --paired A/reads.fq B/reads.fq                  → exit 0; T1's inputs minus -o. Both
                                                      report pairs written beside their inputs.
n2  --paired -o out a/w.fq a/x.fq b/y.fq b/z.fq     → exit 0; 4 _val_ + 8 reports in one dir.
n3  --paired -o out a/w.fq a/dup.fq b/y.fq b/dup.fq → exit 1; correct: two inputs named
                                                      dup.fq would write one report path.
t6a --paired --basename foo r1.fq r2.fq             → exit 0; foo_val_{1,2} + r1/r2 reports.
t6b --paired sample.fq sample.fastq                 → exit 0; sample_val_{1,2} + 4 reports.
```

`n1` is the important one: T1 with `-o out` removed still runs, so the new candidates reject
on the shared *output directory*, not on the shared filename.

## Focus item 4 — the six tests: **all six can fail; both acceptance cases are at the boundary**

**Negative-control spot-check without touching the tree.** Commenting the block out is
impossible on a shared tree, so I re-ran each rejection shape with `--no_report_file` added —
the only input that disables the new block. Each flips to exit 0, proving the rejection is
caused by the new candidates and nothing else:

```
t1  --paired -o out A/reads.fq B/reads.fq                           → exit 1, dup-msg
t2  --paired -o out a/r1.fq b/R1.fq                                 → exit 1, dup-msg
t2b --paired --no_report_file -o out a/r1.fq b/R1.fq                → exit 0, r1_val_1 + R1_val_2
t3b --paired --no_report_file -o out a/r1.fq a/r2.fq b/R2.fq b/x.fq → exit 0, 4 _val_ files
t4b --paired --output-format ubam --no_report_file -o out A/… B/…    → exit 0, reads_val.bam
```

This is *stronger* than §6's comment-out for T1–T4: it also proves the gate is wired to the
same flag as the writers (A2). It confirms §8's claim.

- **T1** is genuinely case-free — primaries differ by the positional digit under any fold, so
  only the reports can collide. **T4** is its uBAM twin; `t4b` shows the FASTQ block cannot be
  what rejects it (single `_val.bam` primary).
- **T2 / T3** depend on case-folding but assert only the *rejection*, so they are
  filesystem-independent as documented. Neither rejects for the wrong reason: T2's primaries
  fold to `r1_val_1`/`r1_val_2`, T3's to four distinct keys, so the reports are the sole
  colliding pair in both.
- **T5**'s `read_dir(&out).count() == 2` is **robust here**: `out` is a fresh subdirectory of a
  per-tag temp dir, and with reports off and FastQC unrequested the run's entire write set
  into `-o` is the two `_val_` files — confirmed by `find -type f` on `t2b` yielding exactly
  two. The two `is_file()` assertions pin *which* two, so a wrong-namer regression cannot
  satisfy the count alone.
- **T6 case 1** (`--basename foo`, no `-o`) catches a candidate list keyed on `basename` like
  the primaries: both reports would then be `foo…` and collide. Real boundary — `t6a`
  confirms the primaries are `foo_val_{1,2}` while the reports keep `r1.fq`/`r2.fq`.
- **T6 case 2** (`sample.fq` + `sample.fastq`) catches a candidate list keyed on the *stem*:
  both stems are `sample`, so a stem-keyed list would reject a legal pair. Real boundary, and
  the only case in the suite separating `report_name` from `single_end_output_name` semantics.

Gaps, all minor: L-1, L-2, L-4.

## Focus item 5 — duplication: see M-2.

---

## Findings

### Medium

**M-1 — the same defect is live on `--clump_only --paired`; not fixed, not scoped out.**

`clumped_paired_output_names` (`io.rs:387-412`) appends the positional discriminator
`_clumped_1`/`_clumped_2` while `clumping_report_name` (`io.rs:469-482`) keys on the full input
filename with no discriminator — structurally identical to `_val_1`/`_val_2` vs `report_name`,
i.e. **exactly** the §2 root cause. `clump_only_paired` writes two per-mate reports gated on
`no_report_file` (`clump_only.rs:534-549`), but the mode dispatches through
`run_specialty_paired` (`main.rs:2462`), whose pre-flight pushes **only** the two names the
closure returns (`:2474-2480`). Reproduced on this branch:

```
$ trim_galore --clump_only --paired -o out A/reads.fq B/reads.fq
exit=0
out/reads_clumped_1.fq
out/reads_clumped_2.fq
out/reads.fq_clumping_report.txt        ← ONE report for two inputs
$ head -3 out/reads.fq_clumping_report.txt
… Input:  B/reads.fq                    ← the survivor is R2's; R1's was overwritten
```

Control: `--clump_only --paired -o out A/p1.fq B/p2.fq` writes both reports, so the mode is
otherwise well-behaved. This is the #388 failure mode verbatim, one mode over. Out of the
plan's literal scope, but §5 A4 records the uBAM-SE residual and is silent here.
**Recommendation:** either give `run_specialty_paired` an optional secondary-name closure and
pass the clumping-report names from the `--clump_only --paired` arm, or file a follow-up issue
and add it to A4 as a known residual. Do not leave it undocumented.

**M-2 — candidate-side report naming is now duplicated across two sites while the writer side
is centralised.**

The two inserted blocks are near-identical; the only difference is the vector pushed into
(`candidates` vs `planned`). This matters more than seven duplicated lines suggest, because
the *writer* side was deliberately centralised for exactly this reason —
`write_paired_reports`'s doc comment says it exists "so we don't drift against
`run_paired_ubam_single_file`'s identical block" (`main.rs:1564-1565`). The candidate side now
carries the drift risk the writer side was refactored to remove: a change to report naming or
to the gate must be mirrored twice, and a partial edit reintroduces #388 on whichever path is
missed. There is also an in-repo precedent for the fix: `planned_secondary_outputs`
(`main.rs:83-104`) already computes precisely these four paths for the SE FASTQ path with the
identical gate. **Recommendation (structure, not a bug):** a small
`fn paired_report_candidates(cli, r1, r2, output_dir) -> Vec<PathBuf>` called from both sites,
or generalise `planned_secondary_outputs` to take an explicit input slice. Low risk, and it
makes the two paths provably identical rather than identical by inspection.

### Low

**L-1 — T2/T3/T4 do not assert that a refused run wrote nothing.** T1 checks
`read_dir(&out).next().is_none()`; T2–T4 assert only exit status and message. The file already
has `assert_rejected_cleanly` (`tests/integration_output_collision.rs:95`), which asserts
non-zero exit **plus** `PREFIX` **plus** an empty residue — and the module docstring at `:5-9`
calls "two reports beside one data file" #383's most misleading artifact, which is the
invariant those three leave unchecked. It would also restore the `PREFIX` assertion the new
tests drop (they match only the `DUP_MSG` fragment, which still passes if the message loses its
"case-insensitive, for APFS/NTFS safety" framing). **Recommendation:** route T1–T4 through
`assert_rejected_cleanly(&out, ok, &stderr, DUP_MSG)` — a drop-in for all four.

**L-2 — T3's comment omits the filesystem caveat T2's states.** T3 relies on `a/r2.fq` vs
`b/R2.fq` folding together, so on a case-sensitive filesystem it is the same #216-style false
positive T2's comment carefully explains, but T3's one-liner ("the cross-pair route") does not
say so. **Recommendation:** one clause, or make T3 case-free — verified working inputs:
`--paired -o out a/w.fq a/dup.fq b/y.fq b/dup.fq` rejects on any filesystem (`n3` above).

**L-3 — one CHANGELOG sentence overstates its scope.** "two genuinely distinct inputs whose
filenames differ only in case are refused even on a case-sensitive filesystem" is true of
paired *reports* but reads generally; SE has been like this since #385. Consider "…are refused
as a paired pair even on…". Wording only — both §3 behaviour changes are disclosed, which is
the part that matters.

**L-4 — T1, the fix's raison d'être, has no acceptance twin.** The module docstring states the
suite's own convention at `:11-14`: "Every rejection case here is paired with an acceptance
case on the same dispatch path and output format." T1's nearest neighbour — the *same* two
inputs with `-o out` removed — is asserted nowhere, and T6's two cases are different shapes.
Without it, a candidate list that keyed reports on the filename while ignoring `output_dir`
would reject a legal run and the suite would stay green. Verified the case passes today (`n1`:
exit 0, `A/reads_val_{1,2}.fq` plus a report pair in each of `A/` and `B/`).
**Recommendation:** add it as a third T6 case — four lines.

### Efficiency / errors — no finding

Four `PathBuf` allocations per pair in a loop that already allocates two to five, once before
any I/O; `preflight_output_collisions` is O(n) over a pre-sized map (`io.rs:103`). No new panic
path: `chunk[0]`/`chunk[1]` sit inside `cli.input.chunks(2)` loops whose length `Cli::validate`
guarantees even, and the `len() == 1` paired case is diverted at `:716`/`:1839` before either
block runs. No new error path: both namers are infallible (`unwrap_or_default` on `file_name`).

---

## Recommendation summary

| ID | Priority | Action |
|---|---|---|
| M-1 | Medium | Fix `--clump_only --paired` too, **or** file a follow-up and record it in PLAN §5 A4 |
| M-2 | Medium | Extract the duplicated candidate block into one helper (drift risk the writers were refactored to remove) |
| L-1 | Low | Route T1–T4 through the existing `assert_rejected_cleanly` |
| L-2 | Low | Note T3's filesystem sensitivity, or make it case-free |
| L-3 | Low | Narrow one CHANGELOG sentence to the paired-report scope |
| L-4 | Low | Add T1's acceptance twin (same inputs, no `-o`) as a third T6 case |

None of these block the diff. M-1 is the one that should not stay silent, because it is the
same bug class the plan says it is closing.
