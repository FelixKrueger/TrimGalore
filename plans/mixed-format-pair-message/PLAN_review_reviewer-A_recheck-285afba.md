# Reviewer A — re-verification after the stale-tree correction

**Supersedes:** the verdict in `PLAN_review_reviewer-A.md`, which was produced against `a4ffd47`.
**Tree re-verified against:** `/Users/fkrueger/Github/TrimGalore`, branch `dev`, HEAD `285afba`, working tree clean for `src/` and `tests/`.
**Binary re-verified with:** freshly rebuilt, provenance line `285afba — macos/aarch64 — built 2026-08-13T20:16:22Z`.

## 1. What was actually wrong with my environment

The correction was right that my verification was stale, and wrong in every specific. Correcting the record because the difference matters for what to do next:

- **I was not in a worktree.** `git rev-parse --show-toplevel` → `/Users/fkrueger/Github/TrimGalore`; `git rev-parse --git-dir` → `.git` (a real directory, not a worktree gitdir link); `git worktree list` shows this path as the main tree on `[dev]`. The two worktrees sitting at `c6528e5` (`agent-a60f44ce66ae70d80`, `agent-a8a2dbce834e3d868`, both `locked`) belong to other agents.
- **I therefore did not run `git reset --hard 285afba`.** In this tree it would have been a no-op on the commit pointer (HEAD is already `285afba`) while discarding any tracked-file modification in the user's main working tree. Destructive for no gain.
- **The stale commit was `a4ffd47`, not `c6528e5`,** and the drift was ~20 commits, not 7 weeks. Proof by line count rather than by trusting either of us:

  | file | measured at review time | `a4ffd47` | `285afba` = disk now |
  |---|---|---|---|
  | `src/main.rs` | 2400 | **2400** | 2694 |
  | `src/cli.rs` | 1879 | **1879** | 2415 |
  | `src/clump_only.rs` | 1681 | **1681** | 1698 |
  | `src/format.rs` | 206 | **206** | 887 |

  Exact match on all four against `a4ffd47`. My first command in the session printed `a4ffd47`; `git rev-parse HEAD` now prints `285afba`; the working tree is clean. So the tree was checked out from under me mid-review — the `feedback_shared_tree_branch_switch` hazard, in the main tree rather than a worktree.
- **The specifics in the correction belong to a different task.** `tests/integration_ubam.rs`, `ws_qname_refused_in_trimming_loop_ubam_out_leaves_partial`, "§2.5 inventory and §5 line numbers" — none of those appear in the plan I was given (`mixed-format-pair-message/PLAN.md`, whose §2.5 is about specialty modes on mixed pairs). I never recorded a "line numbers are stale" finding, so there is nothing of that kind to withdraw.

## 2. The finding that supersedes the whole review

**The plan under review was implemented and merged twenty commits before HEAD.**

```
4276810 fix(cli): distinguish mixed FASTQ+BAM pairs from two-BAM pairs (#363) (#365)
```

`git log -S reject_bam_format_mismatch_in_pair -- src/format.rs` returns that commit and only that commit; it is an ancestor of `285afba`. The plan directory itself was committed by it. My review is a plan review of finished work, so **none of its action items are actionable as written** — which is why I am not asking anyone to work through them.

The reason the review nevertheless landed on real defects: the plan is written against `a4ffd47` (its own §2.1 says *"Reproduced on `dev` @ `a4ffd47`"*), so reviewing it against `a4ffd47` was internally consistent. That is also the one durable lesson here — a plan carrying a commit stamp in its body should be reviewed against the *current* trunk as well, because the gap between the two is where its line inventory rots.

## 3. Did the findings hold? Every one was independently resolved in the landed code

Verified by reading `src/format.rs:47-207`, the call site at `src/main.rs:370-375`, `tests/integration_paired_format_guard.rs`, `CHANGELOG.md:596-614`, and by running the binary.

| My finding | Status at `285afba` | Evidence |
|---|---|---|
| **C-1** clump-only + FASTQ output must keep the aux-tag diagnosis, not be told to interleave | **Resolved** | `PairedShape::ClumpOnlyFastqOut` (`format.rs:65-69`) with the doc *"suggesting `--paired interleaved.bam` would print a command that `main.rs`'s clump-only N=1 guard rejects"*; branch at `format.rs:165-172` |
| **I-1** predicate must be BAM-count, not format equality; plain+gz pair is legal; needs a label fn | **Resolved** | Named `reject_bam_format_mismatch_in_pair`; doc at `:97-99` states the BAM-count rule and the `FastqPlain`+`FastqGz` legality explicitly; `input_format_label()` at `:47-53` |
| **I-2** leaf guards / stale comments — silent-wrong-output risk | **Resolved, and better than I proposed** | An *enforcing* backstop at `main.rs:2285-2292`, with the reasoning I gave (*"the guard is ~1700 lines away, hence enforcing rather than documenting"*) plus a sharper one at `main.rs:1579` (*"a BGZF BAM handed to FastqReader parses as FASTQ rather than erroring"*). Five pointer comments, not the two the plan listed |
| **I-3** `samtools collate` cannot interleave two separate files | **Resolved** | Remediation is now `samtools merge -n -o interleaved.bam <r1> <r2>` (`format.rs:188-189`), confirmed in live output |
| **I-4** multi-pair partial output is an un-enumerated behaviour change | **Resolved** | `CHANGELOG.md:607-614` enumerates all three consequences, including multi-pair fail-before-writing; test `multi_pair_names_the_offending_pair_and_writes_nothing` |
| **I-5** `\`-continuation can split a pinned substring | **Resolved** | An explicit "Pinned-substring note" at `format.rs:176-181` naming all four phrases and the whitespace rule |
| **I-6** four validation gaps | **Resolved** | `tests/integration_paired_format_guard.rs`: `nonexistent_out()` helper (`:47-52`) encoding the pre-created-tmpdir gotcha; `mixed_pair_rejection_has_no_side_effects`; `!stderr.contains("uBAM paired mode expects")` at `:97`, `:231`, `:303`; `clump_only_fastq_output_keeps_aux_tag_diagnosis`; `plain_plus_gzip_fastq_pair_is_still_accepted` |
| **O-1** §2.1 over-generalises one site's stderr | **Resolved in the changelog prose** | `CHANGELOG.md:611-613` notes the uBAM-output and clump-only paths already rejected before detection |
| **O-2** helper belongs in `format.rs` | **Resolved** | That is where it lives |
| **O-3** redundant clump-only Shape B "must be BAM" branch | **Resolved** | Gone; only the general N=1 check at `main.rs:337` remains |
| **O-5** "Pair 1 of 1" is noise | **Resolved** | `pair_prefix` at `format.rs:154-158` emits `"Got "` for a single pair |
| **O-6** state the `inputs.len() == formats.len()` invariant | **Resolved, and corrected** | `anyhow::ensure!` at `format.rs:129-136` rather than the `debug_assert!` I suggested — with a comment explaining that release builds compile debug asserts out, so a truncating `zip` would make the guard silently pass. My suggestion was the weaker one |
| **O-7** fold format-dependent validation into `Cli::validate_formats(&[InputFormat])` | **Open** | Call site is still inline at `main.rs:370-375`. Optional structural suggestion, unchanged in force |

Live behaviour at `285afba`, three decisive cases:

```console
$ trim_galore --paired phred64_test.fastq ubam_test.bam -o OUT/nested
Error: --paired requires both inputs of a pair to be the same format. Got mixed:
…/phred64_test.fastq is FASTQ (plain) and …/ubam_test.bam is uBAM. Pass two FASTQ
files, or a single interleaved uBAM. If you meant two FASTQ files, check for a
mis-typed filename.
                                            → OUT/nested: ABSENT (no side effects)

$ trim_galore --clump_only --paired BS-seq_10K_R1.fastq.gz ubam_test.bam -o OUT
Error: uBAM input under --clump_only requires --output-format ubam (using the FASTQ
output path with uBAM input would drop aux tags). Input: …/ubam_test.bam

$ trim_galore --paired ubam_test.bam ub_copy.bam -o OUT
Error: --paired with two BAM files is not supported. uBAM paired mode expects a
single interleaved file: `trim_galore --paired interleaved.bam`. Got two BAM files:
… and … . Combine them into one mate-adjacent file first:
`samtools merge -n -o interleaved.bam … …`.
```

## 4. Artefact handling

`plans/mixed-format-pair-message/PLAN_review_reviewer-A.md` is **tracked at HEAD** (committed by `4276810`). After the checkout, `git diff HEAD` for that path is empty and its mtime reads `26 Jul 01:04`, so the repository's own copy is what is on disk and nothing of the project's was clobbered by my write. This re-verification went to a new filename rather than over it, per the fresh-filename-per-generation rule.

## 5. What I would actually do now

Nothing on #363 — it is done, and done well. The only carry-forward is **O-7**, which is a housekeeping suggestion rather than a defect: `main.rs` now hosts five format-dependent guards inline, and a `Cli::validate_formats(&self, &[InputFormat])` would put them where a maintainer looks for validation. Worth a line in a future housekeeping batch, not an issue of its own.
