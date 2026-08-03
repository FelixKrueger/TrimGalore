# CODE REVIEW — Reviewer B — #363 accurate per-pair format rejection

**Reviewed:** working tree on `fix/mixed-format-pair-message`, based on `dev` @ `a4ffd47` (nothing committed).
**Files:** `src/format.rs`, `src/main.rs`, `src/clump_only.rs`, `tests/integration_paired_format_guard.rs` (untracked), `CHANGELOG.md`.
**Plan:** `plans/mixed-format-pair-message/PLAN.md` v2, treated as a claim to verify.

## Summary

The change is sound. The core design — a single pure decision function over the already-materialised
`input_formats`, dispatched by a `PairedShape` enum, called once in `main()`'s cheap-validation
window — is correct, and every one of the specific hazards the review brief named checks out under
inspection and under the binary. In particular the `ClumpOnlyFastqOut` asymmetry (the Critical from
the plan's first draft) is implemented correctly and its message is textually identical to the one it
replaces; the predicate really is BAM *count*, not format equality; no pinned substring straddles a
line break; the retired dead check really was unreachable; and the validation-matrix code paths are
not touched at all.

I found one genuine defect (a stale cross-reference left pointing at code this change deleted) plus a
weakly-specific test assertion, and I fixed both. Everything else is a recommendation.

- `cargo test --release` — **441 passed, 0 failed** (before and after my edits)
- `cargo clippy --all-targets --release -- -D warnings` — clean
- `cargo fmt --all -- --check` — clean
- Manual smoke tests across 12 invocations — all as documented

**Process note for the lead, read this first.** Reviewer A is editing the same working tree
concurrently. One of my `Edit` calls returned *"the file had been modified on disk since you last
read it"*, and a second edit failed because A had already applied the identical fix. The CHANGELOG
has also been rewritten under me mid-review. Neither reviewer's "fixes applied" section is therefore
a complete account of the tree — verify the final `git diff` directly before merging.

---

## Verified sound

Stated here because a review that only lists complaints hides the fact that these were checked.

**Guard activation condition — correct, and no paired entry point is missed.**
I enumerated every paired dispatch site with `grep -n 'cli.paired\|chunks(2)' src/main.rs`:
`:259` (N=1 guard), `:286` (the new guard), `:464`/`:471` (clump-only FASTQ), `:535` (clump-only
uBAM), `:673` (paired uBAM single file), `:687` (paired FASTQ trim), `:1753`/`:1766`
(`run_ubam_output`), and `:2388`/`:2405` (`run_specialty_paired`, shared by `--clock`, `--implicon`
and clump-only-FASTQ-paired). Every one is either covered by the guard or deliberately exempt.

Fires as required (**ran the binary** for each): trim, trim + `--output-format ubam`,
`--clump_only` + FASTQ out, `--clump_only` + `--output-format ubam`.

Does *not* fire, behaviour unchanged (**ran the binary** for each):
- `--hardtrim3 10 --paired <fastq> <bam>` → exit 0, both outputs written. And it is genuinely
  harmless: `specialty::hardtrim5`/`hardtrim3` open via `crate::format::open_sync_reader`
  (`specialty.rs:34`), so the BAM side is read by `BamReader`, not misparsed. The test comment
  claiming this is accurate.
- `--clock --paired <fastq> <bam>` → unchanged `Paired-end files have different numbers of reads!`
- `--implicon=8 --paired <fastq> <bam>` → same unchanged error
- N=1 interleaved uBAM → accepted

Non-paired modes cannot reach it: `--demux` + `--paired` is rejected in `Cli::validate` at
`cli.rs:923` (*"Demultiplexing is only allowed for single-end files"*). `--passthrough` keeps its own
earlier, more specific rejection because `main.rs:253` precedes `main.rs:286` — **confirmed by
running** `--paired --passthrough <i1> <fq> <bam>`, which still reports *"--passthrough is not
supported with uBAM input in this release"*. `--retain_unpaired` has no format interaction and
correctly resolves to `PairedShape::Trim`.

**The `ClumpOnlyFastqOut` asymmetry — correct and faithful.**
`format.rs:141-149` returns before the `n_bam == 2` / `n_bam == 1` split, so the arm effectively
rejects on `n_bam >= 1`. The message is character-for-character the text at `clump_only.rs:404-407`,
and `first_bam` reproduces the original loop's *first-BAM-in-order* choice
(`clump_only.rs:401-410` iterated `[input_r1, input_r2]` and bailed on the first hit). **Ran all
three shapes** — `(fq, bam)`, `(bam, fq)`, `(bam, bam)` — and all three produce the identical
aux-tag message naming the correct file. The original `clump_only.rs` check is retained as
defence-in-depth for direct library callers, which is right.

**Predicate is BAM count, not format equality.** `pair_guard_accepts_plain_plus_gzip_fastq`
(`format.rs:392`) plus `plain_plus_gzip_fastq_pair_is_still_accepted`
(`integration_paired_format_guard.rs:363`). The integration test is *not* vacuous: it decompresses a
fixture to a plain `.fastq`, runs the pair, and asserts `out.join("plain_R1_val_1.fq").exists()`.

**Retired code really was unreachable.** The deleted check's condition was
`clump_only ∧ output_format==UBam ∧ paired ∧ N==1 ∧ ¬BAM`; `main.rs:259`'s is
`paired ∧ N==1 ∧ ¬BAM`. Strict subset, and `:259` precedes the clump-only dispatch at `:447` with no
intervening early return. Confirmed by reading both conditions.

**`chunks(2)` indexing is safe today.** `Cli::validate` calls `validate_paired_input("Paired-end")`
at `cli.rs:596` under a plain `if self.paired`, and `cli.rs:506` enforces an even count with only the
`N == 1 ∧ paired ∧ mode_label == "Paired-end"` carve-out at `:503` — which the guard excludes. So an
odd-length slice cannot reach `paths[1]`. See Low-3 for the residual hardening point.

**Backstop placement in `run_paired` is right.** `main.rs:1287-1296` sits after the banner but
*before* both reader-open branches — `open_threaded_reader` at `:1301` and the hard-coded
`FastqReader::open` at `:1333` — and before any `FastqWriter::create` at `:1336`. The
`--cores 1` default path is therefore covered, and nothing is written when it fires.

**Rust line-continuation hazard — clear.** All four pinned substrings render contiguously.
**Verified by running the binary and reading stderr**, not by re-reading the literals:
`--paired with two BAM files is not supported. uBAM paired mode expects a single interleaved file: …`
and `--paired requires both inputs of a pair to be the same format. Got mixed: …`. The three existing
tests pass unedited, and `git diff --stat tests/` was empty before my own test edits.

**Message register.** Plain, declarative, no filler, no second person, no scolding. Both cases read
as grammatical sentences: single-pair *"Got mixed: A is FASTQ (plain) and B is uBAM."*; multi-pair
*"Pair 2 of 2 is mixed: A is FASTQ (gzip) and B is uBAM."* The mixed message's *"Pass two FASTQ
files, or a single interleaved uBAM"* is valid for all three shapes that use it — I checked that a
two-FASTQ pair and a single interleaved uBAM are both legal under `Trim`, `TrimUbamOut` and
`ClumpOnlyUbamOut`.

**Validation-matrix invariant holds.** The SE, `--hardtrim5`, `--clock` and `--demux` code paths are
not touched by the diff at all. On the PE trim path the only change is the backstop's message string;
its control flow (bail iff any input is BAM) is identical, so an all-FASTQ pair takes exactly the
same route. `input_format_label` moved modules but its match arms are byte-identical. **Verified by
reading the complete diff**, plus 441 green tests — not by md5-comparing outputs; PLAN §13 claims an
md5 measurement was performed, which I did not re-run.

**Efficiency.** No hot path touched; one pass over a handful of `Copy` elements with early exit.
Net I/O on accepted paths is unchanged-to-lower: `run_ubam_output`'s deleted loop did 2
`detect_input_format` calls per pair and its leaf did 1 (3 total) versus 2 now; same arithmetic for
clump-only Shape A. Note PLAN §6's claim that *all three* deleted guards stop re-detecting is
slightly off — `run_paired`'s was kept, so the plain trim path still pays 2 detects per pair. That is
the correct trade (Low-2 below) and the cost is per pair, not per read.

---

## Issues by area

### Logic

Nothing wrong found. See Low-2 and Low-3 for two narrow hardening points.

### Errors

**M-1 (FIXED) — stale cross-reference pointing at code this change deleted.**
`src/main.rs`, formerly lines 1991-1993 in `run_ubam_output_paired_two_files`:

```rust
// Two-file paired-BAM rejection lives up-front in `run_ubam_output`
// (code-review round-2 B-NIT-1) — both inputs are guaranteed FASTQ
// here.
```

**Confirmed wrong, not suspected.** `git show HEAD:src/main.rs | grep -n B-NIT-1` returns two hits:
`:1749` (the loop this change deletes) and `:1983` (this comment referring to it). The loop is gone;
the comment still says the rejection *"lives up-front in `run_ubam_output`"*. Worse, it sat three
lines above the new backstop comment, which says the opposite — that the guard is *"in `main()` …
~1700 lines away"*. Two adjacent comments disagreeing about where an invariant is established is
exactly what PLAN §4 Step 6 set out to prevent; it listed four stale comments and this is a fifth
that was missed. Deleted, since the new backstop block below it states and now *enforces* the same
fact.

**M-2 (already fixed by Reviewer A concurrently) — internal-error message named the wrong function.**
`main.rs:1290` said `internal error: uBAM input reached run_paired_end (…)`. The enclosing function is
`main.rs::run_paired` (`:1238`); `run_paired_end` is `trimmer::run_paired_end` (`:1359`), a different
function. A maintainer grepping the reported name would land in `trimmer.rs` and not find the bail.
The other two backstops name their real functions (`clump_only_paired_to_bam_one_pair`,
`run_ubam_output_paired_two_files`), which is what identifies this as a slip rather than a choice —
PLAN §4 Step 4 proposed the wrong name verbatim, so the plan carried it in. Now reads `run_paired`.

**M-3 (already fixed by Reviewer A concurrently) — CHANGELOG attributed a wrong command to the
shipped message.** The entry claimed *"It suggested `samtools collate -O r1.bam r2.bam`"*. It never
did: `git show HEAD:src/main.rs | grep -n samtools` returns only `samtools import` (`:243`) and
`samtools fastq` (`:254`). The old two-BAM message (`git show HEAD:src/main.rs`, lines 1271-1275) ends
at *"Got two BAM files; one of them is {}."* with no remediation command at all. The `collate` form
existed only as a **proposal in the plan's v1 draft**, which was corrected before implementation — so
the CHANGELOG described a bug that was never shipped, in a user-facing file. Now reworded to say the
old message gave no combine command.

### Structure

Naming, placement and doc comments are good. `PairedShape` in `format.rs` is the right home:
`format.rs` owns `InputFormat`, has an existing `#[cfg(test)]` module, and `main.rs` has none.
Promoting `input_format_label` removes a duplicate rather than adding an abstraction.

### Test quality

Real assertions, not vacuous. The specific trap the plan warned about is avoided properly:
`nonexistent_out()` (`integration_paired_format_guard.rs:46`) asserts the path is absent as a
*precondition* before handing it to `-o`, so the post-run `!out.exists()` is a genuine statement
about `ensure_output_dir` never having run. I confirmed empirically that none of the rejected
invocations creates its output directory.

**L-1 (FIXED) — one assertion could pass without the guard having fired.**
`multi_pair_names_the_offending_pair_and_writes_nothing` asserted
`stderr.contains("Pair 2 of 2")`. The per-pair progress banner at `main.rs:754` prints
`=== Pair 2 of 2 ===`, which also satisfies that substring — so a regression in which the guard
stopped firing and the run *proceeded* to pair 2 would still match. (The companion `!out.exists()`
assertion does catch it in practice, since `ensure_output_dir` runs before dispatch, so this was
latent rather than live.) Tightened to `"Pair 2 of 2 is mixed"`, which only the guard can produce,
with a comment saying why.

**L-2 (FIXED) — test doc comment states a requirement the test does not meet.**
The same test's comment said *"Four **distinct** paths are required"*, but only three distinct paths
are passed — `BS-seq_10K_R2.fastq.gz` appears twice, as pair 1's R2 and pair 2's R1. PLAN §9 V8 has
the same wording, so the plan carried it in. The property that actually matters is distinct
*output* paths, since a collision would fail the test on precedence rather than on pair indexing.
Rewritten to say that.

---

## Fixes applied

All three re-verified afterwards: `cargo test --release` 441 passed / 0 failed,
`cargo clippy --all-targets --release -- -D warnings` clean, `cargo fmt --all -- --check` clean.

| # | File | Change |
|---|---|---|
| M-1 | `src/main.rs`, `run_ubam_output_paired_two_files` | Deleted the 3-line comment pointing at the `run_ubam_output` loop this change removed |
| L-1 | `tests/integration_paired_format_guard.rs:346` | `contains("Pair 2 of 2")` → `contains("Pair 2 of 2 is mixed")`, + comment on why |
| L-2 | `tests/integration_paired_format_guard.rs:328` | Rewrote the "four distinct paths" doc comment to state the real invariant (distinct outputs) |

Not mine, applied by Reviewer A during my review: M-2 and M-3 above.

---

## Recommendations

### Medium

**Med-1 — CHANGELOG generalises one path's behaviour to all four.** The bullet reads:

> The check also moved into input validation, ahead of dispatch. Three consequences on rejected
> runs, all improvements but all observable:
> - Adapter auto-detection no longer runs first. It previously scanned up to 1 M reads per input
>   before the guaranteed rejection.

True only of the plain trim path. **Verified by reading the call order:**

| Shape | Did adapter detection run before the old rejection? |
|---|---|
| `Trim` | **Yes** — `setup_trimming` at `main.rs:771`, then `run_paired` → old bail at `:1268` |
| `TrimUbamOut` | No — the deleted loop was at `run_ubam_output`'s top (`HEAD:1745-1761`), `setup_trimming` at `:1806` |
| `ClumpOnlyFastqOut` | No — `--clump_only` performs no adapter detection at all |
| `ClumpOnlyUbamOut` | No — same |

The third bullet (*"Multi-pair runs now fail before any pair is processed"*) splits the same way:
true for `Trim` and `ClumpOnlyFastqOut`, where the old check sat inside the per-pair loop; the two
uBAM-out paths already had up-front loops over all chunks. The middle bullet
(*"`--output_dir` is no longer created"*) **is** true for all four, since `ensure_output_dir`
(`main.rs:315`) precedes every dispatch branch.

This matters because PLAN §2.1 flags this precise generalisation as one of v1's errors — *"v1
generalised one site's stderr to both; corrected here"* — and it has come back in the CHANGELOG.
Minimal fix: scope the first and third bullets, e.g. *"On the ordinary trimming path, adapter
auto-detection no longer runs first…"* and *"Multi-pair trimming runs now fail before any pair is
processed…"*. Left as a recommendation rather than a fix because CHANGELOG wording is editorial.

### Low

**Low-1 — the two leaf backstops check mismatch only, narrower than the check they replace.**
`main.rs:2012` and `clump_only.rs:947` both test
`matches!(fmt_r1, UnalignedBam) != matches!(fmt_r2, UnalignedBam)`, so a **two-BAM** pair reaching
either leaf is not caught. The deleted `main.rs:536-558` caught both cases. The invariant the guard
actually establishes at these leaves is `n_bam == 0`, so `if n_bam != 0` would be the exact backstop.

Two arguments cut the other way, which is why this is Low and a recommendation: (a) at
`clump_only.rs:956` the code branches on `if matches!(fmt, UnalignedBam) { peek_header(r1_path) }`,
i.e. it carries latent two-BAM Shape A support that a stricter backstop would contradict; and (b) a
two-BAM pair would produce plausible output rather than the silently-mixed output the mixed case
produces. But note PLAN §4 Step 5's stated rationale — *"read the source header from R1 only"* —
applies to two-BAM too: R2's `@PG` chain would be silently dropped. Either widen the check or record
in the comment why two-BAM is deliberately out of its scope.

**Low-2 — `reject_bam_format_mismatch_in_pair` is `pub` and asserts only half of what `chunks(2)`
relies on.** `format.rs:112-116` `debug_assert_eq!`s `inputs.len() == formats.len()` (PLAN A7) but
not evenness, and `debug_assert` compiles out in release. An odd-length slice from a future caller
panics on `paths[1]` at `format.rs:139` rather than erroring — and `total_pairs = inputs.len() / 2`
would also under-report (`"Pair 2 of 1"`). Unreachable today (verified above), so this is hardening,
not a bug. Either add `debug_assert!(inputs.len().is_multiple_of(2))` next to the existing assert, or
destructure with `if let [r1, r2] = paths` so a short chunk is skipped rather than fatal.

**Low-3 — multi-pair `ClumpOnlyFastqOut` loses the pair index.** `pair_prefix` is computed at
`format.rs:130-134` but the `ClumpOnlyFastqOut` arm bails at `:142` before using it, so
`--clump_only --paired A.fq B.fq C.fq D.bam` names the offending *file* but not the pair. This
follows directly from PLAN §3.4's decision to preserve that message verbatim (pinned by V13), so it
is a consequence rather than a defect — but the information is cheap to add and the other three
shapes have it. Suggest appending `({pair_prefix}…)` when `total_pairs > 1`, or noting the omission
in the comment at `:137-140`.

**Low-4 — *"Pair 2 of 2 is two BAM files: A and B."*** Grammatical and unambiguous, but "contains
two BAM files" reads better and breaks no pinned substring — `two BAM files is not supported` is a
separate sentence, and the mixed-message test's `!contains("two BAM files")` applies only to the
mixed branch. Purely editorial.

**Low-5 — PLAN §13's deviation list is missing one entry.** §3.4 specified the mixed message as
*"Pass two FASTQ files, or a single interleaved uBAM`{fmt_hint}`."*; the implementation drops
`{fmt_hint}`. Harmless — both legal shapes are already named for every shape that uses the message,
and `{fmt_hint}` was never specified — but §13 lists four deviations and this is a fifth. Worth one
line if the plan is meant to stay an accurate record.

---

## Not findings, checked and dismissed

- **`is_bam(f)` with `f: &&InputFormat`** inside `.filter()` — compiles via deref coercion at the
  call site. No issue.
- **`pair_prefix` allocating a `String` unused on one arm** — an error path taken once per run.
- **Guard now preceding the three collision pre-flights** — a real precedence change, correctly
  identified and justified in PLAN §2.4, and `ubam_out_two_bam_pair_rejected` still passes.
- **Extra `detect_input_format` calls in the two leaves** — one file open plus one BGZF-block
  decompress per pair, on paths that then read the whole file. Net I/O still down (see above).
- **`flate2` used from an integration test** — regular `dependencies` are linked into test targets
  alongside `dev-dependencies`, and it compiles.
