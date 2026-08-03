# Reviewer B — Code Review, feature/clump-only-ubam (v2 uBAM add-on)

Scope: working-tree diff against `feature/clump-only`. Focus: correctness,
hidden assumptions, rejection matrix, resource paths, test quality.

## 1. Verdict

**NEEDS SIGNIFICANT REWORK** — one reproducible panic path plus a
significant test-coverage gap on a stated load-bearing invariant. The
first issue must be fixed before merge; the second should be fixed
alongside.

## 2. Critical findings

### C-1. Panic on `--clump_only --paired <uBAM>` without `--output-format ubam`

**(a) What's wrong.** `Cli::validate` used to reject `--clump_only + --paired + N=1` early (v1 line 705-716 in the pre-diff cli.rs). v2 deleted that guard (cli.rs:706-711, replaced with a comment claiming "dispatch bails with a clear error"). But dispatch only bails for the `OutputFormat::UBam` arm (main.rs:422-436). If the user forgets `--output-format ubam` and passes a single uBAM, `output_format` stays at its `Fastq` default and the flow reaches `run_specialty_paired`. That function unconditionally indexes `chunk[1]` at main.rs:2222 (`output_names(&chunk[0], &chunk[1])`) with `cli.input.chunks(2)` yielding a length-1 chunk — panic with `index out of bounds: the len is 1 but the index is 1`.

**(b) Reproduction.** Trace, all inputs verifiable in-tree:
- `Cli::validate_paired_input("Paired-end")` allows N=1 (cli.rs:503-505).
- The general N=1 guard at main.rs:215-221 skips when input[0] is `InputFormat::UnalignedBam` — so a single uBAM slips through.
- main.rs:354 matches `OutputFormat::Fastq` (default).
- main.rs:357 `if cli.paired` → true → `run_specialty_paired(...)`.
- main.rs:2222 `chunk[1]` on a length-1 chunk → panic.

**(c) Fix.** Restore an early-reject in `Cli::validate` or add a guard in the `OutputFormat::Fastq` arm of the clump_only dispatch (main.rs:355-414):

```rust
if cli.paired && cli.input.len() == 1 {
    anyhow::bail!(
        "--clump_only --paired with a single input requires \
         --output-format ubam (interleaved uBAM). Got FASTQ-output default. \
         Either add --output-format ubam or pass R1 and R2 as two files."
    );
}
```

## 3. High-severity findings

### H-1. Aux-tag round-trip has zero end-to-end test coverage

**(a) What's wrong.** PLAN, CHANGELOG, and docs/clump-only.md all promise
byte-identity of *"id + seq + qual + preserved aux tags"* through the
FASTQ intermediate. The integration file `tests/integration_clump_only_ubam.rs`
never passes `--preserve-tags`. `se_ubam_in_ubam_out_pg_chain` (line 190-191)
does `assert_eq!(sorted_bam_tuples(&input), sorted_bam_tuples(&out))` where
`bam_tuples` reads aux tags (line 54-63) — but `test_files/ubam_test.bam`
contains **no aux tags on any record** (verified via
`samtools view test_files/ubam_test.bam` — columns end at $11). So the
assertion passes trivially even if `--clump_only` silently dropped
every tag it saw.

**(b) Evidence.**
- `tests/integration_clump_only_ubam.rs`: no `--preserve-tags` in any test.
- The trim uBAM path uses `test_files/ubam_test_with_tags.bam` for exactly this coverage (`tests/integration_ubam_out.rs:236, 340, 394`). The clump path does not.
- The three new CI validation steps (`.github/workflows/ci.yml:9-61`) also use `awk '{print $1,$2,$10,$11}'` — dropping aux tags.

**(c) Fix.** Add at minimum one integration test that:
1. Runs `--clump_only --output-format ubam --preserve-tags CB,UB` on `ubam_test_with_tags.bam`.
2. Asserts every input tag appears verbatim on the correctly re-associated output record.

Also strengthen the CI `record parity` steps to include tag preservation when a tag-carrying fixture is available.

## 4. Medium / notable non-blocking

### N-1. Dead code in Shape A BAM branch

`clump_only_paired_to_bam_one_pair` line 771-776 (`clump_only.rs`) checks
`if matches!(fmt, InputFormat::UnalignedBam) { Some(peek_header(...)) }`
for Shape A, but main.rs:459-481 rejects any BAM in a two-file pair.
This branch is unreachable. Not wrong — just confusing. Consider dropping
or turning into a `debug_assert!`.

### N-2. `preflight_collision_bam` on 1-path slice is a no-op

main.rs:443 calls `preflight_collision_bam(&[planned])` for Shape B (single
interleaved uBAM). One element can't collide with itself. Harmless; consider
skipping to make intent clear.

### N-3. Truncated-BAM risk on mid-stream error

`clump_only_single_to_bam` and `clump_only_paired_to_bam_one_pair` propagate
`?` from `next_record()`, `flush_bin_*_to_bam(...)`, and
`writer.write_record(...)`. `BamWriter` does not implement `Drop`, so a
mid-stream error skips `writer.finish()` and produces a BAM without the
BGZF EOF marker. Standard Rust pattern (caller discards the output on
error), but downstream tools that tolerate missing EOF (samtools with the
`ignore-truncation` fallback) could silently see a partial reorder.
Matches the existing trim uBAM path behaviour; note for the follow-up.

### N-4. Per-input sanity check gap on multi-input

`sanity_check_any` runs only on `cli.input[0]` (main.rs:176). Multi-input SE
and multi-pair PE clump_only runs proceed without sanity checking pairs 2+.
Pre-existing pattern from v1 FASTQ; the trim uBAM path DOES sanity-check
each pair (main.rs:1690-1692). Consider parity.

### N-5. `preserved_tags` populated on FASTQ input

`ClumpOnlyStats.preserved_tags` (clump_only.rs:619, 831) is unconditionally
set from `&cli.preserve_tags`, so a `--clump_only --output-format ubam
--preserve-tags CB,UB fastq.gz` run would print `Preserved tags: CB,UB` in
the report even though no tag ever existed. The main.rs:198-203 guard
should catch all-FASTQ + uBAM output + preserve-tags, so this is unreachable
in practice — but the coupling is worth documenting or hardening.

## 5. What the implementation does well

The dispatch topology in main.rs is well-organised: format-guards run
before any reader is opened, collision pre-flight is unified via
`preflight_collision_bam`, and the four sub-shapes (SE-BAM, PE-Shape-A,
PE-Shape-B, and the existing FASTQ variants) all reach the same
`clump_only_*_to_bam` entry point. `PairedInputSetup` struct is a clean
sidestep of the type-complexity clippy warning. The
`stats_shape_from` extension and the report label fallbacks preserve v1
byte-identity of the FASTQ report shape while cleanly adding the uBAM
variant. `@PG` chain preservation delegates fully to `bam::build_output_header`
(existing, reviewed), keeping the load-bearing invariant on the trim
uBAM path — good code reuse.

---

Review file: `/Users/fkrueger/Github/TrimGalore/plans/clump-only/CODE_v2_review_reviewer-B.md`
