# v2 uBAM code review — Reviewer A

## Verdict
**APPROVE WITH REVISIONS.** No correctness bugs that would produce wrong output on the happy path. The invariant-critical code (flag bits, mate-adjacent write order, `@PG` chain, aux-tag round-trip) is sound and reuses the trim-path primitives verbatim. Blockers are all on the test/observability side: several planned tests were not written, and one that was written does not exercise the intended code path.

## Critical findings

### C-1. `rejects_two_bam_paired` never reaches the two-BAM Shape A guard
`tests/integration_clump_only_ubam.rs:322-347` passes the same `ubam_test.bam` twice. `Cli::validate` → `validate_paired_input("Paired-end")` (src/cli.rs:513-521) rejects the R1==R2 duplicate first, so the two-BAM Shape A rejection at `src/main.rs:456-467` is never executed by the suite. The test's own comment concedes this ("Cli::validate detects the R1==R2 case first"). Result: the main two-BAM guard is dead-code-coverage-wise.

Fix: copy the fixture under a second name inside the temp dir and pass R1=orig, R2=copy — as the multi-pair test already does at lines 391-392.

### C-2. Planned tests missing (§7 of the plan)
None of the 4 unit tests the plan committed to (`test_clump_only_single_to_bam_permutation`, `test_clump_only_paired_to_bam_lockstep`, `test_clump_only_bam_deterministic_records`, `test_clump_only_ubam_in_ubam_out_tag_roundtrip`) appear in `src/clump_only.rs` (grep returns 11 test fns, all v1). Integration tests cover permutation and interleave shape, but **aux-tag round-trip through the FASTQ intermediate is exercised nowhere** — the `sorted_bam_tuples` equality in `se_ubam_in_ubam_out_pg_chain` compares aux via `format!("{:?}", value)` on both sides read the same way, which does not prove the write-path re-encoded them correctly (they could both be missing and still match). This is the load-bearing v2 invariant (A/Z/i/f preserved, B/H rejected) and it has no dedicated assertion.

Also missing per plan §7: `rejects_mixed_format_paired` (guard at `src/main.rs:468-473` untested) and `pe_bam_collision_preflight_case_folded` (`preflight_collision_bam` untested).

### C-3. Multi-pair PE-BAM regresses observability vs. v1
The new PE-BAM dispatch (`src/main.rs:479-521`) calls `clump_only_paired_to_bam_one_pair` in a bare `for chunk in cli.input.chunks(2)` loop with no per-pair `eprintln!("=== ... pair N of M ===")` and no `.with_context(||"processing pair N of M …")`. The sibling paths (v1 FASTQ via `run_specialty_paired` at `src/main.rs:2237-2256`, and the trim uBAM path at `src/main.rs:1683-1713`) do both. On failure at pair 3 of 5 the user now sees just `Failed to create BAM output: /tmp/.../foo_clumped.bam` with no indication which pair. This is a real UX regression from v1 shape.

### C-4. PE-BAM multi-pair skips per-pair `sanity_check_any`
Only `cli.input[0]` is sanity-checked (`src/main.rs:176`). The trim uBAM PE path pointedly re-runs `sanity_check_any(&chunk[0])` for pairs > 0 and `sanity_check_any(&chunk[1])` for every pair (`src/main.rs:1689-1692`). The clump uBAM PE path skips both. An empty BAM at pair 2 will produce an empty header-only BAM output instead of the intended "uBAM with no records" bail.

### C-5. Misleading comment on empty-record handling
`src/clump_only.rs` (in `clump_only_single_to_bam`, comment above `writer.finish()`) says "BamWriter::finish handles zero records (verified during plan-review at bam.rs:1641-1647); no explicit empty-record path needed here". `bam.rs:1641-1671` (`bam_writer_synthesises_header_for_fastq_input`) writes ONE record — it does not test zero records. There is no `writer.finish()`-with-zero-records assertion anywhere in the tree. Either add a real empty-input test or remove the false attribution.

## Notable but non-blocking

### N-1. `has_trim_galore_pg` matches only key literal `b"trim_galore"`
`tests/integration_clump_only_ubam.rs:97`. If the input BAM already contains a `trim_galore` `@PG` (re-processing a clumped BAM), `Programs::add` disambiguates to `trim_galore-trim_galore` and the assertion fails even though provenance is correctly preserved. Latent fragility — no test currently triggers it, but worth matching by prefix.

### N-2. `--dont_gzip + --output-format ubam` rejection is a behavior change for the trim uBAM path
`src/cli.rs:547-553` hoists the reject to the shared §3.4a block, which the diff comment acknowledges: "the trim uBAM path (which previously accepted this silently — a real gap)". Correct call, but scripts relying on the old silent-accept will now error. Confirm this is mentioned in the CHANGELOG entry (not reviewed here beyond noting).

### N-3. `input_bytes` on Shape A is R1+R2 combined, printed against R1's path in the report
Same shape as v1's paired-FASTQ report. Report line `Input: R1.bam (uBAM, 47001200 bytes)` where the number is R1+R2. Minor semantic quirk; matches v1 for consistency.

### N-4. `input_format_label` returns `"uBAM"` for input vs `"uBAM (BGZF)"` for output
Aesthetic asymmetry — both are BGZF. Fine.

### N-5. `detect_input_format` runs 2× per input on the PE-BAM Shape A path
Once in the guard loop (`src/main.rs:461-462`), once in `clump_only_paired_to_bam_one_pair` (via `open_sync_reader`). Trivial cost; only note if refactoring later.

## What the implementation does well
The v2 diff cleanly bolts the uBAM path onto the v1 bin dispatcher: `flush_bin_single_to_bam` / `flush_bin_paired_to_bam` reuse `sort_single_by_key` / `sort_paired_by_key` unchanged, so v1's byte-identity invariant carries through verbatim. Mate-adjacent output is enforced by construction (single `sort_paired_by_key` then interleaved `write_record(_, Some(1|2))`). Format guards for two-BAM Shape A / non-BAM Shape B / mixed-format Shape A all live BEFORE `preflight_collision_bam`, which itself lives BEFORE any reader opens — so bad-input classes fail early without side-effects. The `PairedInputSetup` extraction sidesteps the `type_complexity` clippy warning without hiding intent. Cross-family compression-ratio (FASTQ.gz → BAM.bgzf) is correctly gated on `input_compressed && output_compressed && output_bytes > 0` per §Resolved decision 2. Report writer's fallback (`if input_format_label.is_empty() { … }`) preserves v1 on-disk report byte-identity on the FASTQ path — nice.
