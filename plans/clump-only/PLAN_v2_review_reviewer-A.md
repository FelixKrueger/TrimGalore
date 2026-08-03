# Reviewer A — PLAN_v2_ubam.md review

**Verdict:** APPROVE WITH REVISIONS

The plan is well-scaffolded and its verified references check out (`BamWriter::create` @ src/bam.rs:528 with 4 params matches; `BamReader::open_paired_interleaved_with_tags` @ src/bam.rs:146 exists; `format::open_sync_reader` returns `Box<dyn RecordSource>` @ src/format.rs:89-101; `estimated_record_bytes` counts `rec.id.len()` in full @ src/clump.rs:250-252; `hardtrim5_to_bam` signature @ src/specialty.rs:106 matches). But there are four blocking contradictions/gaps and a handful of smaller issues.

## Critical findings

### C1. Direct self-contradiction on `--fastqc` behavior
(a) §Rejection matrix / §Signature / §Validation-test-8 say `--fastqc` + `--output-format ubam` → **warn-and-skip**, no fastqc/fastqc_args params.
(b) §Resolved decisions #1 says run FastQC **directly on the BAM** (fastqc-rust 1.0.1 accepts `--format bam`), no warn-and-skip.

These are incompatible. If OQ1 is truly "locked" to option (b), the function signatures at lines 152-173 must add `fastqc: bool` + `fastqc_args: Option<&str>`, the rejection matrix line at 107 must be removed, validation test #8 at line 404 must be inverted (assert `_fastqc.html` **is** produced, no warning). Pick one and align the whole plan.

### C2. PE dispatch for N=1 with non-BAM input is silently broken
Current v1 CLI rejects `--clump_only --paired` with `input.len()==1` at src/cli.rs:701 with a precise message. Plan §Implementation-step-1 removes this rejection; §Implementation-step-3 unconditionally routes N=1 to `BamReader::open_paired_interleaved_with_tags(&inputs[0], …)`. If the user passes `--clump_only --paired one.fastq.gz`, this now fails inside `noodles::bam` with a cryptic BAM-magic error rather than a clean "you need two FASTQ files or one interleaved uBAM" message. Fix: format-detect `inputs[0]` before the N=1 branch and bail with a clear message if not `UnalignedBam`. Also consider whether `Cli::validate_paired_input("Paired-end")` at src/cli.rs:503-504 (which already lets N=1 through for `--paired`) is the right place to gate this.

### C3. False claim about `--dont_gzip` under trim uBAM path
Plan §Rejection matrix at line 100 asserts that rejecting `--dont_gzip` + `--output-format ubam` "matches how the trim uBAM path treats this". It does not — src/cli.rs:543-581 (§3.4a) accepts `--dont_gzip` silently on the trim uBAM path. Two options: (i) drop the parity claim and add the rejection clump-locally; (ii) better — add the rejection to §3.4a for all `OutputFormat::UBam` invocations (correct fix, catches an actual latent hole in the trim path too). The plan should be explicit about which one it picks; leaving parity claims wrong is a hazard once someone reads the code.

### C4. PE-BAM path bypasses `run_specialty_paired`'s collision pre-flight
Plan §Implementation-step-2 dispatches PE-BAM directly (`clump_only::clump_only_paired_to_bam(&cli.input, …)`), bypassing `run_specialty_paired` (src/main.rs:2061). The plan then claims: "PE-BAM's pre-flight uses `run_specialty_paired`'s existing collision check via the naming closure" — but `run_specialty_paired` is never called on this path, so this claim is false. `run_specialty_paired` also uses `cli.input.chunks(2)` (line 2076), which is structurally incompatible with Shape B (N=1). Fix: write an explicit pre-flight loop inside `clump_only_paired_to_bam` (or a wrapper) that computes the one output path per chunk (Shape A) or the single output path (Shape B), normalises via `naming::norm_path`, and rejects collisions before opening any writer. Same for the multi-SE-BAM loop (§Implementation-step-2 shows a bare `for input in &cli.input` — the v1 pre-flight loop must be preserved and switched to `clumped_bam_output_name` under `OutputFormat::UBam`).

## Notable but non-blocking

- **Mixed-format PE Shape A** (FASTQ R1 + uBAM R2, or vice versa): plan §Implementation-step-3 only peeks a header from `inputs[0]` and drops R2's `@PG` chain silently. Either reject the mixed case explicitly or document the "R1 wins for header" rule in Behavior.
- **`--cores` on BAM output**: CLAUDE.md states uBAM output is always single-threaded, `--cores` silently ignored. Plan's signatures accept `cores: usize` without documenting the silent-ignore semantics. Add a one-line note in §Signature.
- **Shape B report filename**: For single-interleaved-BAM PE input, "the input stem" for `<stem>_clumping_report.txt` isn't defined. Specify (e.g., strip `.bam` from `inputs[0]`).
- **`clumped_paired_bam_output_name(input_r1, _input_r2: Option<&Path>, …)`**: the second param is `Option`, but §Implementation-step-2's dispatch code doesn't call it — the plan should show how Shape B (N=1, no r2) versus Shape A wires up the caller.
- **Cross-run body-byte-identity claim**: identical body bytes require identical *input file bytes AND* identical sort tiebreaks. The claim is fine but should note the sort is stable via content cascade at src/clump.rs:271-280 — that's what makes it hold.
- **Test #9 "empty-input BAM"** is claimed to be a unit test but requires exercising `format::open_sync_reader` on a real BAM file; may need a fixture, not pure synthesis. Non-blocking.
- **Assumption #5 dead flag**: §Assumptions says "`--paired` + N=1 with uBAM input is legal (existing carve-out in `Cli::validate::validate_paired_input`)" — true (src/cli.rs:503), but the current `--clump_only` block at src/cli.rs:701 *overrides* that carve-out. Plan step 1 must explicitly remove line 701's gate, which it does implicitly by saying "Remove the current 'Input-format: uBAM input rejected' rejection" — but that phrasing doesn't match what's actually in the code (there is no such rejection in CLI; the uBAM rejection lives in `clump_only.rs::reject_ubam`). The plan should name the two exact removals: `clump_only.rs::reject_ubam` calls (lines 264, 391-392) and the `cli.rs:701` N=1 rejection and the `cli.rs:816-820` UBam-output rejection.

## What the plan does well

Verified references are precise and correct against source; template selection (`hardtrim5_to_bam`) is the right one; v1 primitive reuse (SingleBin/PairedBin/sort/resolve_layout) preserves the invariants Phase-1 already secured; the Behavior section's contract carefully separates record-body byte-identity from whole-file byte-identity via the `@PG` carve-out, matching CLAUDE.md's stated pattern for the trim uBAM path; validation section names concrete failure modes for each check.
