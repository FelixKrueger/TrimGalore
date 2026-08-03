# PLAN v2 review — Reviewer B

## Verdict
**APPROVE WITH REVISIONS** — the shape is sound and the API references are accurate,
but four cross-section contradictions and one silent-regression need to be fixed
before implementation kicks off.

## Critical findings

### C1. Self-contradiction on `--fastqc` handling
- **What.** Three sections of the plan (Rejection matrix, lines 106-107; Signature note,
  line 184; Validation #8, line 404) all specify **warn-and-skip** for `--fastqc` under
  `--output-format ubam`. The **Resolved decisions §1** (line 413) locks the OPPOSITE:
  "run FastQC on the BAM directly. fastqc-rust natively supports .bam/.ubam/.sam."
- **Evidence.** Plan lines 106-107 vs 413; Signature line 184 says "The BAM path skips
  FastQC per Open Question 1"; Validation #8 (line 404) tests `no _fastqc.html produced`.
- **Fix.** Pick one and rewrite the other three sections consistently. If the resolved
  decision stands (run FastQC on BAM), then: (a) drop the "warn-and-skip" bullet from
  the Rejection matrix, (b) add `fastqc: bool` + `fastqc_args: Option<&str>` back to the
  new function signatures, (c) rewrite Validation #8 as a positive assertion (fastqc
  artifacts ARE produced). The Signature line 184 comment needs a full rewrite too.

### C2. Shape-A silently accepts two-BAM inputs (PE dispatch gap)
- **What.** Plan §3 dispatches PE by `inputs.len()`: N=1 → interleaved-BAM via
  `open_paired_interleaved_with_tags`; N=2 → two independent readers via
  `open_sync_reader`. When `inputs.len() == 2` and BOTH inputs happen to be uBAMs,
  the code will read them as two independent BAMs (assuming a synchronized ordering
  that need not exist), silently producing broken output.
- **Evidence.** Plan §3 impl (lines 258-278); the trim uBAM path REJECTS this same
  configuration explicitly at `main.rs:1502-1514` with a clear diagnostic pointing at
  `samtools collate` / single-interleaved input.
- **Fix.** Add the same rejection to `clump_only_paired_to_bam`: after detecting formats,
  if `inputs.len() == 2` and either input is `InputFormat::UnalignedBam`, bail with the
  trim path's message. Or upstream this into `Cli::validate` alongside the existing
  §3.4a checks so `--clump_only` inherits the same guard.

### C3. Multi-pair PE silently regresses under `--output-format ubam`
- **What.** v1 supports N=4+ FASTQ inputs (2+ pairs) via `run_specialty_paired`'s
  `chunks(2)` loop (`main.rs:333-353`). The plan's `clump_only_paired_to_bam(inputs:
  &[PathBuf])` matches `inputs.len() in {1, 2}` and bails on everything else
  (line 277-278). So `--clump_only --paired --output-format ubam R1a R2a R1b R2b` errors
  out where v1 accepted 2 pairs.
- **Evidence.** Plan §Signature (lines 165-172) + §3 impl (line 260-278) vs
  `main.rs:333` v1 path. Trim uBAM path handles multi-pair (`main.rs:1538` loop).
- **Fix.** Either (a) wrap the two-shape dispatch in a `chunks(2)` loop over `cli.input`,
  producing one interleaved BAM per pair (matches trim uBAM path); or (b) explicitly
  reject `input.len() > 2` at the CLI layer with a message pointing users at multiple
  invocations. Silently bailing during dispatch is the worst of both.

### C4. CLI validation edit is under-specified
- **What.** Plan §1 says: "Remove the current 'Input-format: uBAM input rejected'
  rejection." That rejection isn't in `cli.rs` — it lives in `clump_only.rs::reject_ubam`
  (line 233) invoked from `clump_only_single` / `clump_only_paired`. What IS in `cli.rs`
  is the `--paired && input.len() == 1` rejection at lines 701-707 that blocks Shape B.
- **Evidence.** `cli.rs:701-707` vs plan §1 bullet 1; `clump_only.rs:264,391-392` for the
  actual runtime `reject_ubam` calls.
- **Fix.** Rewrite §1 to specify: (a) remove `reject_ubam(input)?` from
  `clump_only_single`/`clump_only_paired` in `clump_only.rs`; (b) make the N=1 rejection
  at `cli.rs:701` conditional on `--output-format` (allow N=1 only if uBAM output, so
  Shape B is reachable); (c) remove the `--output-format ubam` deferral at
  `cli.rs:816-819`; (d) then add the new `--dont_gzip` rejection.

### C5. `--dont_gzip` rejection: false claim of "matches trim uBAM path"
- **What.** Plan §1 and Rejection matrix (lines 100-101) claim `--dont_gzip +
  --output-format ubam` rejection "matches how the trim uBAM path treats this."
  The trim uBAM path (`cli.rs::validate` §3.4a at `cli.rs:546-580`) does NOT reject
  `--dont_gzip` — only `--clumpify + --dont_gzip` is rejected (`cli.rs:650-652`).
- **Evidence.** Grep `dont_gzip` in cli.rs — no cross-check with `OutputFormat::UBam`.
- **Fix.** Either (a) drop the "matches trim" claim and justify on its own merits, or
  (b) hoist the new rejection into the shared §3.4a block so it applies to ALL uBAM
  output paths uniformly (recommended — the concept "BGZF is always compressed" is
  path-agnostic).

## Notable but non-blocking

- **N2. Mixed input formats in Shape A** (FASTQ R1 + uBAM R2, or vice versa). Plan §3
  peeks the header from `inputs[0]` only, so if R1 is FASTQ and R2 is uBAM, R2's `@PG`
  chain is silently dropped. Aux tags on R2 still fold into `FastqRecord.id` and thus
  reach the output — so provenance is inconsistent. Consider: warn, or require both
  sides same format, or preserve the union.
- **N3. Aligned-BAM per-record rejection** not explicitly mentioned. Inherited from
  `BamReader::next_record` but worth a one-line note in §Behavior for reviewers.
- **N4. `estimated_record_bytes` audit (§6) is a no-op** — `clump.rs:250-252` already
  counts `rec.id.len()` in full, so aux-tag tails are naturally included. Plan is
  correct to expect this; the "verify at implementation time" language is fine but the
  outcome is decidable now.
- **N5. Empty-input BAM invariant depends on unverified assumption** (§Assumptions
  bullet 6). The BamWriter test at `bam.rs:1641-1647` (empty writer path via no
  `write_record` calls) already covers this — cite it as evidence rather than deferring.
- **N6. Cross-run BAM determinism** (Validation #4) assumes noodles emits aux tags in a
  stable order across runs. Testing via `samtools view | sort` covers text-level order
  but if aux-tag emission order jitters, this could still be flaky. Low risk; flag for
  monitoring.
- **N7. `BamWriter::create`'s `_preserve_tags` param is unused** (`bam.rs:531`, doc says
  "kept for API symmetry"). Plan threads it through anyway; harmless but worth knowing.

## What the plan does well

Signatures accurately match the real code (`BamWriter::create` at `bam.rs:528`,
`peek_header` at `bam.rs:813`, `open_paired_interleaved_with_tags` at `bam.rs:146`,
`build_output_header` at `bam.rs:619`, `open_sync_reader` at `format.rs:89`). The
architectural placement — early dispatch from `main.rs::main`, sibling to
`hardtrim5_to_bam`, reusing v1's `SingleBin`/`PairedBin`/sort primitives — is exactly
right. The validation matrix (§V1-V9) covers the load-bearing invariants (record
byte-identity, pair lockstep, `@PG` chain, cross-run determinism, tag round-trip, PE
shape dispatch, empty input). The peak-memory bonus observation (BAM's streaming write
sidesteps v1's per-bin gzip Vec) is a real and correct win.
