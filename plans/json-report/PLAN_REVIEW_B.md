# Plan Review B — JSON Trimming Report

## Verdict: REVISE

The plan is well-structured, correctly identifies the right integration points, and the schema design is mostly sound. However, there are several issues ranging from a critical type mismatch that will cause a compile error, to missing fields that exist in the current text report, to gaps in the manual JSON writing strategy that could produce invalid output. These must be addressed before implementation.

---

## Critical Issues

### 1. Type mismatch: `write_json_report` takes `config: &TrimConfig` but there are TWO `TrimConfig` types

The plan specifies the function signature as:
```rust
pub fn write_json_report<W: Write>(
    w: &mut W,
    config: &TrimConfig,
    stats: &TrimStats,
    ...
)
```

But the codebase has **two separate** `TrimConfig` structs:
- `report::TrimConfig` (src/report.rs, line 104) — holds string-typed fields (`adapter: String`, `version: String`, etc.), used for report writing.
- `trimmer::TrimConfig` (src/trimmer.rs, line 14) — holds byte-vec fields (`adapter: Vec<u8>`) and operational fields (`clip_r1`, `clip_r2`, `max_n`, etc.), used for trimming.

The plan must explicitly state which `TrimConfig` it uses. Since this function lives in `report.rs`, it should take `report::TrimConfig`, which is what `main.rs` already constructs for the text report (lines 324-344 and 496-516). This is likely the intent, but it needs to be stated clearly because the implementer could easily reach for `trimmer::TrimConfig` and then struggle with type mismatches.

### 2. `--no_report_file` interaction is not specified

The plan says "The JSON report is always emitted (no opt-in flag needed)" (line 206-208). But `--no_report_file` (src/cli.rs, line 122-124) currently suppresses the text report. The `if !cli.no_report_file` guard on main.rs lines 321 and 491 wraps the entire report-writing block, including where the JSON write call would be inserted.

The plan must decide: does `--no_report_file` suppress the JSON report too? Almost certainly yes — a user who says "no reports" expects no report files of any kind. But the plan is silent on this. If the JSON write call goes inside the existing `if !cli.no_report_file` block (which is the natural placement), it inherits this behavior by default. This should be explicitly stated as a design decision.

### 3. JSON string escaping for filenames and adapter sequences is not addressed

The plan acknowledges manual JSON writing but provides no strategy for string escaping. JSON requires escaping of `\`, `"`, control characters, and certain Unicode codepoints. Real-world filenames can contain:
- Backslashes (Windows paths, though rare in bioinformatics)
- Double quotes (pathological but legal)
- Non-ASCII characters (international users)

Adapter sequences are safe (DNA alphabet only), but `command_line` (line 65 of the schema) can contain arbitrary user input — paths with spaces, quotes, backslashes, etc. The `input_filename` field can also contain special characters.

The plan must require a `json_escape_string()` helper function that handles at minimum: `\` -> `\\`, `"` -> `\"`, newline -> `\n`, tab -> `\t`, and other control characters. Without this, the JSON output will be invalid for some users and silently break MultiQC parsing.

---

## Recommendations

### 4. Missing field: `discard_untrimmed` is not in `parameters`

The `TrimStats` struct has `discarded_untrimmed` (src/report.rs, line 43), and the text report conditionally displays it (lines 184-187, 293-296). The JSON schema includes it in `read_processing.reads_discarded_untrimmed` (line 92), which is good. However, the **triggering parameter** `discard_untrimmed: bool` is not listed in the `parameters` object. A consumer seeing `reads_discarded_untrimmed: 500` cannot tell whether the user explicitly requested this behavior. Add `"discard_untrimmed": false` to the `parameters` section.

### 5. Missing fields from `report::TrimConfig` that are not in JSON `parameters`

The `report::TrimConfig` (src/report.rs, lines 104-126) includes `gzip: bool` which is not in the JSON schema's `parameters` section. While `gzip` is less relevant for MultiQC, for schema completeness (and since the plan says "has everything"), it should be included or explicitly noted as omitted-by-design.

### 6. `adapter_r2` in schema does not match the real data structure

The JSON schema shows `"adapter_r2": "AGATCGGAAGAGC"` as a string (line 69). But in `report::TrimConfig`, `adapter_r2` is `Option<String>` (line 109). For SE mode or when no R2 adapter is explicitly set, this is `None`. The schema should specify that this field is `null` when no R2 adapter is configured, not just show it populated.

### 7. `adapter_trimming.type` is hardcoded as "regular 3'" but the codebase supports more

The schema shows `"type": "regular 3'"` (line 102). Looking at the text report's Cutadapt section (src/report.rs, line 317), the type is also hardcoded to `"regular 3'"`. This is consistent, so the hardcoding is fine — just confirming this is intentional and not a gap.

### 8. `pair_validation` placement: plan says "R2 only" but this creates a consumer trap

The plan says `pair_validation: null` for R1/SE (line 155), and that it's populated for R2 in paired-end mode (line 132-133). This mirrors the Perl behavior where pair validation stats appear only in the R2 text report (confirmed in main.rs, line 524: `if idx == 1`).

However, this means a consumer parsing only R1's JSON report in paired-end mode would see `pair_validation: null` and have no way to get pair stats without also finding the R2 report. The plan should at least document this behavior clearly for MultiQC implementers. Consider: should `pair_validation` be written into BOTH reports in paired-end mode? MultiQC might want it from either file.

### 9. `read_number` for single-end should be documented

The plan says `read_number: 1` for SE (line 158). This is reasonable but should be made explicit: the implementer needs to know to pass `read_number = 1` in `run_single_file()` and `read_number = 1` / `read_number = 2` in `run_paired()`.

### 10. `--basename` interaction with JSON report naming is not discussed

When `--basename` is used (src/cli.rs, line 116), the trimmed output filenames change (e.g., `custom_R1_val_1.fq.gz` instead of `sample_R1_val_1.fq.gz`). But `report_name()` in src/io.rs (line 96-109) always uses the **input filename**, not the basename. So the JSON report naming should follow the same pattern — `{input_basename}_trimming_report.json` where `input_basename` is the original input filename, not the `--basename` override. The plan's naming convention (line 47-48) implies this, but it should explicitly state that `--basename` does NOT affect report filenames (matching current text report behavior).

### 11. `--output_dir` interaction needs explicit mention

The text report already respects `--output_dir` via `naming::report_name()` (src/io.rs, line 96-109). The JSON report should use the same logic. The plan should state: "use the same `naming::report_name()` pattern, but with `.json` extension instead of `.txt`." Consider adding a `json_report_name()` function to `src/io.rs`, or parameterizing `report_name()` with an extension.

### 12. `length_distribution` keys start at 1, but `adapter_length_counts` is 0-indexed

The `TrimStats.adapter_length_counts` Vec (src/report.rs, line 29) is indexed by match length, with index 0 being unused (the text report skips it: `if length == 0 || count == 0 { continue; }` at line 331). The JSON schema correctly shows keys starting at "1" (line 107). The plan should note this: skip index 0 and skip entries with count 0 (sparse representation).

### 13. "No serde" approach is reasonable but needs clear implementation guidance

The plan correctly identifies that the JSON is shallow (max 2 levels deep) and serde would be overkill. However, the two suggested approaches — "helper functions" vs "build the JSON as a String" — should be narrowed to one. I recommend the helper-function approach with `write!` macros and a `BufWriter`, matching the existing text report pattern. Building the entire JSON as a `String` first wastes memory for large length distributions and doesn't match the existing code style.

Also: the helpers need to handle trailing comma management in JSON. The last field in each object must NOT have a trailing comma. This is a common source of bugs in manual JSON writing. The plan should specify a strategy (e.g., a `json_field` helper that takes a `first: &mut bool` flag, or building fields into a Vec and joining with commas).

### 14. Test plan needs more specificity

The test plan (lines 228-235) has four items but is vague:
- "Unit tests: produces valid JSON" — how? The codebase has no JSON parser dependency. The test would need to either add `serde_json` as a dev-dependency for validation, or use a manual validation approach. The plan should specify.
- "Integration: Run TrimGalore on test data" — which test files? The `test_files/` directory has specific files; the plan should reference them.
- "Validation: Parse the JSON output with Python/jq" — this belongs in a CI script or manual testing guide, not as a unit/integration test. Clarify where this lives.

Missing test cases:
- Test with special characters in filenames (spaces, quotes)
- Test that `--no_report_file` suppresses JSON too
- Test with `--output_dir` to verify JSON lands in the right directory
- Test paired-end to verify R1 has `pair_validation: null` and R2 has it populated
- Test that `length_distribution` is sparse (zero counts omitted)
- Test with all-zero stats (empty input file edge case)

---

## Observations

### 15. The schema is well-designed for forward compatibility

`schema_version: 1` is a good practice. Sections with zero values always present (line 165-166) is the right call — it prevents consumers from needing to handle both "key missing" and "key present with zeros" cases.

### 16. Explicit `null` for `pair_validation` in SE mode is correct

This is better than omitting the key. It lets consumers use a single schema definition without optional fields.

### 17. No risk to existing text reports

The plan correctly identifies that no changes are needed to the text report (lines 215-223). The JSON is a new file alongside the existing `.txt`. The MultiQC fallback strategy (prefer JSON, fall back to text) is clean.

### 18. `basepair_processing` section is missing `adapter_trimmed_bp`

The text report does not explicitly track "basepairs removed by adapter trimming" as a separate field, so the JSON correctly omits it. However, this means `total_bp_processed - quality_trimmed_bp - total_bp_written` is NOT equal to adapter-trimmed bp (it also includes length-filtered reads). This is fine but worth noting for MultiQC implementers: the bp accounting does not sum to total_bp_processed when reads are discarded entirely.

### 19. Disk space impact is negligible

A JSON report for a typical run will be 1-3 KB. Even with length distribution, it will rarely exceed 5 KB. The plan's assessment that it's "small" (line 207) is correct.

### 20. The plan does not mention `merge()` on TrimStats for the parallel path

In `src/parallel.rs`, worker threads produce per-batch `TrimStats` that get merged via `TrimStats::merge()` (src/report.rs, line 48). The `adapter_length_counts` merge (lines 59-64) is correct — element-wise addition. The JSON report is written from the merged stats, so this should work transparently. No action needed, just confirming correctness.

### 21. Consider adding `output_filename` to the JSON

The text report does not include the output filename, but the JSON could. Knowing that `sample_R1.fastq.gz` was trimmed to `sample_R1_val_1.fq.gz` would help MultiQC link reports to output files. This is optional and could be a follow-up.
