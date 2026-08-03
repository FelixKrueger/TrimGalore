# Plan Review A — JSON Trimming Report

## Verdict: REVISE

The plan is well-structured and makes sound high-level design decisions (schema versioning, sparse length distribution, always-present sections). However, there are several gaps in field coverage, ambiguities in the proposed function signature vs. actual types, and missing consideration of important CLI interactions. These need to be addressed before implementation.

---

## Critical Issues

### 1. Type mismatch: `TrimConfig` in plan vs. actual codebase

The plan's proposed function signature references `report::TrimConfig` (the report-layer config struct), but the JSON schema includes fields that do not exist on `report::TrimConfig`:

- **`mode`** (line 65 of plan: `"paired-end"`): Not a field on `report::TrimConfig`. Must be derived from `config.paired`.
- **`read_number`** (line 66): Not available on any config struct — this is a caller-supplied parameter. The plan correctly shows it as a function parameter (`read_number: u8`), so this is fine.
- **`command_line`** (line 67): Available on `report::TrimConfig.command_line` — OK.

This is not a blocker, but the plan should explicitly note that `mode` is derived from `config.paired` to avoid ambiguity during implementation.

### 2. Missing fields from `TrimStats` that exist in the text report

The JSON schema omits `discarded_untrimmed` from the `read_processing` section. The text report includes "Reads discarded as untrimmed" (report.rs line 184-187, and the Cutadapt section line 293-296). The field exists on `TrimStats` (line 43: `pub discarded_untrimmed: usize`).

The plan's JSON schema (line 92) lists `reads_discarded_untrimmed: 0` — wait, it IS there. Let me re-check... Yes, it is present at line 92 of the plan. This is correct. I retract this point.

### 3. Missing parameters from the JSON `parameters` section

The following trimming parameters exist in `trimmer::TrimConfig` and affect the output but are absent from the plan's `parameters` object:

- **`clip_r1`**, **`clip_r2`**, **`three_prime_clip_r1`**, **`three_prime_clip_r2`** (hard-clipping offsets): These are important for reproducibility. A consumer looking at the JSON should be able to reconstruct what clipping was applied. They are in `trimmer::TrimConfig` lines 25-28.
- **`max_n`** (N-filtering threshold): This parameter controls the `too_many_n` filter. Without it in the parameters section, a consumer cannot determine whether N-filtering was applied or what the threshold was. In `trimmer::TrimConfig` line 23.
- **`discard_untrimmed`** (boolean): Controls whether reads without adapters are discarded. In `trimmer::TrimConfig` line 36. While the count is in `read_processing`, the parameter itself should be in `parameters` for reproducibility.
- **`consider_already_trimmed`** (threshold): This parameter can cause adapter trimming to be skipped entirely. Not capturing it means a consumer cannot distinguish "no adapter trimming happened because reads were already trimmed" from "no adapter trimming happened because no adapters were detected."

**Recommendation:** Add these to `parameters`:
```json
"clip_r1": null,
"clip_r2": 2,
"three_prime_clip_r1": null,
"three_prime_clip_r2": null,
"max_n": null,
"discard_untrimmed": false,
"consider_already_trimmed": null
```

### 4. `pair_validation` placement is on R2 only — mismatch with per-file JSON

The plan says (line 132): "For R2 in paired-end mode, `pair_validation` is populated." This matches the current text report behavior (main.rs line 523-525: pair validation stats go in R2 report only).

However, this is a poor design for JSON. The pair validation stats belong to the **pair**, not to R2. A consumer parsing `sample_R1.fastq.gz_trimming_report.json` will see `"pair_validation": null` and have no way to know that pair validation happened — they must discover and open the R2 JSON file to find pair stats. This is a leaky abstraction from the text-report era.

**Recommendation:** Include `pair_validation` in **both** R1 and R2 JSON reports. The data is identical for both and the duplication is trivial. This makes each JSON file self-contained, which is the whole point of switching to structured output. Alternatively, add a `"pair_report"` field to R1 pointing to the R2 report filename, so consumers know where to look.

### 5. No `--no_report_file` interaction specified

The plan says (line 206-208): "The JSON report is always emitted (no opt-in flag needed)." But `--no_report_file` exists (cli.rs line 123-124) and currently suppresses the text report (main.rs lines 321, 491). The plan does not state whether `--no_report_file` should also suppress the JSON report.

If `--no_report_file` suppresses text but not JSON, that is a surprising asymmetry. If it suppresses both, that needs to be stated. Either way, this must be explicitly decided.

**Recommendation:** `--no_report_file` should suppress both text and JSON reports, since the flag's semantic is "do not write report files." State this explicitly in the plan.

---

## Recommendations

### 6. JSON string escaping must be specified

The plan mentions (line 193): "Helper: `json_string(w, key, value)`" but does not address JSON string escaping. The `command_line` field will contain user-supplied paths and arguments that may include:
- Backslashes (Windows paths, though unlikely for this tool)
- Double quotes (e.g., `--fastqc_args "--nogroup"`)
- Unicode characters in filenames

The `input_filename` and `adapter` fields could also contain characters needing escaping. The plan should specify that a `json_escape_string()` helper is needed that handles at minimum: `\`, `"`, `\n`, `\r`, `\t`, and control characters (`\u00xx`).

Without this, the tool will produce invalid JSON for certain inputs. This is a correctness issue.

### 7. `adapter_r2` in `parameters` but adapter_trimming only has one adapter

The JSON schema has `adapter_r2` in `parameters` (line 71) but only a single `adapter_trimming` section (lines 100-111). In paired-end mode, R1 and R2 may have different adapters (e.g., Small RNA, BGI). Each per-file JSON report should show the adapter that was actually used for **that** read's trimming, not the R1 adapter.

Looking at the code: in `run_paired` (main.rs line 499), the `adapter_seq` passed to both R1 and R2 `report_cfg` is the same string (the R1 adapter). The R2 adapter is `adapter_r2_seq`. But in `trimmer.rs` line 83-86, the actual adapter used for R2 is `config.adapter_r2.unwrap_or(&config.adapter)`.

**Recommendation:** The `adapter_trimming.sequence` field in R2's JSON should be the R2 adapter (or the R1 adapter if no R2-specific adapter was set). The plan should explicitly state how to resolve the adapter for each read.

### 8. `adapter_trimming.type` is hardcoded to "regular 3'"

The plan shows `"type": "regular 3'"` (line 103). This matches the current behavior (all adapters are searched as 3' adapters via `find_3prime_adapter` in trimmer.rs line 92). However, if future versions add 5' adapter support or linked adapters, this field becomes wrong. This is acceptable for now but should be noted as a future consideration.

### 9. Consider adding a `tool` or `generator` top-level field

The schema has `trim_galore_version` but no explicit tool identifier. If MultiQC discovers a `.json` file, it needs to know this is a TrimGalore report and not some other tool's JSON. Adding `"tool": "Trim Galore"` or `"report_type": "trim_galore_trimming_report"` would make file identification unambiguous.

### 10. `length_distribution` keys include non-zero entries only — document the invariant

The plan says (line 152): "sparse representation avoids emitting hundreds of zero entries." This is good, but the plan should explicitly state: **keys with count 0 are omitted**. This prevents an implementer from conditionally including some zeros.

Looking at the code in `write_cutadapt_section` (report.rs line 330-331): `if length == 0 || count == 0 { continue; }` — this skips both index 0 (unused) and zero-count entries. The JSON writer should follow the same logic.

### 11. JSON filename should respect `--output_dir` and `--basename`

The plan specifies the filename convention (lines 47-54) but does not discuss interaction with `--output_dir` or `--basename`. Looking at the text report naming:
- `report_name()` in io.rs (line 96-109) uses the input filename, not `--basename`, and respects `--output_dir`.
- The text report is always `{input_basename}_trimming_report.txt`.

The JSON report should parallel this exactly: `{input_basename}_trimming_report.json`, in the same directory as the text report. The plan should add a `json_report_name()` function to `io.rs` or explicitly state the naming follows `report_name()` with `.json` instead of `.txt`.

### 12. `basepair_processing` is missing `adapter_trimmed_bp`

The text report's Cutadapt section has "Quality-trimmed" bp, but it does not have a separate "adapter-trimmed bp" field. The JSON schema's `basepair_processing` has `quality_trimmed_bp` and `total_bp_written` but not `adapter_trimmed_bp`. This is consistent with the current text report, but it means the bp removed by adapter trimming can only be inferred (`total_bp_processed - quality_trimmed_bp - total_bp_written - bp_lost_to_other_trimming`). For a machine-readable format, it would be a good addition, but this is a nice-to-have rather than a requirement.

### 13. Error handling for JSON file creation

The plan does not discuss what happens if the JSON file cannot be created (disk full, permissions, etc.). The text report write currently uses `File::create(&report_path)?` with `?` propagation (main.rs line 346, 518), which will abort the entire run. The JSON write should follow the same pattern (fail-fast), but the plan should state this explicitly.

---

## Observations

### 14. "No serde" approach is reasonable

The schema is flat (max 2 levels), types are all primitives, and the only collection is `length_distribution` (a sparse map). Manual JSON writing with `write!` macros is straightforward and avoids adding ~15s to compile time (serde + serde_json). The codebase already has zero serde usage (confirmed in Cargo.toml). This decision is sound.

### 15. Schema versioning is well-designed

`schema_version: 1` with the rule "bump on field additions" is the right approach for a format consumed by an external tool (MultiQC). This allows MultiQC to handle schema evolution gracefully.

### 16. "Always emit all sections with zeros" is the right call

Conditional section omission is the #1 source of bugs in JSON consumers. Having `poly_a_trimming`, `poly_g_trimming`, and `rrbs` always present (even with zeros) means MultiQC can use a fixed parser without null checks. This is good.

### 17. The text report remains unchanged — no backwards compatibility risk

The plan correctly identifies that the JSON is purely additive: a new file alongside the existing text report. Existing MultiQC users will see no change. This is confirmed by the plan's "What Changes in the Text Report" section: "Nothing."

### 18. Test plan is thin

The test plan (lines 228-234) lists 4 categories but lacks specifics:
- No mention of testing with special characters in filenames
- No mention of testing `--no_report_file` interaction
- No mention of testing paired-end reports (R1 vs R2 field differences)
- No mention of testing `--basename` / `--output_dir` naming
- No mention of edge cases: zero reads, all reads filtered, adapter not found
- No mention of validating JSON structure with a JSON parser (not just "basic validator")

**Recommendation:** Add specific test cases:
1. SE run: verify JSON fields match text report values
2. PE run: verify R1 and R2 JSONs have correct `read_number`, adapter sequences
3. `--no_report_file`: verify no JSON is created
4. `--output_dir`: verify JSON lands in the correct directory
5. Input filename with spaces or special characters: verify valid JSON
6. Empty input file: verify JSON is valid with all-zero stats
7. Roundtrip: parse JSON with `serde_json` in test (dev-dependency only) to validate structure

### 19. `write_run_stats` is called nowhere in main.rs

The `write_run_stats()` function in report.rs (line 174) appears to be dead code — `main.rs` calls `write_report_header`, `write_cutadapt_section`, and `write_run_footer`, but never `write_run_stats`. This is not a problem for the JSON plan, but it is worth noting that the text report currently does NOT include the `=== Summary ===` section that `write_run_stats` would produce (it relies on the Cutadapt section's summary instead). The JSON should be comprehensive regardless.

---

## Summary of Required Changes

| # | Severity | Action |
|---|----------|--------|
| 3 | Critical | Add missing parameters: `clip_r1`, `clip_r2`, `three_prime_clip_r1`, `three_prime_clip_r2`, `max_n`, `discard_untrimmed`, `consider_already_trimmed` |
| 4 | Critical | Decide: include `pair_validation` in both R1 and R2 JSON, or document why R2-only is acceptable |
| 5 | Critical | Specify `--no_report_file` interaction with JSON report |
| 6 | High | Specify JSON string escaping requirements |
| 7 | High | Clarify which adapter sequence goes in R2's `adapter_trimming.sequence` |
| 11 | Medium | Document `--output_dir` and `--basename` interaction with JSON filename |
| 9 | Low | Consider adding a `tool` identifier field |
| 18 | Medium | Expand test plan with specific test cases |
