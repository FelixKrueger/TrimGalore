# Plan: JSON Trimming Report for MultiQC

## Goal

Emit a structured JSON report alongside the existing text report, as suggested
by @ewels in MultiQC/MultiQC#3529. This gives MultiQC a clean, regex-free
format to parse, while the text report stays unchanged for backwards compat.

---

## Context

Phil's proposal (MultiQC/MultiQC#3529):
> How about we add a new file format produced by v2 which is easier to parse
> (eg. JSON) and has everything. The old backwards compatible log file just
> needs something new in the first few lines to tell us to skip those files.

**Current state:** MultiQC parses TrimGalore reports via its Cutadapt module,
requiring a "This is cutadapt" shim, regex-based extraction, and various
version-sniffing heuristics. This is fragile and shows "Cutadapt 4.0" in the
Software Versions panel.

**Target state:** MultiQC gets a native TrimGalore module that reads JSON.
When JSON is found, it skips the text report. Old text reports (v0.6.x and
v2.0 with the Cutadapt shim) still work via the existing Cutadapt parser.

---

## Deliverables

| # | Where | File | Purpose |
|---|-------|------|---------|
| 1 | TrimGalore | `src/report.rs` | New `write_json_report()` function + `json_escape_string()` helper |
| 2 | TrimGalore | `src/io.rs` | New `json_report_name()` function |
| 3 | TrimGalore | `src/main.rs` | Call `write_json_report()` after text report |
| 4 | MultiQC (upstream) | New TrimGalore module | Parse the JSON (PR to MultiQC, coordinate with Phil) |

This plan covers #1, #2, and #3 only. #4 is a separate MultiQC PR once we
agree on the schema with Phil.

---

## JSON Schema

### Filename Convention

```
{input_basename}_trimming_report.json
```

Examples:
- SE: `sample.fastq.gz_trimming_report.json`
- PE: `sample_R1.fastq.gz_trimming_report.json`, `sample_R2.fastq.gz_trimming_report.json`

Parallel to the existing `*_trimming_report.txt` naming.

**`--output_dir`**: JSON report goes in the same directory as the text report
(respects `--output_dir` via a new `json_report_name()` in `src/io.rs` that
mirrors `report_name()` but with `.json` extension).

**`--basename`**: Does NOT affect report filenames. Report filenames always use
the original input filename, matching current text report behavior (confirmed:
`report_name()` in `io.rs:96` uses the input path, not `--basename`).

**`--no_report_file`**: Suppresses both text AND JSON reports. The JSON write
call goes inside the existing `if !cli.no_report_file` guard in `main.rs`
(lines 321 and 491). A user who says "no reports" expects no report files of
any kind.

### Schema (v1)

```json
{
  "tool": "Trim Galore",
  "schema_version": 1,
  "trim_galore_version": "2.0.0",
  "input_filename": "sample_R1.fastq.gz",
  "mode": "paired-end",
  "read_number": 1,
  "command_line": "trim_galore --paired --cores 8 sample_R1.fastq.gz sample_R2.fastq.gz",

  "parameters": {
    "quality_cutoff": 20,
    "adapter": "AGATCGGAAGAGC",
    "adapter_r2": "AGATCGGAAGAGC",
    "error_rate": 0.1,
    "stringency": 1,
    "length_cutoff": 20,
    "max_length": null,
    "phred_encoding": 33,
    "trim_n": false,
    "nextseq": false,
    "rrbs": false,
    "non_directional": false,
    "poly_a": false,
    "poly_g": false,
    "clip_r1": null,
    "clip_r2": 2,
    "three_prime_clip_r1": null,
    "three_prime_clip_r2": null,
    "max_n": null,
    "discard_untrimmed": false,
    "consider_already_trimmed": null
  },

  "read_processing": {
    "total_reads": 10000,
    "reads_with_adapter": 5234,
    "reads_written": 9800,
    "reads_too_short": 150,
    "reads_too_long": 0,
    "reads_too_many_n": 50,
    "reads_discarded_untrimmed": 0
  },

  "basepair_processing": {
    "total_bp_processed": 1000000,
    "quality_trimmed_bp": 50000,
    "total_bp_written": 900000
  },

  "adapter_trimming": {
    "sequence": "AGATCGGAAGAGC",
    "type": "regular 3'",
    "length": 13,
    "times_trimmed": 5234,
    "length_distribution": {
      "1": 1000,
      "2": 500,
      "3": 250,
      "5": 120
    }
  },

  "poly_a_trimming": {
    "reads_trimmed": 0,
    "bases_removed": 0
  },

  "poly_g_trimming": {
    "reads_trimmed": 0,
    "bases_removed": 0
  },

  "rrbs": {
    "trimmed_3prime": 0,
    "trimmed_5prime": 0
  },

  "pair_validation": {
    "pairs_analyzed": 10000,
    "pairs_removed": 200,
    "pairs_removed_n": 50,
    "pairs_removed_too_long": 0,
    "r1_unpaired": 10,
    "r2_unpaired": 15
  }
}
```

### Field notes

**`tool: "Trim Galore"`** — unambiguous file identification for MultiQC. When
MultiQC discovers a `*_trimming_report.json`, this field confirms it's a Trim
Galore report (not some other tool's JSON).

**`schema_version: 1`** — allows MultiQC to handle future schema changes
without breaking. If we add fields, we bump the version.

**`mode`** — derived from `report::TrimConfig.paired`: `"paired-end"` or
`"single-end"`. Not a field on the config struct itself.

**`read_number: 1|2`** — tells MultiQC which read this report belongs to
in paired-end mode. For SE, this is `1`. Caller passes `1` in
`run_single_file()`, `1` / `2` in `run_paired()`.

**`adapter_r2`** — `null` when no R2 adapter is configured (SE mode, or PE
without an explicit R2 adapter). Matches the `Option<String>` type on
`report::TrimConfig.adapter_r2`.

**`adapter_trimming.sequence`** — shows the adapter actually used for **this
read's** trimming. For R1/SE: uses the primary adapter. For R2: uses the R2
adapter if one was set (Small RNA, BGI presets), otherwise falls back to the
R1 adapter. This matches `trimmer.rs:83-86` logic.

**`adapter_trimming.type`** — hardcoded to `"regular 3'"` since all current
adapter matching is 3' anchored (`find_3prime_adapter`). If future versions
add 5' or linked adapters, this field will reflect the actual type and
`schema_version` will bump.

**`length_distribution` as object, not array** — sparse representation
(`{"5": 120, "13": 45}`) avoids emitting hundreds of zero entries. Keys
are string-encoded integers (JSON keys must be strings). Index 0 and entries
with count 0 are omitted (matching `write_cutadapt_section` logic at
`report.rs:330-331`: `if length == 0 || count == 0 { continue; }`).

**`pair_validation`** — included in BOTH R1 and R2 JSON reports in paired-end
mode. Unlike the text report (which puts pair validation in R2 only for
historical reasons), the JSON report is self-contained per file. Since v2.0
processes both reads simultaneously, both reports have access to the pair
stats. For SE mode, this field is `null`.

**`parameters` section** — comprehensive: includes all user-configurable
parameters that affect trimming behavior. `clip_r1`, `clip_r2`,
`three_prime_clip_r1`, `three_prime_clip_r2` use `null` when not set
(matching `Option<usize>` types). `max_n` is `null` when not set. This
ensures a consumer can fully reconstruct what processing was applied.

**`consider_already_trimmed`** — `null` when not used, integer threshold
when set. Allows a consumer to distinguish "adapter trimming was suppressed
because reads were already trimmed" from "no adapter was detected."

**Sections with zero values are always present** — `poly_a_trimming`,
`poly_g_trimming`, `rrbs` are always included (with zeros) rather than
conditionally omitted. This makes the schema predictable for consumers.

**No `expect` column in length_distribution** — the expected count
(`total_reads / 4^length`) is trivially computed from `total_reads` and the
length key. No need to store derived data. MultiQC can compute it if needed.

---

## Implementation Details

### `src/report.rs`

New function:

```rust
pub fn write_json_report<W: Write>(
    w: &mut W,
    config: &TrimConfig,              // report::TrimConfig (same struct used by text report)
    stats: &TrimStats,
    pair_stats: Option<&PairValidationStats>,
    read_number: u8,
    trim_config: &crate::trimmer::TrimConfig,  // for clip/max_n/discard params
) -> std::io::Result<()>
```

The function takes `report::TrimConfig` (the report-layer struct, same one
already constructed in `main.rs:324-344` and `main.rs:496-516` for the text
report). It additionally takes `trimmer::TrimConfig` for parameters not on
the report config (clip offsets, max_n, discard_untrimmed, consider_already_trimmed).

**No serde dependency.** The JSON is simple enough to write manually with
`write!`/`writeln!` macros. The schema is flat (max 2 levels deep) and the
types are all primitives (strings, integers, floats, bools, null). Adding
serde + serde_json for this would increase compile time and binary size for
no real benefit.

### JSON writing approach

**Helper functions** (matching the existing text report's `write!` pattern):

```rust
/// Escape a string for JSON output.
/// Handles: \ -> \\, " -> \", newline -> \n, tab -> \t,
/// carriage return -> \r, and other control characters -> \u00xx.
fn json_escape_string(s: &str) -> String

/// Write a JSON key-value pair with string value.
fn json_str<W: Write>(w: &mut W, key: &str, value: &str, comma: bool) -> io::Result<()>

/// Write a JSON key-value pair with integer value.
fn json_int<W: Write>(w: &mut W, key: &str, value: usize, comma: bool) -> io::Result<()>

/// Write a JSON key-value pair with float value.
fn json_float<W: Write>(w: &mut W, key: &str, value: f64, comma: bool) -> io::Result<()>

/// Write a JSON key-value pair with bool value.
fn json_bool<W: Write>(w: &mut W, key: &str, value: bool, comma: bool) -> io::Result<()>

/// Write a JSON key with null value.
fn json_null<W: Write>(w: &mut W, key: &str, comma: bool) -> io::Result<()>
```

**Trailing comma strategy:** Each helper takes a `comma: bool` parameter.
When `true`, a trailing comma is appended after the value. The last field in
each object passes `comma: false`. This is explicit and compile-time
checkable — no runtime state tracking needed.

Pretty-printed with 2-space indentation for human readability.

### `src/io.rs`

New function parallel to `report_name()`:

```rust
/// Generate the JSON trimming report filename.
pub fn json_report_name(input: &Path, output_dir: Option<&Path>) -> PathBuf {
    let input_name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    let report = format!("{}_trimming_report.json", input_name);
    match output_dir {
        Some(dir) => dir.join(&report),
        None => input.parent().unwrap_or(Path::new(".")).join(&report),
    }
}
```

### `src/main.rs`

**Single-end** (`run_single_file()`): Inside the existing `if !cli.no_report_file`
block (line 321), after the text report write:

```rust
// JSON report
let json_path = naming::json_report_name(input, output_dir);
let json_file = File::create(&json_path)?;
let mut jw = BufWriter::new(json_file);
report::write_json_report(
    &mut jw, &report_cfg, &stats,
    None,  // no pair_validation for SE
    1,     // read_number
    &config,
)?;
eprintln!("JSON report: {}", json_path.display());
```

**Paired-end** (`run_paired()`): Inside the existing `if !cli.no_report_file`
block (line 491), within the per-read loop:

```rust
// JSON report
let json_path = naming::json_report_name(input, output_dir);
let json_file = File::create(&json_path)?;
let mut jw = BufWriter::new(json_file);
report::write_json_report(
    &mut jw, &report_cfg, stats,
    Some(&pair_stats),  // BOTH R1 and R2 get pair_validation
    (idx + 1) as u8,    // read_number: 1 or 2
    &config,
)?;
eprintln!("JSON report: {}", json_path.display());
```

Note: in paired-end mode, `pair_validation` is passed to BOTH R1 and R2
reports (unlike the text report which only writes it for R2). The
`adapter_trimming.sequence` for R2 uses the actual R2 adapter from
`report_cfg.adapter_r2` (when set) or falls back to `report_cfg.adapter`.

### Error handling

JSON file creation uses `File::create()` with `?` propagation, matching the
text report pattern (fail-fast on disk errors). If the JSON file cannot be
created, the run aborts — same as for text reports.

### Output file creation

The JSON report is emitted whenever the text report is emitted (i.e., when
`--no_report_file` is NOT set). It's small (1-5 KB) and the information is
already computed. If users don't want it, they can use `--no_report_file` or
ignore/delete it. This is simpler than adding a `--json` flag.

---

## What Changes in the Text Report

**Nothing.** The text report stays exactly as-is, including the Cutadapt
shim section. The JSON report is a new additional file. MultiQC will prefer
the JSON when its TrimGalore module is added; until then, the Cutadapt shim
continues to work.

Phil suggested adding "something new in the first few lines to tell us to
skip those files" — but this is actually better handled on the MultiQC side:
if MultiQC finds a `*_trimming_report.json` next to a `*_trimming_report.txt`,
it uses the JSON and skips the text. No changes needed to the text report.

---

## Test Plan

### Unit tests (in `src/report.rs`)

1. **`test_json_escape_string`** — verify escaping of `\`, `"`, newline, tab,
   control characters, and a clean passthrough for normal strings and DNA
   sequences.

2. **`test_write_json_report_se`** — SE mode: construct a `TrimConfig` +
   `TrimStats` with known values, call `write_json_report()` into a `Vec<u8>`,
   parse the output with `serde_json::Value` (dev-dependency only) to verify:
   - Valid JSON structure
   - All field names present
   - Values match input stats
   - `pair_validation` is `null`
   - `read_number` is `1`
   - `mode` is `"single-end"`

3. **`test_write_json_report_pe`** — PE mode: pass `pair_stats`, verify
   `pair_validation` is populated, `mode` is `"paired-end"`, `read_number`
   matches.

4. **`test_write_json_report_sparse_length_distribution`** — verify that
   index 0 and zero-count entries are omitted from `length_distribution`.

5. **`test_write_json_report_special_characters`** — use a `command_line`
   containing quotes, backslashes, and an `input_filename` with spaces.
   Verify the output is valid JSON.

6. **`test_write_json_report_all_zero_stats`** — empty input edge case.
   Verify valid JSON with all-zero values.

### Integration tests

7. **SE run**: `cargo test` integration test using `test_files/illumina_10K.fastq.gz`.
   Verify both `.txt` and `.json` files are produced. Parse JSON with
   `serde_json` and validate field values match the text report.

8. **PE run**: Same with `test_files/` paired files. Verify both R1 and R2
   JSONs exist, both have `pair_validation` populated, `read_number` is
   correct.

9. **`--no_report_file`**: Verify no `.json` file is created.

10. **`--output_dir`**: Verify JSON lands in the specified directory.

### Dev dependency

Add `serde_json` as a **dev-dependency only** for test JSON parsing:
```toml
[dev-dependencies]
serde_json = "1"
```
This does not affect the release binary or compile time for non-test builds.

---

## Follow-up (out of scope for this plan)

- PR to MultiQC adding a native TrimGalore module that reads the JSON
- Coordinate with @ewels on any schema adjustments
- Post the schema to MultiQC/MultiQC#3529 for review before implementing
  the MultiQC module
