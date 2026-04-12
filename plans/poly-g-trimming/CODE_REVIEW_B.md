# Code Review B: Poly-G Trimming Implementation

**Reviewer:** B  
**Date:** 2026-04-12  
**Scope:** 7 files under `optimus_prime/src/` — poly-G homopolymer trimming feature  
**Plan:** `plans/poly-g-trimming/plan.md`  
**Build status:** Compiles, 74 tests pass, clippy passes (pre-existing warnings only)

---

## Summary

The implementation is solid and complete. The core algorithm generalisation (`homopolymer_trim_index`), adapter-piggybacked detection, standalone fallback scan, CLI wiring, and reporting are all correctly implemented and follow the plan closely. The code is clean, well-documented, and consistent with the existing poly-A pattern. I found no critical bugs. There are two medium issues (plan amendments A2/A5 not fully implemented, off-by-one in scan count) and several low-priority observations.

---

## Issues by Area

### Logic

**M1. Plan amendments A2 and A5 not implemented: poly-G stats not shown when count is zero**  
**Priority: Medium**

Plan amendment A2 says: "Print stats when poly-G was enabled, even if count is 0 (confirms the feature ran)." Plan amendment A5 says: "When `config.poly_g` is true, always print the poly-G stats line in the report (even if 0 reads were trimmed)."

Neither is implemented:
- `report.rs` line 207: `write_run_stats` guards on `stats.poly_g_trimmed > 0`, and the function does not receive a config to check whether poly-G was enabled.
- `main.rs` line 294 (single-end) and lines 445/449 (paired-end): console summary also guards on `> 0`.

This means if poly-G is auto-detected or force-enabled on 4-colour data, the user sees no confirmation in the report or console that the feature actually ran. This is not a correctness bug (trimming still works), but it degrades observability.

**Recommendation:** Pass `config.poly_g` (or a bool) into `write_run_stats`, and change the guard from `if stats.poly_g_trimmed > 0` to `if config.poly_g`. Same for the console `eprintln!` blocks in `run_single` and `run_paired`.

---

**M2. Off-by-one in `reads_scanned` count inflates denominator by 1**  
**Priority: Medium**

Both `autodetect_adapter()` (line 96-100) and `detect_poly_g()` (line 236-240) use this pattern:

```rust
while let Some(record) = reader.next_record()? {
    reads_scanned += 1;
    if reads_scanned > MAX_SCAN_READS {
        break;
    }
    // ... process record ...
}
```

When the file has >= 1M reads, the last iteration increments `reads_scanned` to 1,000,001, checks `> MAX_SCAN_READS` (true), and breaks. But the record at count 1,000,001 was read from the file and counted in `reads_scanned` yet NOT processed (not checked for adapters or poly-G). So `reads_scanned` is 1,000,001 while only 1,000,000 reads were actually scanned.

This is a pre-existing bug (the adapter detection has the same pattern), but it now affects poly-G detection too:
- The threshold calculation `reads_scanned / 10_000` becomes 100 instead of 100 (no practical effect due to integer division).
- The percentage displayed (`poly_g_count / reads_scanned * 100`) is fractionally off: 100/1,000,001 vs 100/1,000,000.

**Recommendation:** Change to `if reads_scanned >= MAX_SCAN_READS { break; }` or move the break check to the beginning of the next iteration. This is a pre-existing issue so it could be fixed separately, but worth noting.

---

### Efficiency

**L1. `detect_poly_g` re-opens and decompresses the entire file**  
**Priority: Low**

When adapter auto-detection is skipped (user specified `--illumina`, `--adapter`, etc.), the standalone `detect_poly_g()` does a full 1M-read scan that re-opens and decompresses the input file. This is architecturally correct (the plan explicitly calls for this fallback), and the cost is acceptable since:
- It only scans 1M reads (a few seconds for gzipped files).
- It only runs when auto-detection was skipped AND neither `--poly_g` nor `--no_poly_g` was specified.

No action needed. Documenting for completeness.

---

### Errors / Edge Cases

**L2. No test for `detect_poly_g` standalone function**  
**Priority: Low**

The plan calls for integration tests (section 3: "Auto-detection test: Create a synthetic FASTQ where ~5% of reads end in 20 G's"). The unit tests for `has_trailing_poly_g` are thorough, but there's no test for the `detect_poly_g` function itself or its interaction with `main.rs` auto-detection logic. This is understandable (it requires FASTQ fixtures), but should be tracked.

---

**L3. Paired-end console summary missing poly-G bases removed count**  
**Priority: Low**

In `main.rs`, the single-end summary (line 296-298) prints both read count and bases removed for poly-G:
```
Reads with poly-G/C trimmed:     N (X.X%)
  Poly-G/C bases removed:        N
```

But the paired-end summary (lines 445-452) only prints read counts, not bases removed:
```
R1 reads with poly-G trimmed:    N (X.X%)
R2 reads with poly-C trimmed:    N (X.X%)
```

This is inconsistent (the poly-A paired-end summary has the same omission, so it's a pre-existing pattern, but still worth noting for the new feature).

**Recommendation:** Add `poly_g_bases_trimmed` lines under each R1/R2 poly-G line in `run_paired`, matching the single-end output format.

---

### Structure / Style

**L4. Clippy warning on `resolve_adapter` return type**  
**Priority: Low**

Clippy flags:
```
very complex type used. Consider factoring parts into `type` definitions
fn resolve_adapter(cli: &Cli) -> Result<(String, String, Option<String>, Option<(usize, usize)>)>
```

The poly-G feature added the `Option<(usize, usize)>` to an already-complex return type. A named struct would improve readability:
```rust
struct AdapterResolution {
    name: String,
    seq: String,
    seq_r2: Option<String>,
    poly_g_scan: Option<(usize, usize)>,
}
```

This is pre-existing complexity that the poly-G feature extended. Not a blocker.

---

**L5. Plan deviation: user-override messages don't show scan percentages**  
**Priority: Low**

The plan (section 4, main.rs) specifies that when `--poly_g` or `--no_poly_g` is used, the message should still show the auto-detection scan percentage so users can sanity-check. The implementation (lines 142-146) prints a simpler message without scan data, because no scan runs when CLI overrides are active.

This is a reasonable simplification (avoids a wasted 1M-read scan just for an informational message). Documenting as a conscious deviation rather than a bug.

---

**L6. `has_trailing_poly_g` could be `pub(crate)` for testing**  
**Priority: Low**

The function is currently `fn` (private). Unit tests are in the same module so this works, but if integration tests need to test detection logic directly, it would need to be made visible. Current visibility is fine for the existing test structure.

---

## Fixes Applied

None. No unambiguous low-risk fixes were identified that could be applied without discussion.

---

## Verification Checklist

| Plan Item | Status | Notes |
|-----------|--------|-------|
| `quality.rs`: Generalise `poly_a_trim_index` to `homopolymer_trim_index` | DONE | Clean refactor with complement mapping; wrapper preserved |
| `quality.rs`: Unit tests for poly-G/poly-C | DONE | 8 tests covering clear tail, no tail, short, exactly-3, errors, all-G, revcomp, post-poly-A |
| `adapter.rs`: `poly_g_count` in `DetectionResult` | DONE | Field added with doc comment |
| `adapter.rs`: Piggyback counting in `autodetect_adapter` | DONE | Counting added to existing scan loop |
| `adapter.rs`: `has_trailing_poly_g` helper | DONE | Private fn with 8 unit tests |
| `adapter.rs`: `detect_poly_g` standalone scan | DONE | Amendment A1 implemented correctly |
| `cli.rs`: `--poly_g` flag with aliases | DONE | `poly-g`, `polyG` aliases; conflicts with `no_poly_g` |
| `cli.rs`: `--no_poly_g` flag with aliases | DONE | `no-poly-g`, `no-polyG` aliases |
| `cli.rs`: A4 — independence note in help | DONE | "This is independent from --nextseq" in help text |
| `trimmer.rs`: `poly_g` in `TrimConfig` | DONE | Bool field alongside `poly_a` |
| `trimmer.rs`: `poly_g_trimmed` in `TrimResult` | DONE | usize field, 0 = none |
| `trimmer.rs`: Step 2.8 after poly-A | DONE | Correct position, correct revcomp logic |
| `trimmer.rs`: Stats collection in `run_single_end` | DONE | Matches poly-A pattern |
| `trimmer.rs`: Stats collection in `run_paired_end` | DONE | Both R1 and R2 stats |
| `report.rs`: `poly_g_trimmed`/`poly_g_bases_trimmed` in TrimStats | DONE | Fields + merge |
| `report.rs`: `poly_g` in report::TrimConfig | DONE | A3 naming correct |
| `report.rs`: Header line when poly_g enabled | DONE | Descriptive message |
| `report.rs`: Stats output | DONE | But only when count > 0 (see M1) |
| `main.rs`: Resolve adapter returns poly-G data | DONE | Tuple with Option for piggybacked data |
| `main.rs`: Standalone detect_poly_g fallback | DONE | Runs when adapter auto-detect skipped |
| `main.rs`: CLI overrides | DONE | `--poly_g` > `--no_poly_g` > auto-detect |
| `main.rs`: On-screen messages | DONE | 4 branches: force-on, force-off, auto-enabled, auto-disabled |
| `main.rs`: Console summary | DONE | But only when count > 0 (see M1) |
| `main.rs`: Wired to TrimConfig | DONE | `config.poly_g = poly_g_enabled` |
| `parallel.rs`: Poly-G stats in `process_pairs` | DONE | Both R1 and R2 |
| `parallel.rs`: Poly-G stats in `process_reads` | DONE | Single-end path |
| A5: Print poly-G stats even when 0 if enabled | NOT DONE | See M1 |
| A2: Console summary always when enabled | NOT DONE | See M1 |

---

## Priority Summary

| Priority | Count | Items |
|----------|-------|-------|
| Critical | 0 | -- |
| High | 0 | -- |
| Medium | 2 | M1 (stats not shown when zero), M2 (off-by-one scan count) |
| Low | 6 | L1-L6 (efficiency, missing integration test, paired summary inconsistency, clippy, plan deviation, visibility) |

**Verdict:** The implementation is correct and ready for use. The two medium issues are both non-blocking (M1 is observability, M2 is cosmetic precision). They can be addressed in a follow-up without risk to correctness.
