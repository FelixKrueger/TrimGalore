# Code Review — Clippy + fmt-drift pass on `src/`

Scope: single-reviewer review of substantive (non-whitespace) changes in the working-tree diff under `src/`. Pure cargo-fmt reformatting (line wraps, import reordering via `anyhow::{Result, bail}`, trailing-comma normalization, etc.) is explicitly excluded from this review.

Method: `git diff HEAD -- src/` inspected item-by-item; where behavior hinged on a subtle semantic, the surrounding context was read from the post-edit file.

---

## Verdicts per expected substantive change

### 1. `src/alignment.rs:~58` — `needless_range_loop` → `iter_mut().enumerate().skip(1)`
Verdict: benign (OK)
Rationale: `for (i, row) in dp.iter_mut().enumerate().skip(1) { row[0] = i; }` produces exactly the same write pattern as the original `for i in 1..=m { dp[i][0] = i; }`. `skip(1)` keeps the enumeration index aligned (starts at i=1), and borrowing is fine since the subsequent nested loop takes `dp` by non-overlapping index reads. No behavior change.

### 2. `src/quality.rs:~135` — `needless_range_loop` → `iter().enumerate().take(n)`
Verdict: benign (OK)
Rationale: `for (i, &base) in sequence.iter().enumerate().take(n)` is equivalent to `for i in 0..n { let base = sequence[i]; }` given `n <= sequence.len()` (the surrounding code establishes `n = sequence.len()`). The bound check is identical. `i` still advances over `0..n` and `best_index` logic is unaffected.

### 3. `src/trimmer.rs:~73` — `doc_lazy_continuation` indent fix
Verdict: benign (OK)
Rationale: Pure doc-comment reformatting. The "2.5. RRBS trim" line is now indented to render as a sub-bullet under step 2 rather than dropping rustdoc's lazy-continuation warning. No code impact.

### 4. `src/trimmer.rs:~170` — `collapsible_if`
Verdict: benign (OK)
Rationale: Nested `if !(config.is_paired && is_r2) { if record.seq.len() >= 2 && had_adapter { ... } }` collapsed to a single `&&`-joined guard. Short-circuit order preserved, guard composition is pure AND — semantically identical.

### 5. `src/trimmer.rs` — `#[allow(clippy::too_many_arguments)]` on `run_paired_end`
Verdict: note (acceptable, not ideal)
Rationale: `#[allow]` is the pragmatic choice here — `run_paired_end` is the top-level paired-end entry point and its argument list reflects genuine call-site needs (two readers, four writers, config, two unpaired-length thresholds, etc.). A config-struct refactor is a defensible future improvement but out of scope for a clippy pass. Flag for someday-cleanup.

### 6. `src/trimmer.rs:~647` — `assert!(!result.adapter_matches.is_empty())`
Verdict: benign (OK)
Rationale: Replacing `assert!(result.adapter_matches.len() >= 1)` with `assert!(!result.adapter_matches.is_empty())` is the idiomatic `len_zero` fix. Behavior identical; error message on failure slightly changes but test intent is preserved.

### 7. `src/parallel.rs` — 3× `#[allow(clippy::too_many_arguments)]`
Verdict: note (acceptable, not ideal)
Rationale: Applied to `run_paired_end_parallel`, `process_paired_batch`, and `process_pairs<W: Write>`. Same reasoning as item 5 — these are the parallel-pipeline plumbing functions where the argument count follows from worker-thread scaffolding, not accidental complexity. A `PairedIoBundle` struct could reduce this, but the `#[allow]` is the correct low-risk call for a clippy-only pass. Noted for future refactor.

### 8. `src/report.rs:~1087` — `range_start <= len - 1` → `range_start < len`
Verdict: benign (OK) — with caveat on the stated rationale
Rationale: The transformation is algebraically equivalent for `len >= 1`, and safer by construction (no `usize` subtraction). The enclosing loop is `for len in 1..=adapter_len`, so `len` is **always ≥ 1** and the original `len - 1` never actually underflowed in any reachable state. The plan's "guards underflow when len=0" framing is slightly overstated — the old code was *defensively* fragile but not buggy. The new form is still preferable: one less subtraction to reason about. Output formatting (the `len - 1` inside `format!`) is unchanged and still safe because that branch is only entered when `len > range_start >= 1`, i.e. `len >= 2`.

### 9. `src/demux.rs:~112` — `write!(w, "...\n")` → `writeln!(w, "...")`
Verdict: benign (OK)
Rationale: `writeln!` appends a single `\n` on all platforms (it is NOT `\r\n` on Windows — `writeln!` uses LF unconditionally), matching the original `\n` literal exactly. Byte-for-byte identical output. Idiomatic clippy fix for `write_with_newline`.

---

## Summary

All 9 substantive changes preserve behavior. Two (`#[allow(clippy::too_many_arguments)]` on `trimmer.rs` and `parallel.rs`) sidestep a real design smell rather than fix it — acceptable for a clippy-cleanup commit, but worth a follow-up ticket to bundle paired-IO arguments into a struct.

The `src/report.rs` change is correct and a small robustness win, though the "underflow on len=0" justification doesn't match the code (loop starts at 1). Worth a one-line commit-message tweak if this lands standalone.

No concerns blocking merge. No hidden behavior changes detected in the reviewed diff.
