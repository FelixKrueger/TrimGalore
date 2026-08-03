# Plan: Multi-file SE input + --discard-untrimmed

## Context

User tested a drop-in command from Perl TrimGalore:
```
trim_galore --cutadapt_args '--discard-untrimmed' -a TCTCCTTGCATAATCACCAACC -j 8 *fastq.gz
```
Two issues surfaced:
1. `*fastq.gz` expands to 3 files — Rust version rejects >1 SE input
2. `--cutadapt_args '--discard-untrimmed'` is silently ignored (deprecated flag)

Both are drop-in compatibility gaps. The Perl version loops over all SE files and passes arbitrary args to Cutadapt.

---

## Feature 1: Multi-file single-end input

### Current behavior
- `cli.rs:268-274`: Rejects `input.len() > 1` when `!paired` (unless clock/implicon)
- `main.rs:182`: Calls `run_single()` once with `cli.input[0]`

### Perl behavior
- Loops over all positional input files (line ~581): `foreach my $filename (@ARGV) { trim($filename, ...) }`
- Each file gets its own adapter auto-detection, trimming, report, and optional demux

### Proposed change
- Remove the `input.len() > 1` guard for single-end mode
- In `main()`, wrap the single-end path in a loop over `cli.input`
- Per-file: adapter auto-detection, poly-G detection, trimming, report, demux
- Adapter auto-detection and poly-G scan must run per-file (different files may have different adapters)
- User-specified adapter (`-a`) skips auto-detection for all files (same as Perl)

### Files changed
- **`cli.rs`**: Remove the single-end >1 file validation (lines 268-274)
- **`main.rs`**: Wrap everything from `resolve_adapter()` through `run_single()` in `for input in &cli.input { ... }`; extract the per-file logic into a helper

### Edge cases
- `--basename` with multiple files: should error (ambiguous output naming). Add validation.
- `--demux` with multiple files: should work (each file gets its own demux)
- Mixed single/paired: not supported (already enforced by `--paired` requiring exactly 2)

---

## Feature 2: `--discard-untrimmed`

### Cutadapt semantics
`--discard-untrimmed`: Only keep reads where at least one adapter was found. Reads without adapter matches are discarded entirely (not written to output).

### Where it hooks in
The trimmer already tracks `had_adapter` per read (`TrimResult.had_adapter`). The filter step in `run_single_end()` (trimmer.rs:240-253) and `run_paired_end()` need a new filter:

```
if discard_untrimmed && !result.had_adapter {
    stats.discarded_untrimmed += 1;
    continue;  // skip writing this read
}
```

For paired-end: discard the pair if NEITHER read had an adapter (matches Cutadapt PE behavior).

### Files changed
- **`cli.rs`**: Add `--discard-untrimmed` flag (bool)
- **`trimmer.rs`**: Add `discard_untrimmed: bool` to `TrimConfig`. Add filter logic in `run_single_end()` and `run_paired_end()`. 
- **`report.rs`**: Add `discarded_untrimmed: usize` to `TrimStats`. Add line to report output.
- **`main.rs`**: Wire `cli.discard_untrimmed` → `config.discard_untrimmed`. Print summary stat.
- **`parallel.rs`**: Pass through `discard_untrimmed` in the parallel paths (uses same `TrimConfig`)

### Report line
Match TrimGalore/Cutadapt style:
```
Reads discarded as untrimmed:    N (X.X%)
```

### Interaction with `--cutadapt_args`
Update the deprecation warning for `--cutadapt_args` to mention that `--discard-untrimmed` is now a native flag.

---

## Implementation order

1. **`--discard-untrimmed`** — smaller, self-contained, more urgent (user's active workflow)
2. **Multi-file SE** — refactor main loop, slightly larger surface area

## Tests

- `cargo test` — all existing tests must pass
- Manual: `trim_galore --discard-untrimmed -a AGATCGGAAGAGC test_files/illumina_10K.fastq.gz` — verify discarded count > 0
- Manual: `trim_galore test_files/illumina_10K.fastq.gz test_files/illumina_10K.fastq.gz` — verify both processed (multi-file SE)
- CI validation unchanged (these are new features, not regressions)
