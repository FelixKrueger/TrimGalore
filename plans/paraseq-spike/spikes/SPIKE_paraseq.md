# SPIKE_paraseq — Does paraseq net out faster than our current FastqReader?

## 1 · Question, success criteria, strategy

**Question.** Does swapping `paraseq` in for our line-by-line `String`-allocating FASTQ parser net out faster, **after** the moment we have to materialize owned records to ship them across the worker-pool channel?

This is the load-bearing question from issue #315 — paraseq's zero-copy is obviously fast in isolation; the real test is whether the saved allocations survive integration into our existing reader-thread → worker-pool architecture, which requires owned records at the channel boundary.

**Success criteria.**
- ≥5% wall-clock improvement on read-parsing throughput (plain FASTQ, 500K × 150 bp ≈ 157 MB) when measured after owned-record conversion.
- All variants produce identical sentinel checksums (read count, sum of seq lengths, sum of last seq byte) to prove no parsing was elided.

**Scope boundary.**
- **In:** parsing throughput on plain (uncompressed) FASTQ. Three variants — baseline, paraseq zero-copy, paraseq + owned conversion.
- **Out:** gzip decompression overhead (would only widen the absolute numbers without changing the recommendation; see §7); paired-end + `--passthrough` (extra cross-stream sync logic); integration into `parallel.rs`; actual trim+compress wall-clock impact.

**Strategy.** Standalone Cargo project under `plans/paraseq-spike/spikes/`. Three (then four) loop bodies over the same 157 MB fixture, each computing the same sentinel value (sum of seq lengths + sum of last-seq-byte) to prevent dead-code elimination. Times reported as median / min / max of 3 runs each.

## 2 · Script path and how to run

```
plans/paraseq-spike/spikes/Cargo.toml
plans/paraseq-spike/spikes/spike_paraseq.rs
```

Build:
```
cd plans/paraseq-spike/spikes
cargo build --release
```

Generate the 500K-read synthetic fixture (deterministic LFSR-driven bases, 150 bp, ~157 MB plain) and benchmark:
```
./target/release/spike_paraseq build-fixture 500000 fixture.fastq
./target/release/spike_paraseq bench fixture.fastq
```

## 3 · Iterations

### Iteration #1 — initial three variants

Initial cut hit two compile errors:
- Trait import path: `paraseq::Record`, not `paraseq::fastx::Record`.
- Method name: `record.qual()` (returns `Option<&[u8]>`), there is no `qual_raw()`.

Fixed and ran. Result (3 runs, median):

| Variant | Median (ms) | vs baseline |
|---|---:|---:|
| (A) baseline (mirrors `src/fastq.rs:322-355`) | 85.5 | — |
| (B) paraseq zero-copy | 32.7 | **2.6× faster** |
| (C) paraseq → owned (naïve: `format!` + `from_utf8.to_string()`) | 101.8 | **19% slower** |

Headline finding: zero-copy is a big win; naïve owned-conversion eats the win **and then some**. Most of the regression in (C) traced to `format!("@{}", ...)` going through `Display`/`Write` machinery for every record.

### Iteration #2 — lean owned conversion

Added variant (D) that prepends `@` manually and uses `String::from_utf8(bytes.to_vec())` — one memcpy + one UTF-8 validation per field, no format-machinery overhead. Same fixture, same 3-run methodology.

| Variant | Median (ms) | vs baseline |
|---|---:|---:|
| (A) baseline | 82.6 | — |
| (B) paraseq zero-copy | 31.9 | **2.59× faster** |
| (C) paraseq → owned naïve | 99.5 | **+20%** |
| (D) paraseq → owned lean | 86.2 | **+4%** (within noise) |

Outcome: confirmed. (D) recovers the `format!` overhead (99.5 → 86.2) but lands at baseline. Three `String` allocations per record is fundamentally the same cost whether they come from `BufRead::read_line` + `trim_end_matches(...).to_string()` or `String::from_utf8(bytes.to_vec())`. The win in (B) lives entirely in the zero-allocation iteration.

### Iteration #3

Not run — the question is answered. The architectural recommendation does not change with additional measurements.

## 4 · Findings summary

1. **paraseq parses ~2.6× faster than our current line-by-line `String`-allocating reader.** That is a real, repeatable, large effect when records can stay borrowed.
2. **Any sensible owned-record conversion erases the win.** Even the lean variant (one memcpy + one UTF-8 validation per field, no `format!`) lands at baseline ±4%. The cost dominator is the *materialization*, not the *parsing*.
3. **Therefore: drop-in paraseq swap in our current architecture is value-negative or neutral.** Our reader thread → bounded channel → worker pool requires owned `FastqRecord` at the channel boundary. The zero-copy benefit cannot survive that hand-off.
4. **The win is only unlocked by an architectural change** — moving parsing *into* the worker threads so records stay borrowed end-to-end (worker pulls a raw byte block from a chunker thread, parses locally, then operates on borrowed slices through trim → filter → compress). That is a substantial rewrite of `parallel.rs`, not a parser swap.
5. **Even with the architectural change, the relative whole-job win is small.** At Buckberry scale the dominant cost is output gzip compression (the reason `DEFAULT_GZIP_LEVEL` is 1 and not 6). 50 ms saved per 500K reads on parsing is a low-single-digit percentage of total wall-clock once decompression + trim + compress are stacked on top.

## 5 · Reference snippets

The paraseq API shape an implementation agent would otherwise re-derive:

```rust
use paraseq::Record;        // trait providing seq_raw / qual / id_str
use paraseq::fastq;

let mut reader = fastq::Reader::from_path(path)?;
let mut record_set = reader.new_record_set();

while record_set.fill(&mut reader)? {
    for record in record_set.iter() {
        let record = record?;
        let id  = record.id_str();          // &str (no '@' prefix)
        let seq = record.seq_raw();         // &[u8]
        let qual = record.qual();           // Option<&[u8]>
        // process — but DO NOT clone into String unless absolutely necessary
    }
}
```

Two non-obvious points worth carrying forward:

- `id_str()` does **not** include the leading `@` — if you want byte-identity with our current `FastqRecord::id` you must prepend it manually.
- `qual()` returns `Option<&[u8]>` (paraseq is shared with FASTA which has no qual line). For FASTQ-only paths it's always `Some`.

## 6 · Recommendation

**Do not do the drop-in paraseq swap.** It is a wash or a regression in our current architecture (one reader thread → owned-record channel → worker pool).

**Park the broader idea** of restructuring workers to parse-in-place from raw byte chunks. The architectural cost is substantial:
- The reader thread becomes a chunker (carve raw bytes on record boundaries, no parsing).
- Each worker owns a chunk, parses with paraseq, then trims/filters/compresses on borrowed slices.
- Three-way header-sync for `--passthrough` becomes harder because per-record sync no longer happens in one place.
- The clumpy minimizer dispatcher (`read_pairs_clumpy`) which routes by canonical minimizer key needs rethinking — the dispatcher *itself* currently parses to compute the key.

…and the upside is modest (single-digit % wall-clock at typical workloads, because compression dominates).

**Suggested resolution for issue #315:** post a "spike result" comment summarising the numbers, close the paraseq sub-thread as "interesting at the parser level, value-negative at the system level given our architecture and gzip dominance." Keep `noodles-bam` / uBAM (separate concern in #316) on the active list.

**When to revisit:** if we ever ship an uncompressed-output mode (or a `--compression 0` profile where output gzip is removed), parsing moves to a larger fraction of total wall-clock and the architectural change becomes more attractive.

## 7 · Limitations

- **Plain FASTQ only.** Production inputs are gzipped (`.fastq.gz`). Decompression typically takes ~3–4× the time of parsing on the same content size, so absolute wall-clocks shift up significantly — but the *relative* gap between (A), (B), and (D) measured here is at the parser layer and would be additive on top of decompression, not erased by it. paraseq supports gzip via `niffler`; a follow-up measurement could quantify the gzip-input case.
- **Single-end only.** Paired-end with `--passthrough` adds per-record three-way `read_id_prefix` sync that would have its own cost in any restructured architecture.
- **Synthetic fixture.** 500K reads of 150 bp with deterministic LFSR-driven bases — read-length-uniform, ID-format-uniform. Real Illumina data has variable-length IDs and occasionally variable seq lengths; could shift the constants but not the qualitative finding (zero-copy win lives in iteration; conversion cost lives in `String` materialization).
- **Macro-level not measured.** This spike only measured the parser. The full TrimGalore wall-clock is parser + trim + compress (and decompress on input). Trim and compress are unchanged across variants and therefore not part of this question.
- **Hardware.** Single run on the author's Darwin/arm64 host. Numbers are relative-comparison-grade, not absolute-throughput-grade.

## 8 · Disposition

- **Outcome:** spike answers the load-bearing question with high confidence. No follow-up spike needed at the parser layer.
- **Action:** see §6. Recommend posting the numbers as a comment on #315 and closing the paraseq sub-thread for now.
- **Throwaway artefacts:** `plans/paraseq-spike/spikes/` (Cargo project, script, fixture). Safe to delete after the recommendation is posted.
