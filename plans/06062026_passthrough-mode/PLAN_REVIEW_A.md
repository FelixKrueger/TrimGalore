# Plan Review A — `--passthrough` mode

**Reviewer:** A
**Plan:** `plans/06062026_passthrough-mode/PLAN.md`
**Date:** 2026-06-06
**Verdict:** Sound overall direction; several non-trivial correctness and ergonomic gaps to close before implementation.

---

## 1. Logic review

### 1.1 Critical: `PairValidationStats` is not `PartialEq` — parity test won't compile as specified

The plan's load-bearing test (Validation §2) says:

> Assert `serial_stats == parallel_stats` (the `PartialEq` derive covers every field)

That works for `TrimStats` (`src/report.rs:8` derives `PartialEq, Default, Debug`) but **not** for `PairValidationStats` — at `src/report.rs:106` the derive is only `#[derive(Debug, Default)]`. The plan's Validation §2 will not compile until either:

- The implementer adds `PartialEq` (and probably `Eq`) to `PairValidationStats`. That's a one-line change but should be in the plan as a Step-3 sub-item; otherwise the parity assertion silently fails to typecheck.
- Or the test asserts field-by-field, which is fragile when new fields land later.

Add it to the plan explicitly. Same applies if the implementer wants to compare the existing `PairValidationStats` against itself (e.g. plain vs clumpy in another test).

### 1.2 Critical: sync check in the parallel reader thread cannot drive the `Discarded` accounting

The plan (Step 6.3, Step 7.3) says the three-way sync check fires in the **reader thread**, then the worker uses `classify_paired` and emits `PairOutcome::Pass / Discarded`. The reader thread does the work of pairing-up `(r1, r2, passthrough)` tuples and shipping them to a worker. That much is fine.

But: the worker thread (`process_pairs`) is the only place where the `Pass/Discarded` decision is made — and the worker needs the **same-row passthrough record** to (a) write it to the gzip member on `Pass`, (b) drop it on `Discarded`, (c) update `passthrough_records_kept/dropped` on `pair_stats`. The plan's `PairedWork` extension (`Option<Vec<FastqRecord>>`) carries the records correctly, and `process_pairs` would have to iterate `r1, r2, pt` in lockstep:

```rust
for ((r1, r2), pt) in reads_r1.iter_mut().zip(reads_r2.iter_mut()).zip(passthrough.iter()) {
    match classify_paired(r1, r2, ...) {
        Pass      => { r1.write_to(...); r2.write_to(...); pt.write_to(writer_pt); pair_stats.passthrough_records_kept += 1; }
        Discarded => { pair_stats.passthrough_records_dropped += 1; }
        Unpaired  => unreachable!()  // validation rejects retain_unpaired
    }
}
```

This is fine — but the plan needs to flag that `process_pairs` currently consumes `&mut [FastqRecord]` for r1/r2 (mutating them in-place during trim) and the passthrough must be `&[FastqRecord]` (no mutation — passthrough is written verbatim). The interleaved iteration above and the lifetime/borrow shape of the worker should be spelled out, because the current `process_pairs` signature (`parallel.rs:421-433`) doesn't have a slot for an immutable third slice.

### 1.3 Critical: `compressed_passthrough` always present, even when feature off — wastes a `Vec<u8>` per batch

The plan's `PairedBatchResult` (signature section) adds `compressed_passthrough: Vec<u8>` (NEW (empty when passthrough disabled)). Allocating an empty `Vec` per batch is fine, but consider `Option<Vec<u8>>` to match `compressed_unpaired_r{1,2}` semantics — those fields are populated only when `retain_unpaired` is on, but the plan now mixes idioms (Option-vs-empty-Vec). Pick one and make it consistent; "empty Vec means disabled" is footgunny when somebody later asks "did we have passthrough on this run?" and gets a confusing answer. I'd recommend `Option<Vec<u8>>` matching the upstream pattern.

### 1.4 Important: the `(None, None)` EOF arm of `read_pairs_round_robin` must drain passthrough — Step 5.5 mentions this for serial but not parallel

Plan Step 5.5 (serial): "The `(None, None)` arm: also consume from the passthrough reader and bail if it returns `Some`." Plan Step 6 (parallel) describes the sync check but doesn't explicitly call out the EOF arm — the implementer needs to handle it in `read_pairs_round_robin` too, otherwise a passthrough file with **more** records than R1/R2 is undetected (it just terminates when R1/R2 EOF without checking that the third file is also at EOF).

Add an explicit Step 6.4.1: "In `read_pairs_round_robin`'s `(None, None)` arm, also call `reader_passthrough.next_record()?`. If it returns `Some`, bail with `'passthrough has more records than R1/R2'`."

### 1.5 Important: header sync check needs to handle missing `@` prefix correctly

The plan's `read_id_prefix()`:

```rust
let trimmed = id.strip_prefix('@').unwrap_or(id);
trimmed.split_ascii_whitespace().next().unwrap_or(trimmed)
```

is fine for well-formed FASTQ. But Trim Galore's `FastqRecord::id` (`src/fastq.rs:32`) stores the **full** header line including the leading `@` (see `record.id.starts_with('@')` check at `src/fastq.rs:418`). The unit test `read_id_prefix("@read1") == "read1"` works, but worth asserting in the plan: "the input to `read_id_prefix` is the raw `id` field from `FastqRecord`, which still contains the leading `@`." Currently the plan documents the function but doesn't make this contract explicit.

Also: `split_ascii_whitespace().next().unwrap_or(trimmed)` is unreachable — `split_ascii_whitespace()` on a non-empty string always yields at least one segment unless the string is entirely whitespace. The function should probably just `unwrap_or("")` (an all-whitespace ID isn't a valid FASTQ header anyway, and `unwrap_or(trimmed)` would return the entire whitespace which defeats the purpose). Minor, but worth a clean spec.

### 1.6 Important: passthrough records must be **written before** the worker can recycle the slot

Plan Step 7.4 says "In the main thread's flush loop... when `output_passthrough.is_some()` write `r.compressed_passthrough` to the passthrough output file." This is correct. But: the parallel path already writes R1/R2 in lockstep (i.e., both written before incrementing `expected`). The plan must ensure the passthrough write happens **in the same loop iteration** as R1/R2, not as a separate sweep — otherwise a panic mid-flush leaves the three files at different progress points and `expected` could advance past unwritten passthrough data. The plan's natural reading does this correctly, but it should be made explicit: "Within the `while let Some(r) = pending.remove(&expected)` block, write R1, R2, *and* passthrough atomically before `expected += 1`."

### 1.7 Important: `--cores >= 2` + plain output (`--dont_gzip`) path

The plan's signature for `process_pairs<W: Write>` should not assume gzip — the existing function accepts `Vec<u8>` directly when `gzip == false` (see `parallel.rs:255-277`). The passthrough must follow the same shape: when `gzip == false`, the third writer is a `&mut Vec<u8>`; when `gzip == true`, it's a `&mut GzEncoder<&mut Vec<u8>>`. The plan's "scoped `GzEncoder` block" language (Step 7.1) is right but the implementer must mirror the existing `if gzip { ... } else { ... }` fork. Easy to miss because there are two equivalent calls to `process_pairs` in the function body. Flag it explicitly in Step 7.

### 1.8 Important: `--basename foo` plus a single-pair `--paired --passthrough` — what gets the basename?

Plan §3 (Behavior, step 3) says the basename applies uniformly: `foo_val_1`, `foo_val_2`, `foo_passthrough`. Assumption #2 in the plan is OK, but `Cli::validate()` already has at `cli.rs:438`:

```rust
if self.paired && self.input.len() > 2 && self.basename.is_some() { bail!(...) }
```

Adding `--passthrough` to a `--paired` invocation does not change `self.input.len()` (passthrough is a separate `Option<PathBuf>` field), so `--basename foo --paired --passthrough` works fine with the existing check. **But:** the plan rejects multi-pair input with passthrough (validation 1.ii) so `--basename` will never combine with multi-pair passthrough. The "what gets the basename" decision is therefore single-pair only, which is fine.

### 1.9 Important: case-folded preflight collides if user has e.g. `R1.fq.gz` and `r1.fq.gz`

The existing `out_paths` HashMap in `main.rs:178-213` is the right place to push the passthrough output (plan Step 2.3). One subtlety: the plan says the passthrough output uses the **passthrough input's** basename (Q&A #3). On a case-insensitive filesystem, `I1.fq.gz` and `i1.fq.gz` would both strip to `I1_passthrough.fq.gz` / `i1_passthrough.fq.gz` and case-fold to the same key — which is the **desired** behaviour (loud failure). But: if the user passes `R1.fq.gz R2.fq.gz --passthrough I1.fq.gz` from the same directory, the candidates are `R1_val_1.fq.gz`, `R2_val_2.fq.gz`, `I1_passthrough.fq.gz` — all distinct, fine. The case that **could** collide is `R1.fq.gz R2.fq.gz --passthrough R1.fq.gz` — but plan validation step 1.ix already catches "passthrough cannot be the same as R1/R2" via path equality. Good.

What's **not** caught: `R1.fq.gz R2.fq.gz --passthrough r1.fq.gz` on a case-insensitive volume — these point to the same file but `PathBuf::eq` says no. Probably out of scope; flag for documentation.

### 1.10 Minor: empty FASTQ — sanity_check rejects it

Plan §Behavior edge case "Empty input files: All three EOF simultaneously → clean exit, zero-records report". But `FastqReader::sanity_check` (`src/fastq.rs:407-438`) **rejects** an empty file with a hard error ("seems to be completely empty"). Step 2 says to run `sanity_check` on the passthrough at startup. So an empty-passthrough run would bail at sanity_check, not reach the lockstep loop. That's the right behaviour, but the plan's edge case statement is wrong — there's no "clean exit, zero-records report" path because sanity_check rejects empty input before the trimming pipeline runs. Update the edge-case description.

---

## 2. Assumptions — surfaced and validated

### 2.1 Hidden assumption: `--paired` is the only multi-input mode, but multi-pair is explicitly disallowed for passthrough

Plan validation 1.ii rejects `--passthrough` + multi-pair (`input.len() != 2`). This is correct given the design, but the user-facing error message should hint why: "Currently --passthrough supports exactly one R1/R2 pair; multi-pair input with parallel passthrough lists is not yet implemented." Otherwise users with shell globs that match multiple samples will read "got N input files" and conclude they need to remove samples, when in reality they should split into separate invocations or wait for v2.

### 2.2 Hidden assumption: passthrough file uses the same Phred encoding as R1/R2

The plan never mentions this, but `--phred33`/`--phred64` doesn't affect the passthrough (it's never quality-trimmed). The `FastqReader` doesn't enforce a Phred encoding — it just reads ASCII lines. So this is fine, but worth a one-liner in Assumptions: "the passthrough file's quality encoding is irrelevant — its quality string is carried verbatim."

### 2.3 Hidden assumption: passthrough file's records are well-formed

`FastqReader::sanity_check` only checks the first record. After that, a passthrough file with corrupt records (e.g. truncated mid-record) will produce a parse error from `next_record()` deep in the lockstep loop. The plan should clarify that mid-file FASTQ-parse errors on the passthrough propagate cleanly via `?` and don't risk inconsistency between R1/R2 outputs and the passthrough output. They probably **do** risk inconsistency: if R1/R2 batches are already mid-flight in workers when a passthrough parse error hits the reader thread, R1/R2 outputs may end up larger than the passthrough output before the error propagates. This is the kind of failure mode that "hard error, no output files written (or empty/partial files)" in Validation §3 papers over — for a malformed passthrough mid-file, the output files **will** be partial, not empty.

This is probably acceptable (users get a clear error and can re-run), but make it explicit in the plan: "On a passthrough read error mid-stream, the output files may be partial. The error is surfaced via `?` propagation through the reader thread → join → result. Recommend documenting this in the `--help` text."

### 2.4 Implicit assumption: order of input files in `--paired --passthrough X R1 R2`

What's the positional order? `trim_galore --paired --passthrough I1.fq.gz R1.fq.gz R2.fq.gz` (passthrough comes from a flag, R1/R2 are positional)? That's the plan's intent (the example in Goal section), but it should be locked down. The plan also implicitly assumes `cli.input.len() == 2` after this — which is correct because `--passthrough` is a flag-with-value, not a positional. Fine, but verify the clap parser doesn't accidentally consume R1 as the value of `--passthrough` if the user writes `--passthrough R1.fq.gz` without an `=` sign (it should — but worth a regression test like `trim_galore --paired --passthrough I1.fq.gz R1.fq.gz R2.fq.gz` runs through Cli::parse_from without I1 being eaten).

### 2.5 Implicit assumption: `--passthrough` is non-repeatable

`#[clap(long = "passthrough")] pub passthrough: Option<PathBuf>` is single-value. Future "v2 multi-pair" would need `Vec<PathBuf>` instead. The plan should note that the v1 → v2 transition will break source compatibility on this struct field. Cheap to defer, but flag it.

---

## 3. Efficiency

### 3.1 Channel pressure

The plan claims (§Efficiency) "no expected throughput regression" because passthrough writes are independent. Verify: each batch now ships **3** vecs of FastqRecord through the channel (R1, R2, passthrough), at `2 * sync_channel(2)` capacity. The reader thread blocks on `send` faster than before. Negligible in practice, but if the user has a small passthrough file (say barcodes are short ~50 bp records), the third vec adds another ~50 KB per batch — well under the 1.2 MB R1/R2 vecs. No concern.

### 3.2 Memory growth

Plan's claim "batch grows from 2.4 MB to 3.6 MB" assumes ~300 bytes per passthrough record. Cell-barcode reads are typically 16–28 bp, so ~80 bytes/record — the actual growth is closer to 2.4 MB → 2.7 MB. Plan over-estimates the memory cost (which is fine, conservative).

### 3.3 GzEncoder count per batch grows from 4 to 5 (or 2 to 3, with no retain_unpaired)

Currently `process_paired_batch` instantiates 2–4 `GzEncoder` instances per batch (R1, R2, optionally unpaired_r1, unpaired_r2). Adding passthrough makes that 3–5. Each `GzEncoder` allocates its own zlib-rs internal buffers (~256 KB by default for typical levels). With `--cores 16`, that's 16 workers × 5 encoders × ~256 KB = ~20 MB of encoder buffers in flight, up from ~16 MB. Fine, but worth noting in the efficiency section so the user can budget for `--memory`.

### 3.4 Sync check cost — plan claims 30 ms / million reads

Plausible but undermeasured. `split_ascii_whitespace()` allocates nothing and is O(prefix length). For 100-byte Illumina headers with ID prefixes ~30 bytes, each call is ~30 ns. Three calls per record × 1 M = 90 M ns = 90 ms. Plan claim of 30 ms is too optimistic by 3x but still well under 0.1% overhead. Update the number for honesty, not because it matters.

### 3.5 BTreeMap memory if passthrough adds large vec

The ordered-flush `BTreeMap<u64, PairedBatchResult>` holds out-of-order batches in memory. With BATCH_SIZE = 4096 and cores = 16, in the worst case ~32 batches might queue up (workers race ahead of main thread). Each `PairedBatchResult` now holds compressed R1 + R2 + passthrough; the passthrough adds ~50 KB compressed per batch. 32 × 50 KB = 1.6 MB transient. Negligible.

---

## 4. Validation sufficiency

### 4.1 §2 is the right invariant — but it's incomplete

The serial-vs-parallel parity test (§2) is correct as the load-bearing invariant given there's no Perl 0.6.11 baseline. But:

- **Test as written asserts decoded-byte-identity of passthrough output, but doesn't assert it for R1/R2 outputs in the same call.** That's a regression on the test's value — if a future bug makes the parallel R1 output differ from serial, the passthrough parity check would pass but the broader parity would silently regress. Add R1+R2 decoded-record equality to the same test.

- **The test doesn't exercise `--retain_unpaired` rejection** — i.e., it doesn't lock down that the worker hits the `debug_assert!(false, "Unpaired outcome impossible")` path on a malformed run. That's an internal contract; a negative test that confirms the validation block rejects it (Validation §5) is sufficient.

- **The test doesn't cover multi-batch behaviour** — 50 records all fit in one batch (`BATCH_SIZE = 4096`). The whole point of the parallel path is ordering across batches. The fixture should have **at least 8192 records** so the test exercises ≥2 full batches plus a partial. Otherwise the ordered-flush `BTreeMap` reassembly logic is untested for passthrough. Concrete suggestion: bump the fixture to ~10K pairs (still tiny on disk), which forces ≥3 batches at `BATCH_SIZE = 4096`.

### 4.2 Missing test: passthrough with all pairs dropped

The plan covers "all pairs dropped" in the §Edge cases prose but no explicit test asserts the passthrough output is a valid (potentially empty) gzip member when 100% of pairs fail `--length 1000` or `--discard_untrimmed` (with reads that have no adapter). Add an explicit test.

### 4.3 Missing test: passthrough with `--cores 1` going through the serial path

Plan §Self-review says "`--cores 1`: Goes through the serial path (`run_paired_end`); same passthrough plumbing applies." The §1 test exercises `run_paired_end` directly. But the actual entry point is `main.rs::run_paired` which routes by `cli.cores > 1 || cli.clumpify`. There's no integration test asserting `cli.cores == 1 && cli.passthrough == Some(...)` produces the expected output. The §1 test calls `run_paired_end` directly, bypassing the dispatcher. Add a CLI-level integration test that runs `trim_galore --paired --passthrough I1.fq.gz R1.fq.gz R2.fq.gz` with default `--cores 1` and verifies all three outputs.

### 4.4 Missing test: gzip vs plain output

The fixture in §1/§2 is presumably gzipped (`.fq.gz`). Add a test for `--dont_gzip` to exercise the plain-output branch of `process_paired_batch` (which uses `Vec<u8>` directly instead of `GzEncoder`). The two paths diverge inside `process_paired_batch`; if a bug only lives in the plain branch, gzip tests miss it.

### 4.5 §3 truncation test should also cover the parallel path's threading

The plan says §3 calls `run_paired_end_parallel`. Good. But verify the test asserts the error message names the file (the plan's behaviour spec says it should). A weaker assertion like "Err contains 'passthrough'" is less specific than checking the exact filename.

### 4.6 §4 shuffled-ID test fragility

The plan says §4 creates `multiome_I1_shuffled.fq.gz` by permuting records. But: the parallel path batches before sync-checking. If the shuffle happens within a single batch and the reader thread doesn't catch the desync at record 0 (instead seeing a mismatch at record 1000 inside a batch), the worker may still process up to that record before the reader propagates the error. The test should assert (a) an error is returned, (b) the error message identifies the **first** mismatching record (not an arbitrary one). The plan says "Hard error at the first mismatching record" — verify the implementation honours that, because the parallel path's reader can see records ahead of where the sync check fires.

### 4.7 No test for the `(Some, None)` passthrough EOF arm in the parallel reader

§3 covers truncated passthrough but only via `run_paired_end_parallel`. Add an explicit test for the converse: passthrough has **more** records than R1/R2 (i.e., R1/R2 EOF first). The plan's Step 6.4 enumerates the error message ("passthrough has more records than R1/R2") but no test confirms the error fires correctly in the parallel path's `(None, None)` arm.

### 4.8 Suggested 6th test: round-trip through `main()` via CLI integration

The existing `cli.rs::tests` use `Cli::parse_from` but don't actually invoke `main()`. None of the proposed tests verify that `--passthrough` survives the CLI → validate → dispatch → run path end-to-end. Add a smoke test in `tests/integration_passthrough.rs` (top-level `tests/` dir) that builds the binary and runs it via `Command::new("trim_galore")` on the fixture, then verifies output files exist. Otherwise the wiring in main.rs Step 8 ("thread `cli.passthrough.as_deref()`...") is only verified by hand.

---

## 5. Alternatives

### 5.1 Header sync check: strip `/[123]` suffix in v1

Plan Q&A #2 leaves this as an "open" decision but defaults to "no strip". Counter-argument: legacy `@read/1` `@read/2` headers are still extremely common in archived data (SRA, ENA, every pre-CASAVA-1.8 dataset). Trim Galore's user base skews toward bisulfite / clock / RRBS — these datasets routinely come from older platforms. Stripping `/[123]` is a one-line change and zero downside (no false-positive risk because `/1` is never part of a real ID prefix). I'd argue for stripping by default in v1.

Counter-counter: if you strip and the user has weird headers like `@read/1.some_metadata`, stripping `/1` from the middle would corrupt the prefix. Safe approach: strip only **trailing** `/1` `/2` `/3` using `trim_end_matches`. That's still one line:

```rust
let trimmed = id.strip_prefix('@').unwrap_or(id);
let head = trimmed.split_ascii_whitespace().next().unwrap_or("");
head.trim_end_matches(|c: char| matches!(c, '1' | '2' | '3'))
    .trim_end_matches('/')
```

Two-line, robust. Worth doing in v1.

### 5.2 Positional vs flag CLI shape

The plan uses `--passthrough <FILE>` (flag). Alternative: a positional, e.g. `trim_galore --paired R1.fq.gz R2.fq.gz I1.fq.gz` where 3 inputs in `--paired` mode is the trigger. Rejected — it conflicts with multi-pair `--paired` semantics (`R1 R2 R1 R2 ...`) and confuses the validator (3 inputs is already an "odd count" error). The flag is the right shape.

What about `--carry` as an alias? `--passthrough` is unambiguous but long. Consider adding `--carry` as a short alias for ergonomic discoverability. Optional.

### 5.3 Stats: new struct vs extending `PairValidationStats`

Plan extends `PairValidationStats` with 3 new u64 fields. Alternative: introduce a new `PassthroughStats` struct, returned as a 4th element of the tuple from `run_paired_end` / `run_paired_end_parallel`. This keeps `PairValidationStats` semantically tight (only pair-validation state, not passthrough state).

Counter: returning 4 elements from a function already returning 3 ripples through every caller. Adding to `PairValidationStats` is less disruptive but conflates two concerns. I'd lean toward a wrapper:

```rust
struct PassthroughStats {
    records_checked: u64,
    records_kept: u64,
    records_dropped: u64,
}
// then in run_paired_end:
pub fn run_paired_end(...) -> Result<(TrimStats, TrimStats, PairValidationStats, Option<PassthroughStats>)>
```

`Option<PassthroughStats>` makes the "feature disabled" case unambiguous — no `passthrough_records_checked == 0` ambiguity (which could also mean a zero-record empty input that did go through passthrough). Cleaner.

This is a real ergonomic improvement and worth doing now rather than refactoring later.

### 5.4 FastQC integration

Plan Step 10 says FastQC runs on all three outputs. Consider: do users want FastQC on the passthrough file? It's an index/barcode read — the per-base quality plot is meaningless because all bases are typically the same (the cell barcode). FastQC will produce a valid report but it'll be useless. Two options:

- Run FastQC on all three (plan's choice) — generates an extra report users probably ignore. Cost: ~50 ms × output_size.
- Skip FastQC on passthrough by default; add `--fastqc-passthrough` to opt in.

Option 2 is more user-friendly but adds CLI surface. Option 1 is simpler. Either is defensible — flag this for the implementer.

### 5.5 Output naming: `_passthrough` suffix

Plan uses `<stem>_passthrough.fq(.gz)`. Alternatives:
- `<stem>_carried.fq.gz` — shorter
- `<stem>.passthrough.fq.gz` — dot-separated, more like `trimmomatic`'s convention
- Keep the original filename verbatim with `_unchanged` suffix

The `_passthrough` suffix is reasonable. Just lock the name early — users will write Snakemake rules against it and renaming later is a breaking change.

---

## 6. Action items (prioritized)

### Critical (must fix before implementation)

1. **C1**: Add `PartialEq` (and probably `Eq`) to `PairValidationStats` (`src/report.rs:106`). Validation §2's parity assertion as written won't compile otherwise. Make this an explicit Step-3 sub-item in the plan. (See §1.1.)
2. **C2**: Spell out the worker-thread interleaving of R1/R2/passthrough in `process_pairs`, including the `passthrough: &[FastqRecord]` (immutable) shape vs the existing `&mut [FastqRecord]` for r1/r2. The current plan signature for `process_pairs` doesn't show the third slice's mutability. (See §1.2.)
3. **C3**: Decide between `Option<Vec<u8>>` and "empty `Vec<u8>` means disabled" for `compressed_passthrough` in `PairedBatchResult`. Currently the plan is inconsistent vs existing `compressed_unpaired_*` fields. (See §1.3.)
4. **C4**: Add Step 6.4.1 to the plan: handle the EOF arm of `read_pairs_round_robin` to drain passthrough and detect "passthrough longer than R1/R2". (See §1.4.)

### Important (close before implementation start, but don't block the design)

5. **I1**: Make the §2 parity test fixture ≥ 8192 records so multi-batch ordering is exercised. Also assert R1+R2 decoded-byte-identity in the same test, not just passthrough. (See §4.1.)
6. **I2**: Add a test for `--cores 1` going through `main()` end-to-end via a CLI integration test (or via the `run_paired` dispatcher), not just direct `run_paired_end` invocation. (See §4.3.)
7. **I3**: Add a test for `--dont_gzip` passthrough output. The plain branch of `process_paired_batch` diverges from the gzip branch. (See §4.4.)
8. **I4**: Add a test for passthrough-longer-than-R1/R2 (the converse of §3). (See §4.7.)
9. **I5**: Use a positional alternative or document the parse shape — confirm clap doesn't eat R1 as `--passthrough`'s value when no `=` is used. Add a unit test for the parse. (See §2.4.)
10. **I6**: Strip trailing `/[123]` in `read_id_prefix()` for v1 — legacy headers are extremely common in archived data. Two lines. (See §5.1.)
11. **I7**: Make the parallel reader's flush loop write R1+R2+passthrough atomically per `expected` increment. Plan implies this but doesn't say it. (See §1.6.)
12. **I8**: Update the §Edge case "Empty input files" to reflect that `sanity_check` rejects empty input — there's no "clean exit zero-records" path. (See §1.10.)

### Optional (nice-to-have / future-proofing)

13. **O1**: Consider `Option<PassthroughStats>` as a 4th return value instead of extending `PairValidationStats`. Cleaner separation of concerns. (See §5.3.)
14. **O2**: Decide FastQC behaviour for the passthrough output — run by default, or skip with opt-in `--fastqc-passthrough`? (See §5.4.)
15. **O3**: Update the efficiency claim from "30 ms / 1 M reads" to ~90 ms / 1 M reads for honesty (still well under 0.1%). (See §3.4.)
16. **O4**: Document mid-stream passthrough parse errors leave partial output files. Mention in `--help` text. (See §2.3.)
17. **O5**: Add `--carry` as a short alias for `--passthrough` for ergonomics. (See §5.2.)
18. **O6**: Validation 1.ii error message: hint that multi-pair passthrough is "not yet implemented" rather than just "got N files". (See §2.1.)

---

## Summary

The plan is **architecturally sound** and aligns well with existing patterns (sink-agnostic `classify_paired`, ordered-flush `BTreeMap`, output-collision preflight). The headline gap is the missing `PartialEq` on `PairValidationStats` (C1), which would have caused the load-bearing parity test to fail compilation on day 1. The worker-thread plumbing for the third stream needs more concrete spec (C2). Validation coverage has real gaps around multi-batch ordering, plain output, and CLI-level integration — fix them or accept blind spots. Otherwise this is a 200-LoC feature with a clean blast radius and low regression risk on existing flags.
