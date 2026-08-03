# CODE review — Reviewer A — #363 mixed-format pair message

**Branch:** `fix/mixed-format-pair-message` (uncommitted working tree) off `dev` @ `a4ffd47`
**Reviewed:** `src/format.rs`, `src/main.rs`, `src/clump_only.rs`, `tests/integration_paired_format_guard.rs` (untracked), `CHANGELOG.md`
**Plan:** `plans/mixed-format-pair-message/PLAN.md` v2 (§13 treated as a claim, verified independently)

## Summary

The implementation is sound. The `PairedShape` design does what the plan says it does, the guard fires on exactly the intended set of invocations, the `ClumpOnlyFastqOut` message is byte-identical to the one it replaces, the BAM-count predicate correctly accepts plain+gz pairs, and none of the four pinned substrings is broken by a line continuation. I verified every one of these by running the binary, not by reading the literals.

Verification performed:

- `cargo build --release`, `cargo test --release` → **441 passed, 0 failed** (369 lib + 12 + 15 + 1 + **11 new** + 2 + 8 + 23).
- `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` clean (forced a re-check by `touch`ing the four changed files first, so this is not a cached result).
- **Byte-identity invariant (V11) confirmed empirically, not inferred.** Built `dev` @ `a4ffd47` in a throwaway worktree with a separate `--target-dir`, ran both binaries over PE, SE, `--hardtrim5`, `--clock`, `--demux` and `--clump_only`, and md5'd all 15 resulting `*.fq`/`*.fq.gz` outputs: **identical**. Worktree removed afterwards.
- 17 manual binary invocations covering all four `PairedShape` values × mixed/two-BAM, both argument orders, multi-pair, all four specialty modes, `--clumpify`, `--passthrough`, `--demux`, `--retain_unpaired`, N=1 interleaved uBAM, plain+gz.

Two problems found and fixed (one Medium, one Low). Two recommendations worth the maintainer's attention, one of which is a **pre-existing docs statement that flatly contradicts the error message this PR is polishing**.

---

## Point-by-point on the items flagged for scrutiny

### Guard activation condition — correct

`src/main.rs:286-301`. Verified by running the binary, one invocation per case:

| Invocation | Result | Correct? |
|---|---|---|
| `--paired <fq> <bam>` | mixed message | ✅ |
| `--paired <bam> <fq>` | mixed message, labels swapped | ✅ |
| `--paired <bam1> <bam2>` | two-BAM message | ✅ |
| `--paired --output-format ubam` × both | shape-appropriate, remediation echoes `--output-format ubam` | ✅ |
| `--clump_only --paired` (FASTQ out) × mixed and two-BAM | aux-tag message both times | ✅ |
| `--clump_only --paired --output-format ubam` × both | mixed vs two-BAM distinguished, remediation echoes both flags | ✅ |
| `--clumpify --cores 2 --paired <bam1> <bam2>` | two-BAM message (guard does fire) | ✅ |
| `--retain_unpaired --paired <fq> <bam>` | mixed message | ✅ |
| `--hardtrim5 10 --paired <fq> <bam>` | exit 0, writes both `*.10bp_5prime.fq` | ✅ unchanged |
| `--clock --paired <fq> <bam>` | `Paired-end files have different numbers of reads!` | ✅ unchanged |
| `--implicon --paired <fq> <bam>` | same read-count error | ✅ unchanged |
| `--paired <interleaved.bam>` (N=1) | exit 0, `_val_1`/`_val_2` written | ✅ not caught |
| `--paired --passthrough <i1> <fq> <bam>` | `--passthrough is not supported with uBAM input` | ✅ keeps its own earlier message |
| `--demux <sheet> --paired <fq> <bam>` | `Demultiplexing is only allowed for single-end files` | ✅ `Cli::validate` wins |

No paired entry point is missed and no non-paired one is wrongly caught. `--demux` cannot reach the guard: `cli.rs:922-925` rejects `--demux` + `--paired` inside `Cli::validate()`, which runs at `main.rs:166`. `--passthrough` is rejected against any BAM at `main.rs:253`, i.e. 33 lines earlier — precedence preserved as §2.4 intends.

Paired dispatch coverage is complete. The four paired sinks are `run_paired` (`main.rs:773`), `run_paired_ubam_single_file` (`:675`, N=1 only), `run_ubam_output` (`:666`) and the clump-only branches (`:471`, `:535`). Every N>1 non-specialty route passes through `main.rs:286`.

### `ClumpOnlyFastqOut` asymmetry — correctly implemented, message faithful

`src/format.rs:169-176`. The arm keys on `n_bam >= 1` implicitly: the loop `continue`s only on `n_bam == 0`, so reaching the `matches!(shape, ClumpOnlyFastqOut)` branch already means at least one BAM. Confirmed live for all three shapes (fq+bam, bam+fq, bam+bam) — all three produce the aux-tag message.

I diffed the literal against `src/clump_only.rs:403-408` character by character: identical, including the parenthetical and the `Input: {}` suffix. The "which file to name" rule also matches — the original loops `[input_r1, input_r2]` and reports the first BAM; the new code picks `if is_bam(&fmts[0]) { r1 } else { r2 }`. Same file in every case.

Note the original check at `clump_only.rs:399-410` was **not** deleted, so this path now has two identical guards with `main()`'s firing first. Harmless and arguably good (it protects `clump_only_paired` as a library entry point).

### The predicate is BAM count, not format equality — confirmed

`src/format.rs:141` counts `UnalignedBam` and nothing else. `--paired <plain.fastq> <R2.fastq.gz>` exits 0 and writes `plain_R1_val_1.fq` — verified by running the binary and by `plain_plus_gzip_fastq_pair_is_still_accepted`. The unit test `pair_guard_accepts_plain_plus_gzip_fastq` pins it at the function level too.

### The two belt-and-braces checks — reachable and correct, but narrower than what they replaced

`main.rs:2011-2020` and `clump_only.rs:947-957`. Both compare `matches!(fmt_r1, UnalignedBam) != matches!(fmt_r2, UnalignedBam)`, i.e. they catch exactly the hazard described: a mixed pair would take `source_header` from R1 alone while opening each side by per-file detection, emitting a BAM that mixes FASTQ-derived and BAM-derived records. That hazard is guarded.

However, the code each one replaced was **broader**. The deleted loop in `run_ubam_output` (`dev`'s `main.rs:1748-1761`) and the deleted block in the clump-only Shape A branch (`dev`'s `:536-558`) both rejected *any* BAM, so a two-BAM pair was also stopped at the leaf. The new backstops let a two-BAM pair through. See Low-2 — it is not currently reachable and the outcome would probably be acceptable output rather than corruption, but the stated intent ("so a future addition to the exemption list cannot reopen the silent-wrong-output path") is only half met.

### The internal-invariant `bail!` in `run_paired` — correctly placed

`main.rs:1287-1296` sits above the `if cli.cores > 1 || cli.clumpify` split at `:1298`, so it covers the parallel branch and the sequential `FastqReader::open` branch at `:1338-1339` alike. Placement is right.

### Rust string-continuation hazard — none of the four substrings broken

Verified by running the binary and reading the emitted single-line messages, not by inspecting the literals:

```
Error: --paired with two BAM files is not supported. uBAM paired mode expects a single
interleaved file: `trim_galore --paired interleaved.bam`. Got two BAM files: A and B.
Combine them into one mate-adjacent file first: `samtools merge -n -o interleaved.bam A B`.

Error: --paired requires both inputs of a pair to be the same format. Got mixed: A is
FASTQ (plain) and B is uBAM. Pass two FASTQ files, or a single interleaved uBAM. If you
meant two FASTQ files, check for a mis-typed filename.
```

`two BAM files is not supported`, `uBAM paired mode expects`, `single interleaved` and `same format` all appear intact. `git diff --stat -- tests/` is empty, so the three pre-existing tests pass unedited — confirmed by the suite run.

### Retired code — genuinely unreachable

`main.rs:536-541` (comment replacing `dev`'s `:498-512`). The retired condition was `clump_only ∧ output_format==UBam ∧ paired ∧ N==1 ∧ ¬UnalignedBam`. The guard at `main.rs:259-265` fires on `paired ∧ N==1 ∧ ¬UnalignedBam` and runs ~240 lines earlier in `main()`, before the clump-only dispatch at `:447`. Strict subset, so unreachable. No test pinned its wording (`grep`ed `tests/` for its distinctive phrases — no hits).

### Message quality — good, with one editorial nit

Plain declarative register throughout, no filler, no second-person scolding. `If you meant two FASTQ files, check for a mis-typed filename.` is second-person but imperative-helpful, not scolding — reads like the rest of the codebase's error text.

Both the single-pair and multi-pair forms are grammatical: *"Got mixed: A is FASTQ (plain) and B is uBAM."* / *"Pair 2 of 2 is mixed: A is … and B is …"*. The `"Got "` prefix (deviation 2 in §13) was the right call — the empty prefix would have produced a fragment.

I also checked that every remediation the messages offer is actually accepted by the binary, per shape:

| Shape | Advice given | Accepted? |
|---|---|---|
| `Trim` | `--paired interleaved.bam` | ✅ exit 0 |
| `TrimUbamOut` | `--paired --output-format ubam interleaved.bam` | ✅ |
| `ClumpOnlyUbamOut` | `--clump_only --paired --output-format ubam interleaved.bam` | ✅ (also two FASTQ ✅) |
| `Trim` + `--clumpify` | `--paired interleaved.bam` | ✅ (`--clumpify --paired interleaved.bam` exits 0) |
| `ClumpOnlyFastqOut` | add `--output-format ubam` | ✅ |

No shape is handed a command the binary rejects — the class of bug #363 is about.

Nit: `Got mixed:` is telegraphic. `Got a mixed pair:` reads better and still satisfies the pinned `contains("mixed")`. Maintainer's call.

### Test quality — good; the named trap is closed

`tests/integration_paired_format_guard.rs:46-50` `nonexistent_out` asserts its own precondition (`!p.exists()`) before handing the path to `-o`, so the four `!out.exists()` assertions are real. The vacuity trap the plan calls out is closed, and closed the right way (assertion, not comment).

The two negative stderr assertions in `mixed_pair_rejection_has_no_side_effects` are non-vacuous — I confirmed both target strings exist (`main.rs:1036` `"Auto-detecting adapter type..."`, `main.rs:1252` `"Trimming (paired-end):"`) **and** that an accepted paired run prints both (grep count = 2). If the guard were moved back below dispatch, both would fail.

Tag collision check: all 12 `tempdir(tag)` values are distinct, and the `tg_pfg_` prefix does not collide with `tg_int_ubam_out_*` / `tg_clump_only_ubam_*` used by the other test binaries.

One test is near-vacuous — see Low-4.

---

## Fixes applied

Both re-verified afterwards with `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings` and the full `cargo test --release` (441 passing).

### FIX 1 (Medium) — `CHANGELOG.md` asserted a bug that never shipped

`CHANGELOG.md:152-157` (pre-fix) claimed:

> The two-BAM message's remediation was also corrected. It suggested `samtools collate -O r1.bam r2.bam`, which does not combine two files …

**This is false as a statement about released behaviour, and I confirmed it rather than suspected it.** `git show a4ffd47:src/main.rs` — the three shipped two-BAM messages (`:543`, `:1271`, `:1754`) end with *"Got two BAM files; one of them is {}."* and contain **no `samtools` suggestion at all**. The `samtools collate` hint the entry attributes to them lives in `src/bam.rs:57`, in a completely different error (`GROUPED_INPUT_ERR`, for a single interleaved BAM with grouped rather than adjacent mates) — and there it is in the **correct** one-input form, `samtools collate -O input.bam tmp > interleaved.bam`.

The destructive two-input form only ever existed in v1 of the plan. Shipping this paragraph would tell users TrimGalore once printed a silently-lossy command, and could lead a reader to distrust the *correct* hint in `bam.rs:57`.

Rewrote the paragraph to describe what actually changed (the message now names a combine command where it previously named none) and to keep the `collate` caveat as a forward-looking note rather than a claim about history. Also fixed the last sentence, which said `--clock`/`--implicon` *"accept mixed pairs today"* — they do not; both fail with `Paired-end files have different numbers of reads!` (verified live), which is precisely why §2.6/A4 files them as out-of-scope rather than correct.

### FIX 2 (Low) — internal-invariant message named the wrong function

`src/main.rs:1290`. The backstop inside `fn run_paired` (`main.rs:1238`) said:

```
"internal error: uBAM input reached run_paired_end ({}); …"
```

`run_paired_end` is a different function — `trimmer::run_paired_end` at `src/trimmer.rs:381`, one of the two public trimmer entry points named in `CLAUDE.md`. A maintainer or bug reporter following this string lands in the wrong file. The sibling backstops name their host functions correctly (`run_ubam_output_paired_two_files`, `clump_only_paired_to_bam_one_pair`), so this one was the odd one out. Changed to `run_paired`. No test pins the string (grepped `tests/` and `src/`).

---

## Recommendations

### High-1 — `docs/src/content/docs/quickstart.md:57` contradicts the message this PR is polishing

Pre-existing on `dev`, not introduced here, so strictly outside the diff — but it is the same subject and a user reading the quickstart is exactly the user who hits this error.

```
uBAM input is auto-detected — no flag needed. Paired reads may come as two BAM
files or a single interleaved BAM (samtools `sort -n` / `collate` / Picard /
fgbio order):
```

"Paired reads may come as two BAM files" is false — the binary rejects it, and this PR makes that rejection more emphatic. The code block immediately below only shows `trim_galore --paired interleaved.bam`, so the sentence is a leftover. Verified by reading the file and by running the invocation it promises.

Suggested minimal correction, in the existing register:

```
uBAM input is auto-detected — no flag needed. Paired reads must arrive as a
single interleaved BAM with mates adjacent (samtools `sort -n` / `collate` /
Picard / fgbio all produce this order); two separate BAM files are not
supported and are rejected up-front.
```

Not applied: it is maintainer-owned prose outside the reviewed diff, and #363 is arguably the right issue to close it under rather than a silent drive-by.

### Medium-1 — `debug_assert_eq!` leaves the guard's own precondition unchecked in release

`src/format.rs:130-134`. `[profile.release]` in `Cargo.toml:58-62` does not set `debug-assertions`, so it defaults off and this assertion compiles out of every shipped binary. Two distinct failure modes if `inputs.len() != formats.len()` ever holds:

1. **Silent skip** (the worse one). `inputs.chunks(2).zip(formats.chunks(2))` truncates to the shorter iterator. A short `formats` means later pairs are never examined and the function returns `Ok(())` — a validation guard that passes *without having checked*. No panic, no message.
2. **Panic.** An odd-length `inputs` gives a final 1-element chunk and `paths[1]` panics with index-out-of-bounds instead of erroring cleanly.

Neither is reachable today: `input_formats` is built by mapping over `cli.input` (`main.rs:184-188`), and `Cli::validate()` → `validate_paired_input("Paired-end")` (`cli.rs:596`, under a plain `if self.paired`) rejects odd N for every paired mode including `--clump_only`. I confirmed both. So this is defence-in-depth, not a live bug.

Given that the entire purpose of this function is to not let a bad pair through silently, a hard check is cheap insurance:

```rust
anyhow::ensure!(
    inputs.len() == formats.len() && inputs.len().is_multiple_of(2),
    "internal error: reject_bam_format_mismatch_in_pair got {} inputs and {} formats",
    inputs.len(),
    formats.len()
);
```

and then `chunks_exact(2)`. Recommending rather than fixing: the plan explicitly specified `debug_assert` (§5, A7), so changing it is a design decision, not a defect fix.

### Low-1 — belt-and-braces checks are narrower than the code they replaced

Detailed above. Both new backstops test *mismatch*; the deleted code tested *any BAM*. If a future exemption ever let a two-BAM pair reach `run_ubam_output_paired_two_files` or `clump_only_paired_to_bam_one_pair` Shape A, it would pass. The consequence is milder than for a mixed pair (both sides would be real BAM readers, so records are consistent; only R2's `@HD`/`@PG`/tag dictionary is discarded), which is arguably why it was not flagged — but if the goal is "the invariant is enforced, not documented", the enforced predicate should match the documented one. One extra clause, or reuse `reject_bam_format_mismatch_in_pair` itself with the leaf's own shape.

### Low-2 — §6's "less I/O on all runs" is two-thirds true

`main.rs:1287-1296` retains a per-pair `detect_input_format` call on both inputs for **every** accepted paired FASTQ run — it is the backstop, so this is deliberate and correct. But the plan's §6 claim 2 and the general framing ("*removes* I/O") describe three deleted re-detections; only two were deleted. The cost is two `open` + one BGZF-block decompress per pair against a run that trims millions of reads, so this is a documentation nit, not a performance finding. Worth a one-line correction in `PLAN.md` §13 if the plan is kept as the record.

### Low-3 — `PLAN.md` §13/V8 misdescribes the multi-pair test

`tests/integration_paired_format_guard.rs:336-344` uses **three** distinct paths (`R1`, `R2`, `R2`, `bam`), not the four the plan's V8 row insists on. The reuse is benign here — pairs `(R1,R2)` and `(R2,bam)` are not duplicates, within-pair R1≠R2 holds, and the four planned outputs (`R1_val_1`, `R2_val_2`, `R2_val_1`, `ubam_test_val_2`) do not collide, so the collision pre-flight never fires and the test does exercise pair indexing rather than precedence. I confirmed independently with four genuinely distinct paths (`BS-seq_10K_R1.fastq.gz`, `BS-seq_10K_R2.fastq.gz`, `clock_10K_R1.fastq.gz`, `ubam_test.bam`) that the message is still `Pair 2 of 2`. **No test change needed** — only the plan's claim is loose.

### Low-4 — one near-vacuous assertion

`tests/integration_paired_format_guard.rs:178-181`:

```rust
!stderr.contains("samtools collate -O test_files/ubam_test.bam")
```

This guards a string that never existed in any shipped message (see FIX 1) and that can only reappear if someone types it deliberately, fixture path and all. `!stderr.contains("samtools collate")` would be a genuine, still-true assertion of the intended property ("do not suggest collate for combining two files") and would survive a fixture rename.

### Low-5 — undocumented deviation: `{fmt_hint}` dropped from the mixed message

`PLAN.md` §3.4 specifies `Pass two FASTQ files, or a single interleaved uBAM{fmt_hint}.`; the implementation (`src/format.rs:181`) has no `{fmt_hint}`. Not listed among §13's four deviations. Behaviourally harmless — I verified the bare advice is valid for all three shapes that emit it — but §13 is otherwise a careful record, so it is worth adding.

### Low-6 — style: `anyhow::bail!` in `clump_only.rs`

`src/clump_only.rs:948`. The module imports `bail` at `:37` and the adjacent Shape B backstop at `:914` uses the bare `bail!`. Use `bail!` for consistency. (`main.rs` has no such import, so `anyhow::bail!` there is correct.)

---

## Things I checked and found sound — no action

- **CHANGELOG placement.** Appended to the single existing `#### Fixes` block at `CHANGELOG.md:84` under `### Unreleased`; no second block introduced. The three enumerated observable consequences are each accurate (verified: no output dir created, no `Auto-detecting adapter` line, no partial multi-pair output).
- **No `-D warnings` exposure from the deletions.** `InputFormat` and `detect_input_format` remain used at `main.rs:28`, `:49`, `:66`, `:187`, `:259`, `:2009-2010`, `:2022`. Clippy clean on a forced re-check.
- **Efficiency.** One pass over an already-materialised `Vec<InputFormat>` of `Copy` elements with early exit, O(N/2). `pair_prefix` is `format!`-allocated before the `ClumpOnlyFastqOut` arm that does not use it — an error path, so irrelevant.
- **No external runtime dependency introduced.** The `samtools merge -n` hint is text in an error message; nothing shells out.
- **uBAM detection stays content-based.** The guard consumes `detect_input_format` output and adds no filename logic of its own.
- **N=1 routing.** `--paired` + N=1 goes to `run_paired_ubam_single_file` (`main.rs:673-685`) or `run_ubam_output_paired_single_file` (`:1752`), never to `run_paired`, so the backstop's blanket "any BAM" rejection cannot catch the legal interleaved case. Verified live.
- **`run_ubam_output`'s deleted loop.** The N=1 branch returns at `main.rs:1752-1762` before the `if cli.paired` block, so deleting the loop cannot expose an N=1 path.
- **Unit tests.** Nine tests in `format.rs`, all asserting real content including both negative properties (`!contains("uBAM paired mode expects")` on mixed, `!contains("interleaved.bam")` on `ClumpOnlyFastqOut`). `paths(n)`/`err()` helpers keep them readable without fixtures. Placement in `format.rs` rather than `main.rs` is right — `main.rs` still has zero `#[cfg(test)]` modules (`cargo test` reports `0 passed` for the bin target).
