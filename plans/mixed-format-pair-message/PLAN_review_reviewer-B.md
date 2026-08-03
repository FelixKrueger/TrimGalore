# PLAN review — Reviewer B

**Plan:** `plans/mixed-format-pair-message/PLAN.md` (#363)
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev` @ `a4ffd47`, worktree clean for `src/`, `tests/`, `CHANGELOG.md`
**Method:** every load-bearing `file:line` claim re-read in the source; ten invocations reproduced against a freshly-built `./target/release/trim_galore`; scratch output written outside the repo. Where I state "verified live" I ran the binary; where I state "read" I read the code; anything I could not confirm is labelled as such.

**Verdict:** the plan is unusually well-grounded — its factual claims about the two defective sites, the guard position, `Cli::validate()`'s preconditions, and the existing test assertions all check out. But it contains **one real regression** (a paste-and-fail remediation command for `--clump_only` on the FASTQ-output path) that its own validation table would not catch, and it leaves two implementation details unspecified (`{fmt}` rendering, helper placement) whose obvious resolutions are wrong or off-convention.

---

## 1. What I confirmed correct

Stated up front so the findings below are read in proportion.

| Plan claim | Status | How |
|---|---|---|
| §2.1 — `main.rs:1268-1278` and `:1750-1762` bail on the *first* BAM, so a mixed pair gets the two-BAM message | **Confirmed** | Read both sites; reproduced live in both argument orders |
| §2.1 — rejection happens after adapter detection, poly-G detection, the banner, and `ensure_output_dir` | **Confirmed** | Live: `-o outA` left the directory behind and printed all four |
| §2.1/§2.2 — these are the **only** two defective sites | **Confirmed** | `grep -rn 'two BAM files\|same format\|mixed' src/` → `main.rs:1271`, `:1754` (defective) and `:543`/`:553` (the correct clump-only pair). No other paired-input format guard emits either message |
| §2.2 — `main.rs:537-557` is already correct, not a defect | **Confirmed** | Read `:536-558`; live `--clump_only --paired --output-format ubam <fq> <bam>` → *"requires both input files to be the same format. Got mixed: A and B."* |
| §2.3 — `:264` sits after the N=1 check (`:257-263`) and before `gzip` (`:272`) / `ensure_output_dir` (`:277`) | **Confirmed** | Read `:250-280`. Nothing at all between `:263` and `:272` |
| §2.3 — the guard precedes every dispatch branch | **Confirmed** | hardtrim5 `:321`, hardtrim3 `:346`, clock `:371`, implicon `:385`, clump-only `:409`, uBAM-out `:653`, N=1 paired trim `:661`, paired trim `:675`, SE `:785` — all below `:277` |
| §2.3 — `input_formats` already exists at `:182-186`, guard needs no new I/O | **Confirmed** | Read; comment at `:178-181` says exactly what the plan quotes |
| §2.4 — all three existing assertions are loose substrings satisfied by §3.3 | **Confirmed** | Read `integration_clump_only_ubam.rs:397-403` (`"single interleaved" ‖ "uBAM paired mode expects"`), `:614-618` (`"same format" ‖ "mixed"`), `integration_ubam_out.rs:459-463` (`"two BAM files is not supported" ‖ "interleaved file"`). §3.3's strings contain all five |
| §2.4 — nothing covers site `:1268` | **Confirmed** | No test in `tests/` invokes `--paired <bam> <bam>` without `--output-format ubam` |
| §2.5 — `--hardtrim5 10 --paired <fq> <bam>` silently accepts | **Confirmed live** | Wrote both `BS-seq_10K_R1.10bp_5prime.fq.gz` and `bam1.10bp_5prime.fq.gz`. (Mechanism: `:321-345` loops `for input in &cli.input` and never consults `cli.paired`, so pairing is genuinely inert there — the plan's characterisation is right) |
| §2.5 — `--clock --paired <fq> <bam>` → *"Paired-end files have different numbers of reads!"* | **Confirmed live** | Exactly that message |
| §3.1 — field types | **Confirmed** | `cli.rs:372` `hardtrim5: Option<usize>`, `:377` `hardtrim3: Option<usize>`, `:383` `clock: bool`, `:391` `implicon: Option<usize>`, `:208` `clump_only: bool`, `:89` `paired: bool` |
| §3.3 — mixed-batch-across-a-batch is a shipped feature that the wording must not contradict | **Confirmed** | `integration_ubam_out.rs:340 ubam_out_preserve_tags_mixed_batch_allowed` is single-end, so a paired-only guard cannot touch it. "both inputs of **a pair**" is the right qualifier |
| A1 — content-based detection, `bgzip`-framed FASTQ → FASTQ | **Confirmed** | `format.rs:40-77`, decompress-and-probe for `BAM\1` |
| A2 — even count `cli.rs:506`, R1≠R2 `cli.rs:513`, N=1 carve-out `cli.rs:503` | **Confirmed**, and stronger than stated | Also verified the *reachability*: `validate_paired_input("Paired-end")` is invoked at `cli.rs:596` under plain `if self.paired`, so it covers `--clump_only --paired` too — which A2 relies on implicitly but does not say |
| A3 — `--paired` + N=1 legal only for uBAM, rejected otherwise at `:257-263` | **Confirmed** | Read; and `--paired <interleaved.bam>` with FASTQ output succeeds live (produced `ubam_paired_test_val_{1,2}.fq` + 4 reports), so §3.3's "or a single interleaved uBAM" hint is valid **on the plain trim path** |
| V10's "421 tests" | **Confirmed** | `grep -c '#\[test\]'` across `src/` + `tests/` sums to exactly 421 |
| §6 — no measurable cost | **Confirmed** | One pass over an in-memory `Vec` of ≤ N elements |

**Entry points the plan did not enumerate, which I checked and which are correctly non-applicable:**

- `--demux` — rejected with `--paired` at CLI level (verified live: *"Demultiplexing is only allowed for single-end files"*). Not a paired entry point.
- `--passthrough` — `main.rs:251-256` rejects it against *any* BAM input, i.e. **before** `:264`. So `--paired --passthrough pt.fq <fq> <bam>` keeps its own, more specific message. Correct precedence; see O5.
- `--retain_unpaired` — no format interaction (only `--output-format ubam` rejects it, `cli.rs:586`); mixed pairs route through the ordinary paired path and are covered.
- `run_specialty_paired` (`main.rs:2347`) — the shared multi-pair driver for `--clock`, `--implicon`, **and** `--clump_only --paired` FASTQ output. It has no format guard of its own; the clump-only leg's guard lives in the callee (`clump_only.rs:399-410`). This is where finding C1 comes from.

---

## 2. Logic review — findings

### C1 (Critical) — the remediation command is wrong for `--clump_only` on the FASTQ-output path

**Confirmed wrong, verified live.**

§3.1 deliberately includes `cli.clump_only` in the guard. §5 defines `fmt_flag` as `" --output-format ubam"` **when uBAM output is requested, else `""`**. Combine the two and the FASTQ-output clump-only case produces a remediation the tool itself rejects:

| Invocation (all FASTQ output) | Today (verified live) | Under the plan |
|---|---|---|
| `--clump_only --paired <bam1> <bam2>` | `uBAM input under --clump_only requires --output-format ubam (…would drop aux tags). Input: bam1.bam` | `--clump_only --paired with two BAM files is not supported. … expects a single interleaved file: `trim_galore --clump_only --paired interleaved.bam`. …` |
| `--clump_only --paired <fq> <bam>` | same message | `--clump_only --paired requires both inputs of a pair to be the same format. … Pass two FASTQ files, or a single interleaved uBAM. …` |

Both new messages are wrong in the same way. I verified live that `trim_galore --clump_only --paired <interleaved.bam>` — the *exact* command the two-BAM message would print — fails:

```
Error: --clump_only --paired requires two FASTQ input files. Single-file
interleaved uBAM input needs --output-format ubam (add `--output-format ubam`
to the command line).
```

(`main.rs:426-432`.) So the plan converts a message that names the one flag the user needs into a paste-and-fail command, on a fix whose entire premise (§1, goal 1) is *"must not be offered a remediation that does not apply"*. The mixed-pair variant is the same defect in softer form: "or a single interleaved uBAM" is not an available option in that mode.

Today's `clump_only.rs:403-409` message is genuinely the better diagnosis for this shape: the user's actual problem is not "your pair is heterogeneous", it is "you have BAM input on a FASTQ-output path". The format mismatch is secondary.

**Fix (smallest correct change):** derive `fmt_flag` as `" --output-format ubam"` when `matches!(cli.output_format, OutputFormat::UBam) || cli.clump_only` — i.e. suggest the flag when uBAM output *is* requested (so the remediation reproduces the user's mode) **or** when `--clump_only` is set and it is not (so the remediation adds what is missing). Note this makes `fmt_flag`'s semantics "the flag the *remediation* needs", not "the flag the user passed" — worth saying so in the doc comment, because the two coincide everywhere except here.

**Alternative fix:** exclude `cli.clump_only && output_format == Fastq` from the hoisted guard and let `clump_only.rs:399-410` keep its better message. Costs the early-firing benefit on that one path; I prefer the `fmt_flag` fix.

### C2 (Critical) — nothing in the plan's validation, or in the test suite, would catch C1

**Confirmed.** §9's V5 exercises only `--clump_only --paired --output-format ubam`. The FASTQ-output clump-only paired path — the one that regresses — appears nowhere in §9, and §11's edge-case list mentions `main.rs:426` only as "unaffected" (it *is* unaffected; the affected code is `clump_only.rs:399-410`, one call frame further in, which the plan never mentions).

The test suite does not cover it either. `tests/integration_clump_only.rs:303 rejects_ubam_input_without_ubam_output` asserts `stderr.contains("--output-format ubam")` and would be the natural guard — but I read it and it is **single-end** (no `--paired`), so it exercises `clump_only.rs:265`, not `:402`. The paired variant has zero coverage.

**Action:** add to §9 and to `tests/integration_clump_only.rs`:
- `--clump_only --paired <fq.gz> <bam>` (no `--output-format ubam`) → non-zero, stderr contains `--output-format ubam`;
- `--clump_only --paired <bam1> <bam2>` (no `--output-format ubam`) → same.

Both assertions survive the C1 fix and fail without it. As a bonus they close the §2.4 gap for the clump-only leg the same way §5's new test closes it for `:1268`.

### I1 (Important) — `{fmt1}`/`{fmt2}` is unspecified, and the obvious implementation leaks internal enum names

**Confirmed by reading.** §3.3's mixed message interpolates `{fmt1}`/`{fmt2}` and V1 requires the message to name "both files with their formats", but `InputFormat` (`format.rs:24-34`) derives only `Debug` — no `Display`, no label method. The path of least resistance is `{:?}`, which puts **`FastqGz`** and **`UnalignedBam`** into user-facing error text on a widely-used tool.

A suitable label function already exists — but it is **private to `clump_only.rs`**:

```rust
// src/clump_only.rs:685
fn input_format_label(fmt: InputFormat) -> &'static str {
    match fmt {
        InputFormat::FastqPlain => "FASTQ (plain)",
        InputFormat::FastqGz    => "FASTQ (gzip)",
        InputFormat::UnalignedBam => "uBAM",
    }
}
```

**Action:** the plan should specify the rendering. Promote `input_format_label` to `format.rs` (or `impl std::fmt::Display for InputFormat`) and reuse it, so the message reads *"… is FASTQ (gzip) and … is uBAM"*. Pairs with I5: if the helper lives in `format.rs`, the label lives next to it.

One sub-decision to make explicitly: whether the mixed message should distinguish plain from gzipped FASTQ. It must not imply that a plain+gzip pair is an error — it is accepted today. Either collapse both to `"FASTQ"` in this message, or keep the distinction but ensure the sentence attributes the rejection to BAM-vs-FASTQ only.

### I2 (Important) — §7 understates the multi-pair behaviour change

**Confirmed live.** §7 lists "no output directory created on rejection" as the notable side-effect change. The larger one is missing. Today, with the offence in a later pair:

```console
$ trim_galore --paired R1.fq.gz R2.fq.gz R2.fq.gz bam1.bam -o o10b
=== Pair 1 of 2 ===
… Pairs analyzed: 10000  Pairs removed: 4 …      # pair 1 fully trimmed AND WRITTEN
=== Pair 2 of 2 ===
Error: processing pair 2 of 2 … --paired with two BAM files is not supported…
```

Pair 1's `_val_1/_val_2` files and its trimming reports are on disk. After the hoist, **no pair is processed and nothing is written**. That is the right behaviour, but it is a bigger change than a missing directory: a pipeline that consumed whatever pairs succeeded before the failure now gets nothing. Step 8's CHANGELOG entry and §7's bullet list should say "multi-pair runs now fail before any pair is processed, so no partial output is written" rather than only mentioning the output directory.

### I3 (Important) — V8's command does not test pair-indexing on the current binary

**Confirmed live.** V8 proposes `--paired <fq> <fq> <fq> <bam>` (N=4). With the natural fixture choice the third argument repeats the first, and that trips the **output-collision pre-flight** at `main.rs:718-728` first:

```
Error: Output path collision (case-insensitive, for APFS/NTFS safety):
 o10/BS-seq_10K_R1_val_1.fq.gz and o10/BS-seq_10K_R1_val_1.fq.gz …
```

Post-fix the hoisted guard at `:264` fires *before* that pre-flight, so V8 would flip from "collision" to "pair 2 of 2 … mixed" and be recorded as a pass — while actually demonstrating guard-vs-preflight ordering, not pair indexing. Two consequences:

1. **Use four distinct paths.** `--paired R1.fq.gz R2.fq.gz R2.fq.gz bam1.bam` works: I verified it reaches the pair-2 error today (`processing pair 2 of 2 (R1=…R2.fastq.gz, R2=…bam1.bam)`), so the before/after comparison isolates the message and the index.
2. **State the new precedence.** The guard now runs ahead of the paired collision pre-flight (`:683`), ahead of `preflight_collision_bam` on the uBAM-out path (`:1763`), and ahead of the clump-only pre-flight (`:559`). For an input that is *both* collision-prone and format-mixed, the reported error changes. §7 should record it; it is a defensible ordering (format is the more fundamental error) but it is a change.

### I4 (Important) — Step 3 replaces a runtime invariant with a comment

**Confirmed by reading; the risk is inferred, not observed.** `main.rs:1308-1312`:

```rust
// Sequential path (--cores 1) — paired-BAM rejected above, so both
// are FastqReader.
let mut reader_r1 = FastqReader::open(input_r1)?;
```

"Rejected above" is the guard at `:1268` that Step 3 deletes. `--cores` defaults to `1` (`cli.rs:329`), so this is the **default** path. After the change the only thing standing between a BAM and `FastqReader::open` is a guard ~1000 lines away in `main()`, documented by a comment. A BGZF BAM handed to `FastqReader` decompresses to binary and is parsed as FASTQ — exactly the "fail silently and produce wrong results" class `CLAUDE.md` calls out.

The project's own convention is to keep the callee-side check. `clump_only.rs:920-928`:

```rust
// Belt-and-braces guard: main.rs::dispatch should have caught this;
// return a clear error if a caller ever slips a non-BAM N=1 through.
```

**Action:** at `:1268`, keep a cheap `bail!` (or `debug_assert`) with *internal-invariant* wording — e.g. "internal error: paired-BAM reached run_paired; the format guard in main() should have rejected this" — rather than only a pointer comment. Distinct wording keeps §11's `grep -n 'two BAM files' src/` completion check returning exactly one hit. The alternative (thread `&[InputFormat]` into `run_paired`) is cleaner but a wider diff; either is fine, a bare comment is not.

### I5 (Important) — Q2 is answerable now: `main.rs` has no test module

**Confirmed.** `grep -n 'cfg(test)' src/main.rs` → no match. `grep -c '#[test]'` → `src/main.rs:0`, `src/lib.rs:0`; all 421 tests live in library modules or `tests/`. A `#[cfg(test)] mod tests` inside a `[[bin]]` target *would* be compiled and run by `cargo test`, so the plan's stated blocker ("a `[[bin]]`-only helper cannot be unit-tested from `tests/`") is technically about integration tests only — but the crate has **zero** bin-target unit tests, so putting the first one in `main.rs` is off-convention.

**Action:** decide in the plan rather than at implementation time — put `reject_mismatched_pair_formats` in `src/format.rs`. It owns `InputFormat`, already has 7 unit tests, and would then also host the label function from I1. `anyhow` is a normal dependency (`Cargo.toml:32`), so the `Result` signature ports unchanged. Then Step 5's four unit tests (all-FASTQ, mixed both orders, two-BAM, multi-pair index) follow the existing pattern with no fixtures.

---

## 3. Assumptions

| # | Verdict | Note |
|---|---|---|
| A1 | Sound | `format.rs:40-77` read |
| A2 | Sound, and understated | Also holds for `--clump_only --paired` because `cli.rs:596` gates on plain `self.paired`. Worth adding — the guard covering clump-only depends on it |
| A3 | Sound | `:257-263` read; N=1 interleaved uBAM with FASTQ output verified working live |
| A4 | Sound reasoning, one factual sharpening | The claim "mixed pairs are harmless in hardtrim, which processes each file independently" is right for a better reason than stated: `--hardtrim5/3` never consults `cli.paired` at all (`:321-345`, `:346-370`), so `--paired` is inert there. The `--clock`/`--implicon` half of the exemption is different in kind — those *do* pair, and today they fail with a misleading message (`--clock`, verified live) or would proceed on equal read counts. The plan is right not to fix that here, but A4's justification only really covers hardtrim; the clock/implicon exemption is "don't change behaviour in this PR", which is fine but should be said |
| A5 | Sound by inspection | Substrings verified against all three assertion sites; still re-verify by running, as the plan says |
| A6 | Conclusion right, technical claim overstated | "Technically both `FastqReader` and `BamReader` implement `RecordSource`, so a mixed pair *could* be processed" is only true on the `--cores > 1` branch (`:1283-1284` `open_threaded_reader`). The default serial branch (`:1311-1312`) hard-codes `FastqReader::open`. Reword so a future reader doesn't take A6 as "it nearly works already" — that misreading is how I4 becomes a bug |

**Unstated assumption worth adding:** the guard reads `input_formats`, computed once at `:182`. Deleting the callee-side re-detections (`:1269`, `:1752`, `:537-538`) makes the *single* classification at `:182` authoritative for both validation and reader dispatch. That is an improvement — the guard and the reader can no longer disagree about a file — and it is the real argument for the hoist, stronger than the line-count argument in §6.

**No `-D warnings` risk from the deletions** (I checked, since CI is `-D warnings`): after removing the three `detect_input_format` uses, the `format::{InputFormat, detect_input_format}` import at `main.rs:16` is still used at `:28`, `:49`, `:66`, `:185`, `:257`, `:502`, `:1851`, `:1995`. No unused-import failure.

---

## 4. Efficiency

Nothing to add against the plan's §6 — the analysis is correct and the change is a strict win. Two amplifications:

- The deletions remove real I/O, not just lines: `:1269`, `:1752` and `:537-538` each re-open the file and decompress a BGZF block (`format.rs:64-72`) for inputs already classified at `:182`. That is 2 redundant opens + decompresses per pair on the rejected paths, and on the two-file uBAM-out path it happens before every pair.
- Rejected runs skip the up-to-1M-read adapter scan (`adapter.rs`) and, on multi-pair input, skip *processing earlier pairs entirely* (I2) — which is where the real saving is.

`Vec<InputFormat>` is `Copy` per-element and N is a handful; no allocation concerns.

---

## 5. Validation sufficiency

The table is well constructed — V3 (no side effects), V9 (tests unedited), V11 (validation matrix) are exactly the three that would quietly rot, and the plan is right to fence them. Gaps:

| Gap | Severity | Fix |
|---|---|---|
| No clump-only **FASTQ-output** paired case (the C1 regression path) | Critical | Add the two invocations from C2, asserting stderr still contains `--output-format ubam` |
| V8's command tests the wrong thing (I3) | Important | Use four distinct paths; `R1.fq.gz R2.fq.gz R2.fq.gz bam1.bam` verified to reach pair 2 today |
| No check that multi-pair rejection writes **no partial output** (I2) | Important | Extend V3: run the V8 command with `-o freshdir`, assert the directory is absent — i.e. that pair 1's `_val_*` files were never written |
| V6/V11 have no mechanism | Optional | Before touching code, md5 `BS-seq_10K_R{1,2}_val_{1,2}.fq.gz` from a PE run and the SE/hardtrim5/clock outputs; re-compare after. Turns "confirm by inspection" into a check. The *reasoning* behind V11 is sound (the guard only ever rejects), so this is belt-and-braces — but V11 is on the plan's own do-not-skip list |
| No assertion on what the message must **not** contain | Optional | V1 says "**no** bare interleaved-file imperative" in prose; make it an assertion in the new integration test (`!stderr.contains("expects a single interleaved file")` for the mixed case). Otherwise the one substring whose *absence* is the point of #363 is unpinned, and a later wording merge could reintroduce it silently |
| Rendering of `{fmt1}`/`{fmt2}` unpinned (I1) | Optional | Assert `stderr.contains("uBAM")` and `!stderr.contains("UnalignedBam")` |

Two things I checked that are **not** problems: `CARGO_BIN_EXE_trim_galore` (`integration_clump_only_ubam.rs:16`) means `cargo test` always builds the binary it exercises, so V9/V10 cannot pass against a stale build; and the §9 preamble's `cargo build --release` note correctly covers the *manual* V1–V8 runs against `./target/release/trim_galore`.

---

## 6. Alternatives

1. **Message-only, guards left in place.** Already rejected with the user, and correctly — it keeps the output directory, the adapter scan, and (I2) the partial multi-pair output.
2. **Guard inside `Cli::validate()`.** Impossible without file I/O in `validate()`; the `--phred64` precedent documents exactly this at `main.rs:206-210`. The plan's siting is right.
3. **Accept and process mixed pairs.** Dismissed by A6; I agree. Note it is further from working than A6 implies (see A6 above).
4. **Shared helper in `format.rs` + belt-and-braces asserts at the callees.** This is my recommended shape: it resolves I5 (unit-testable, on-convention), I1 (label lives beside `InputFormat`), and I4 (callees keep an enforced invariant, not a comment) in one move, at the cost of one new `pub fn` in `format.rs`.
5. **Keep clump-only FASTQ output out of the guard** (C1's alternative fix). Simpler than reworking `fmt_flag`, but loses the early-firing benefit on that path and leaves the two-BAM clump-only-FASTQ case reporting via a per-file message that doesn't mention the pair. I'd take the `fmt_flag` fix.

---

## 7. Action items

### Critical

1. **Fix the clump-only remediation.** Derive `fmt_flag` as `" --output-format ubam"` when `matches!(cli.output_format, UBam) || cli.clump_only`, and document that `fmt_flag` means "the flag the remediation needs", not "the flag the user passed". Without this, `--clump_only --paired <bam> <bam>` (FASTQ output) prints a command that the binary rejects — verified live. (§3.3, §5; `main.rs:426`, `clump_only.rs:399-410`)
2. **Add the two missing validation cases** that make C1 visible: `--clump_only --paired <fq> <bam>` and `--clump_only --paired <bam1> <bam2>`, both without `--output-format ubam`, asserting stderr contains `--output-format ubam`. `tests/integration_clump_only.rs:303` is single-end and does not cover these. (§9, §2.4)

### Important

3. **Specify `{fmt1}`/`{fmt2}` rendering.** Promote `clump_only.rs:685 input_format_label` to `format.rs` (or add `Display`) and reuse it; `{:?}` would print `FastqGz`/`UnalignedBam` to users. Decide whether the message collapses plain/gzip FASTQ. (§3.3)
4. **Resolve Q2 in the plan:** put the helper in `src/format.rs`. `main.rs` has zero `#[test]`s and no `cfg(test)` module — verified — so a first bin-target unit test would be off-convention. (§4 Step 5, §10 Q2)
5. **Keep an enforced invariant at `:1268`,** not just a pointer comment: the default `--cores 1` path at `:1311` hard-codes `FastqReader::open` on the strength of the guard being deleted. Use internal-invariant wording so §11's grep check still returns one hit. (§4 Step 3)
6. **Record the multi-pair change** in §7 and in Step 8's CHANGELOG entry: previously-written partial output from earlier pairs is no longer produced. Verified live. (§7, §4 Step 8)
7. **Fix V8's command** to use four distinct paths (`R1 R2 R2 bam1` verified to reach pair 2 today), and note in §7 that the guard now precedes all three output-collision pre-flights (`:683`, `:1763`, `:559`), so a doubly-invalid input reports a different error than before. (§9 V8, §7)
8. **Extend V3 to the multi-pair case** — assert no `_val_*` files from pair 1, not just an absent output directory.

### Optional

9. **§2.3 wording:** "Nothing has been created or read when the guard runs" is false — `sanity_check_any(&cli.input[0])` (`:176`) and `detect_input_format` over all inputs (`:182-186`) both read from disk first. §7's "after `sanity_check_any`" is right. Only the "created" half is load-bearing; drop "or read".
10. **Line-number drift** (cosmetic, but the plan's precision is otherwise a strength): §2.2's `:543`/`:552` → the `bail!`s are at `:542`/`:551`, block spans `:536-558` (§4 Step 4 says `:537-557`); §2.3's "specialty dispatch at ~:300" → `:321`; §4 Step 4's `clump_only.rs:940-944` → the comment is at `:942-945`.
11. **Step 8:** a `#### Fixes` block already exists under `### Unreleased` — append to it rather than adding a second.
12. **Sharpen A6** so nobody later reads it as "mixed pairs nearly work": only the `--cores > 1` branch is format-polymorphic; the default serial branch is `FastqReader`-only.
13. **Sharpen A4:** its "processes each file independently" justification covers hardtrim (which ignores `cli.paired` entirely) but not `--clock`/`--implicon`, which do pair and fail misleadingly today. Their exemption is "out of scope for this PR", which is fine — say that.
14. **Add the negative assertion** to the new mixed-pair test (`!stderr.contains("expects a single interleaved file")`). The absent substring is the substance of #363 and is currently unpinned.
15. **Note the `--passthrough` precedence** in §7: `:251-256` rejects passthrough against any BAM before `:264`, so that invocation keeps its own message. Correct as-is; documenting it prevents a future "inconsistency" fix.
16. **§6 could claim more:** the deletions remove redundant `detect_input_format` calls (file open + BGZF block decompress ×2 per pair), and make the single classification at `:182` authoritative for both validation and reader dispatch — a correctness argument, not just a line-count one.
