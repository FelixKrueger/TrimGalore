# PLAN v2 — Accurate per-pair format rejection (#363)

**Issue:** [#363](https://github.com/FelixKrueger/TrimGalore/issues/363) — *Mixed FASTQ+BAM pair under `--paired` reports "two BAM files" when only one input is a BAM*
**Branch:** `fix/mixed-format-pair-message` off `dev`
**Reviews incorporated:** `PLAN_review_reviewer-A.md`, `PLAN_review_reviewer-B.md` (dual independent, both verified against the code and by running the binary)

**Scope decisions taken with the user:**

- Message **and** guard placement, not message-only.
- Specialty modes (`--hardtrim5/3`, `--clock`, `--implicon`) **excluded**; follow-up issue instead.
- **Not** release-gated — the v2.4.0 cut is parked pending feedback on the usefulness of `--clump_only`.
- Reviewer A's **mode enum** adopted over Reviewer B's minimal `fmt_flag` fix (§10 R1).
- Already-dead code at `main.rs:498-512` **retired** in the same pass (§4 Step 7).

See §12 for the revision history and what v1 got wrong.

---

## 1. Goal

Two user-visible outcomes:

1. **Correct diagnosis per input shape.** A pair containing exactly one BAM must be told that both inputs of a pair have to share a format, naming which input is which, and must *not* be offered the single-interleaved-file remediation. A pair of two BAMs keeps that remediation, because for it the advice is right. `--clump_only` on the FASTQ-output path must keep its own, better diagnosis (the problem there is BAM-on-a-FASTQ-path, not heterogeneity).
2. **Reject before doing work.** The rejection must fire in the cheap-validation window — before adapter auto-detection, before the trimming banner, before the output directory is created, and before any pair is processed.

Non-goal: changing *whether* mixed pairs are rejected. Rejection is right. Only the diagnosis and its timing change.

---

## 2. Context

### 2.1 The two defective sites

Both loop over the pair's inputs and `bail!` on the **first** BAM encountered, so the two-BAM message is emitted for any pair containing at least one BAM:

| Site | Path | Reached by |
|---|---|---|
| `src/main.rs:1268-1278` | `run_paired` — trim, FASTQ output | `--paired <fastq> <bam>` |
| `src/main.rs:1750-1762` | `run_ubam_output` — paired Shape A | `--paired --output-format ubam <fastq> <bam>` |

Reproduced on `dev` @ `a4ffd47`, both argument orders. On the trim path the error arrives only after adapter auto-detection has run, poly-G detection has run, the banner has printed planned output paths, and the output directory has been created. Both reviewers reproduced this independently.

On the uBAM-output path (`:1750`) the error fires earlier — only `ensure_output_dir` and the uBAM startup NOTEs precede it. v1 generalised one site's stderr to both; corrected here.

### 2.2 The issue mis-attributes a third site

`#363` lists `src/main.rs:543` as a cause site. It is not — it is the **reference implementation**. The block at `main.rs:536-558` already branches correctly, with `bail!` calls at `:542` (two-BAM) and `:551` (mixed):

```rust
if r1_is_bam && r2_is_bam { anyhow::bail!("--clump_only --paired with two BAM files …") }   // :542
if r1_is_bam != r2_is_bam { anyhow::bail!("--clump_only --paired requires both input …") }  // :551
```

Confirmed live by both reviewers and by me. It also has regression coverage (§2.5). **The fix therefore touches two defective sites, not three**, and this one becomes redundant once the guard is hoisted (§4 Step 5). The issue will be corrected with a comment rather than silently.

### 2.3 Where the new guard goes, and why there

Target: `src/main.rs`, immediately after the N=1 uBAM check at `:257-263`, before the `gzip`/`output_dir` block at `:265-277`.

Four properties of that position, each load-bearing:

- **`input_formats` already exists.** `:182-186` detects every input's format once, explicitly so validation and dispatch can share it. The guard needs no new I/O — and removes some, since two of the three replaced guards each re-called `detect_input_format` (§6).
- **It is the established home for format-dependent validation.** The `--phred64` guard at `:211` sits here for the reason documented at `:206-210`: *"validate() runs before input format detection and cannot see `any_bam`."* Same constraint.
- **Preconditions are already proven.** By `:264`, `Cli::validate()` has enforced an even input count (`cli.rs:506`) and R1≠R2 per pair (`cli.rs:513`), and `main.rs:257` has rejected `--paired` + N=1 + non-BAM. So the guard sees either N=1 (a proven uBAM, correctly exempt) or an even N≥2 with distinct members.
- **It precedes every write.** `naming::ensure_output_dir` is at `:277`; dispatch branches are all below it — hardtrim5 `:321`, clump-only `:409`, uBAM-out `:653`, paired trim `:675`. Nothing is *created* when the guard runs.

Note (corrected from v1): the guard is **not** the first thing to touch the disk. `sanity_check_any(&cli.input[0])` reads at `:176` and `detect_input_format` reads every input at `:182-186`. Only the *created-nothing* half of the claim is load-bearing. A mixed pair whose R1 is an empty BAM therefore still reports the empty-BAM error, which is correct precedence.

### 2.4 New precedence the hoist establishes

The guard now runs ahead of all three output-collision pre-flights — paired trim (`:683`), clump-only (`:559`), uBAM-out (`:1763`). For an input that is *both* collision-prone and format-mixed, the reported error changes from the collision to the format. That is the right precedence (format is the more fundamental error) but it is a change, and it is why V8's command had to be rewritten (§9).

`--passthrough` keeps its own precedence: `main.rs:251-256` rejects it against *any* BAM input, i.e. **before** `:264`. So `--paired --passthrough pt.fq <fq> <bam>` keeps its more specific message. Correct as-is; documented so a future "consistency" fix does not undo it.

### 2.5 Existing test coverage, and two gaps

| Test | Asserts stderr contains | Exercises |
|---|---|---|
| `integration_clump_only_ubam.rs:363` `rejects_two_bam_paired` | `"single interleaved"` ∨ `"uBAM paired mode expects"` | clump-only uBAM, two distinct BAMs |
| `integration_clump_only_ubam.rs:594` `rejects_mixed_format_paired` | `"same format"` ∨ `"mixed"` | clump-only uBAM, mixed pair |
| `integration_ubam_out.rs:443` `ubam_out_two_bam_pair_rejected` | `"two BAM files is not supported"` ∨ `"interleaved file"` | site `:1750`, two distinct BAMs |

All three assert on loose substrings, which is what makes a consolidated guard viable: §3.4's message set retains every one, so **no existing test needs editing**. If any requires a change, the wording drifted too far — treat it as a review finding, not a test to update. §3.4 additionally reflows the text so no pinned substring straddles a Rust line-continuation (§3.5).

**Gap 1:** nothing covers site `:1268` — plain `--paired <bam> <bam>` without `--output-format ubam`.
**Gap 2 (the one that made v1 wrong):** nothing covers `--clump_only --paired` on the **FASTQ-output** path. `integration_clump_only.rs:303 rejects_ubam_input_without_ubam_output` looks like the guard but is **single-end**, so it exercises `clump_only.rs:265`, not `:403`.

### 2.6 Explicitly out of scope

Both to be filed as a separate issue, not fixed here:

- `--hardtrim5 10 --paired <fastq> <bam>` **silently accepts** the mixed pair and writes output. Mechanism confirmed: `main.rs:321-345` loops `for input in &cli.input` and never consults `cli.paired`, so pairing is genuinely inert there.
- `--clock --paired <fastq> <bam>` fails with `Paired-end files have different numbers of reads!`, misleading in the same way, and would have *proceeded* had the record counts matched.

The guard must therefore **not** fire for the four specialty modes; those invocations work today and this change must not remove them. See §8 A4 for the sharpened justification.

---

## 3. Behavior

### 3.1 Guard activation

Runs when **all** hold:

- `cli.paired`
- `cli.input.len() > 1` — N=1 is single-file interleaved uBAM, already validated at `:257-263`
- no specialty mode active: `cli.hardtrim5.is_none() && cli.hardtrim3.is_none() && !cli.clock && cli.implicon.is_none()`

`cli.clump_only` is **not** excluded; it is handled by its own shape (§3.2). Field types verified: `hardtrim5`/`hardtrim3`/`implicon` are `Option<usize>` (`cli.rs:372`, `:377`, `:391`), `clock`/`clump_only`/`paired` are `bool` (`:383`, `:208`, `:89`). `cli.rs:805-814` additionally makes `--clump_only` mutually exclusive with all four specialty modes, so the exemption cannot accidentally swallow clump-only.

### 3.2 Shape determination and per-pair decision

The rejection differs by mode, so the helper takes a shape enum resolved once by the caller:

| Shape | Condition |
|---|---|
| `ClumpOnlyUbamOut` | `cli.clump_only && output_format == UBam` |
| `ClumpOnlyFastqOut` | `cli.clump_only` (FASTQ output) |
| `TrimUbamOut` | `output_format == UBam` |
| `Trim` | otherwise |

For each `chunk` of 2 over the paired inputs, with `n_bam` = count of `InputFormat::UnalignedBam`:

**`ClumpOnlyFastqOut` — predicate is `n_bam >= 1`, not mismatch.** The user's problem in this mode is BAM input on a FASTQ-output path, whether one input or both is a BAM. Emits the aux-tag message (§3.4), preserving today's behaviour exactly. This is the case v1 got wrong.

**All other shapes:**

| `n_bam` | Action |
|---|---|
| 0 | continue — all-FASTQ pair |
| 1 | mixed-format message |
| 2 | two-BAM message |

First offending pair wins; iteration stops there. The message names the 1-based pair index **when there is more than one pair** — "Pair 1 of 1" is noise on the common case.

### 3.3 Critical: the predicate is BAM-count, not format equality

`InputFormat` has **three** variants (`format.rs:25-34`): `FastqPlain`, `FastqGz`, `UnalignedBam`. A naive `formats[0] != formats[1]` would reject a plain+gzipped FASTQ pair, which works today — verified: `--paired plainR1.fastq R2.fastq.gz` exits 0 and writes `plainR1_val_1.fq` + `BS-seq_10K_R2_val_2.fq`.

The helper is therefore named for the BAM-count predicate, not for format equality, and §9 V12 pins the plain+gz acceptance.

### 3.4 Message text

Format names are rendered by a promoted `input_format_label` (§4 Step 2) — `InputFormat` has no `Display`, so `{:?}` would leak `FastqGz`/`UnalignedBam` into user-facing text. Substrings pinned by existing tests are marked ✚.

**Two BAM files** (`Trim`, `TrimUbamOut`, `ClumpOnlyUbamOut`) — semantics unchanged, plus both paths and a *verified* remediation:

```
{mode} with ✚two BAM files is not supported.
✚uBAM paired mode expects a ✚single interleaved file:
`trim_galore {mode}{fmt_flag} interleaved.bam`.
{Pair i of n is }two BAM files: {r1} and {r2}.
Combine them into one mate-adjacent file first:
`samtools merge -n -o interleaved.bam {r1} {r2}`.
```

Retains ✚`two BAM files is not supported`, ✚`single interleaved`, ✚`uBAM paired mode expects`. Names **both** files, replacing v1's *"one of them is {R1}"* — naming a single file when both are BAM was never informative.

**The `samtools` command is verified, not inferred.** v1 proposed `samtools collate -O r1.bam r2.bam`, copied from `bam.rs:52-58`'s `GROUPED_INPUT_ERR`. Tested on samtools 1.21 with two 10 000-record uBAMs: it exits **0** and emits **10 000** records — `collate`'s second positional is the temp-file *prefix*, so `r2.bam` is silently consumed as a prefix and half the library is dropped with no warning. `samtools merge -n -o interleaved.bam r1.bam r2.bam` produces all 20 000 records, mate-adjacent (flags 69/133 alternating), and TrimGalore accepts the result. Reviewer A's `merge -n | collate` pipeline also works; `merge -n` alone is simpler and sufficient, because `merge -n` presumes name-sorted inputs and any FASTQ-derived uBAM is name-sorted. If a future report shows non-name-sorted inputs in the wild, `| samtools collate -O - tmp` is the verified fallback.

**Mixed format** (`Trim`, `TrimUbamOut`, `ClumpOnlyUbamOut`) — the new case:

```
{mode} requires both inputs of a pair to be the ✚same format.
{Pair i of n is }✚mixed: {r1} is {fmt1} and {r2} is {fmt2}.
Pass two FASTQ files, or a single interleaved uBAM.
If you meant two FASTQ files, check for a mis-typed filename.
```

Retains ✚`same format` and ✚`mixed`. Deliberately **no** bare interleaved-file imperative — for a mis-typed filename that advice misdirects, which is the substance of #363. §9 V14 pins its absence.

The sentence attributes the rejection to BAM-vs-FASTQ only, so a plain-vs-gzip pair can never read as an error (§3.3).

**`ClumpOnlyFastqOut`** — preserve today's text from `clump_only.rs:403-409` verbatim:

```
uBAM input under --clump_only requires --output-format ubam
(using the FASTQ output path with uBAM input would drop aux tags).
Input: {first_bam}
```

Verbatim preservation is deliberate: it means this mode's observable behaviour is unchanged except for firing earlier, so the new tests in §9 V13 pin today's message and would fail under v1's design.

Wording constraint from a shipped feature: mixed formats *across a batch* are deliberately allowed (`integration_ubam_out.rs:340`, single-end). Hence "both inputs of **a pair**" throughout. Register per `feedback_docs_register`: plain statements, no filler, no second-person scolding.

### 3.5 Rust line-continuation hazard

`\` + newline in a Rust string literal consumes the newline **and all leading whitespace on the next line**. A break between `"a single"` and `"interleaved file"` therefore yields `"a singleinterleaved file"`, silently breaking `rejects_two_bam_paired`'s pinned substring. Existing code avoids this by ending each line with a space before the `\` (`main.rs:1271-1274`).

**Implementation rule:** none of the four pinned substrings — `two BAM files is not supported`, `single interleaved`, `uBAM paired mode expects`, `same format` — may straddle a source-line break. V9 is the backstop, not the primary defence.

### 3.6 Side effects

None. No reader opened, no adapter detection, no trimming banner, no output directory, no pair processed. The version banner at `main.rs:172-174` prints unconditionally on every run and is unaffected. Exit is non-zero via `anyhow::bail!`.

---

## 4. Implementation outline

1. **Add the shape enum and helper to `src/format.rs`** (not `main.rs` — §10 Q2 resolved). `format.rs` owns `InputFormat`, already has a unit-test module including the load-bearing bgzip-vs-BAM case, and a pure decision function over `&[InputFormat]` is its remit. `main.rs` has **zero** `#[cfg(test)]` modules (verified), so a first bin-target unit test would be off-convention. `anyhow` is a normal dependency (`Cargo.toml:32`), so the `Result` signature ports unchanged.

2. **Promote `input_format_label`** from `clump_only.rs:685` (currently private) to `format.rs`, or add `impl Display for InputFormat`. Update `clump_only.rs`'s use site. This is what keeps enum variant names out of user-facing text.

3. **Call the helper at `src/main.rs:264`**, after the N=1 check, before the `gzip` binding at `:272`. Gate per §3.1, resolve the shape per §3.2. Comment it in the register of the `--phred64` guard above: why here (formats known, nothing created yet), and why specialty modes are exempt.

4. **Replace the defective site at `:1268-1278` with an enforced internal-invariant check**, not a bare comment. The default `--cores 1` path at `:1311` hard-codes `FastqReader::open(input_r1)`, and a BGZF BAM handed to `FastqReader` decompresses to binary and parses as FASTQ — silent wrong output, not an error. Use internal-invariant wording (e.g. *"internal error: uBAM input reached run_paired; the format guard in main() should have rejected it"* — note the enclosing function is `main.rs::run_paired`, **not** `trimmer::run_paired_end`) so §11's `grep -n 'two BAM files' src/` completion check still returns exactly one hit.

5. **Delete the defective site at `:1750-1762` and the now-redundant clump-only branch at `:536-558`**, each replaced by a pointer to the hoisted guard. Keep a cheap belt-and-braces check at the two leaf functions that read the source header from **R1 only** and then open each side independently — `run_ubam_output_paired_two_files` (`main.rs:1993-1995`) and `clump_only_paired_to_bam_one_pair` Shape A (`clump_only.rs:946-958`). A mixed pair reaching either would emit a BAM mixing FASTQ-derived records (no aux tags, no source header) with BAM-derived ones, silently. The project's own convention supports this: `clump_only.rs:920-928` already carries a *"Belt-and-braces guard: main.rs::dispatch should have caught this"* check.

6. **Update all four stale comments**, not the two v1 listed:
   - `main.rs:1266-1267` — *"handled in main() before this fn is called"* (aspirational in v1; the hoist makes it true — state where).
   - `main.rs:1309-1310` — *"paired-BAM rejected above, so both are FastqReader"* — "above" ceases to exist.
   - `main.rs:1993-1994` — *"At this point both inputs are FASTQ"* — an invariant now maintained ~1 700 lines away.
   - `clump_only.rs:942-945` — cites *"Format-guards in main.rs::dispatch"*; the location changes.

7. **Retire the dead check at `main.rs:498-512`** (user-approved). Its condition — `clump_only ∧ UBam ∧ paired ∧ N=1 ∧ ¬BAM` — is a strict subset of `main.rs:257-263`'s, so it is already unreachable. Same family of redundancy this change retires, one branch from code already being edited.

8. **Unit tests in `format.rs`'s test module** (no fixtures needed): all-FASTQ pair passes; plain+gz pair passes (§3.3); mixed rejected in both orders with the mixed message; two-BAM rejected with the two-BAM message; `ClumpOnlyFastqOut` rejects on `n_bam == 1` *and* `n_bam == 2` with the aux-tag message; multi-pair offence reports the right index; single-pair omits the index.

9. **Integration tests** per §9, covering both §2.5 gaps.

10. **Verify the three existing tests pass unedited.** Any failure means the wording drifted — fix the wording, not the test.

11. **CHANGELOG.md** — append to the **existing** `#### Fixes` block under `### Unreleased` (`CHANGELOG.md:84`); do not add a second. Cite #363. Cover: the wording change, the earlier firing, that no output directory is created on rejection, and that **multi-pair runs no longer write partial output for pairs preceding the offending one** (§7).

12. **File the follow-up issue** for §2.6 (specialty modes on mixed pairs), cross-referencing #363.

---

## 5. Signature

```rust
/// Which paired-input shape is being validated. Determines both the predicate
/// and the remediation, which differ by mode: `ClumpOnlyFastqOut` rejects any
/// BAM in the pair (its problem is BAM-on-a-FASTQ-path, not heterogeneity),
/// while the others reject only a format mismatch or a two-BAM pair.
pub enum PairedShape {
    /// `--paired`, FASTQ output.
    Trim,
    /// `--paired --output-format ubam`.
    TrimUbamOut,
    /// `--clump_only --paired`, FASTQ output. Must keep the aux-tag diagnosis
    /// from `clump_only.rs:403-409`: for this mode the fix is a flag, not a
    /// re-shaped input, and suggesting `--paired interleaved.bam` produces a
    /// command that `main.rs:426-432` rejects.
    ClumpOnlyFastqOut,
    /// `--clump_only --paired --output-format ubam`.
    ClumpOnlyUbamOut,
}

/// Reject paired inputs whose two members disagree on being BAM.
///
/// The predicate is a **BAM count per pair**, not format equality:
/// `InputFormat` has three variants, and a `FastqPlain` + `FastqGz` pair is
/// legal (verified). Only BAM-vs-FASTQ within one pair is an error.
///
/// Operates on the formats detected once in `main()` (`input_formats`), so it
/// performs no I/O and runs before any output directory is created — which is
/// the point: the sites this replaces fired only after adapter auto-detection
/// had scanned the inputs, the output directory existed, and (on multi-pair
/// input) earlier pairs had already been written.
///
/// Callers must skip specialty modes (`--hardtrim5/3`, `--clock`,
/// `--implicon`), which accept mixed pairs today, and must not call this for
/// `N == 1` (single-file interleaved uBAM, validated at `main.rs:257-263`).
///
/// Assumes `Cli::validate()` has already enforced an even input count
/// (`cli.rs:506`) and R1 != R2 per pair (`cli.rs:513`), and that
/// `inputs.len() == formats.len()` — the last is `debug_assert`ed, since the
/// `chunks(2)` indexing depends on it.
pub fn reject_bam_format_mismatch_in_pair(
    inputs: &[std::path::Path],
    formats: &[InputFormat],
    shape: PairedShape,
) -> anyhow::Result<()>
```

The shape enum replaces v1's two opaque `&str` parameters. That is not merely tidier: it is what makes the `ClumpOnlyFastqOut` remediation expressible at all, and it removes the risk of a caller passing a `mode` string that disagrees with the actual flags.

---

## 6. Efficiency

No measurable cost: one pass over an already-materialised `Vec<InputFormat>` of a handful of `Copy` elements, with early exit. O(N/2).

Strictly *less* work than today, in three ways — the third is a correctness argument, not a performance one:

1. **On rejected runs**, the whole adapter-detection scan (up to 1 M reads, `adapter.rs`) is skipped, along with poly-G detection and `ensure_output_dir`.
2. **On all runs**, two of the three replaced guards re-called `detect_input_format` (`:1752`, `:537-538`) — a file open plus a BGZF-block decompress, twice per pair — on inputs already classified at `:182`. The third (`:1269`) is **retained** as an internal-invariant backstop (§4 Step 4), so the ordinary paired FASTQ path still pays two detections per pair. That is the deliberate trade for enforcing the invariant rather than documenting it, and the cost is per pair, not per read.
3. **The single classification at `:182` becomes authoritative** for both validation and reader dispatch. The guard and the reader can no longer disagree about a file's format. This is the strongest argument for the hoist; v1 rested on line count.

On multi-pair input the largest saving is behavioural: earlier pairs are no longer processed before the failure (§7).

---

## 7. Integration

**Reads:** `cli.input`, `cli.paired`, the specialty flags, `cli.clump_only`, `cli.output_format`, and the pre-computed `input_formats`. **Writes:** nothing.

**Order:** after `Cli::validate()`, after `sanity_check_any` (`:176`), after format detection (`:182`), after the `--phred64` (`:211`), `--preserve-tags` (`:233`), `--passthrough` (`:251`) and N=1 (`:257`) guards — before `ensure_output_dir` (`:277`), before all three collision pre-flights (§2.4), and before every dispatch branch.

**Behaviour changes visible to users:**

1. Mixed pairs get an accurate message instead of a false one — the point of the fix.
2. **Multi-pair runs no longer write partial output.** Today the trim-path guard is *inside* the per-pair loop, so earlier pairs complete first. Verified on `dev`: `--paired R1.fq.gz R2.fq.gz R2.fq.gz ubam_test.bam` leaves pair 1's `BS-seq_10K_R{1,2}_val_{1,2}.fq.gz` **and both trimming reports** on disk before failing on pair 2. After the hoist that run produces nothing. Fail-fast is right — failing in a second beats failing after hours on pair 97 of 100 — but a pipeline that consumed whatever pairs succeeded now gets nothing. This is the largest change and belongs in the CHANGELOG.
3. Rejected paired runs no longer create the output directory or run adapter detection.
4. Rejection errors lose the `processing pair N of M` context wrapper, replaced by an explicit pair index inside the message (omitted when there is one pair).
5. For input that is both collision-prone and format-mixed, the reported error changes from collision to format (§2.4).

No accepted invocation changes. Specialty modes are untouched. `--passthrough` keeps its own earlier rejection.

**No `-D warnings` exposure from the deletions:** the `format::{InputFormat, detect_input_format}` import at `main.rs:16` remains used at `:28`, `:49`, `:66`, `:185`, `:257`, `:502`, `:1851`, `:1995`.

---

## 8. Assumptions

- **A1.** `detect_input_format` is authoritative and content-based (`format.rs:40-77`), so `bgzip`-framed FASTQ classifies as FASTQ, not BAM — guarded by `detect_bgzipped_fastq_is_fastq_not_bam`. The guard adds no format logic of its own. Fixed rule.
- **A2.** `Cli::validate()` guarantees an even input count for paired mode (`cli.rs:506`) with the N=1 carve-out (`cli.rs:503`), and R1≠R2 per pair (`cli.rs:513`). Reachability also verified: `validate_paired_input("Paired-end")` is invoked at `cli.rs:596` under plain `if self.paired`, so it covers `--clump_only --paired` too — which the guard's coverage of clump-only depends on.
- **A3.** `--paired` + N=1 is legal only for uBAM and is rejected otherwise at `main.rs:257-263`. Verified: `--paired <interleaved.bam>` with FASTQ output succeeds, so §3.4's "or a single interleaved uBAM" hint is valid on the trim path. It is **not** valid on the clump-only FASTQ path (`main.rs:426-432` rejects it), which is why `ClumpOnlyFastqOut` exists.
- **A4.** Excluding specialty modes is a user-confirmed departure from #362's *"rejected uniformly rather than per-mode"* reasoning, and the justification differs by mode:
  - `--hardtrim5/3` — genuinely harmless. `main.rs:321-370` never consults `cli.paired`, so pairing is inert and each file is processed independently. Confirmed live.
  - `--clock` / `--implicon` — **not** harmless. These do pair, and `--clock` today fails with a misleading read-count error and would have proceeded on equal counts. Their exemption is *"out of scope for this PR"*, not *"correct as-is"*. The follow-up issue (§4 Step 12) is doing real work and is not optional.
- **A5.** Existing loose `contains` assertions (§2.5) are satisfied by §3.4's wording — verified by substring inspection, and re-verified by running (V9), not by re-reading the strings. §3.5 is the specific way this can still fail.
- **A6.** Rejection remains correct for mixed pairs, and the accepting path is further from working than v1 implied. Both `FastqReader` and `BamReader` implement `RecordSource`, but only the `--cores > 1` branch is format-polymorphic (`:1283-1284` `open_threaded_reader`); the **default** serial branch hard-codes `FastqReader::open` (`:1311`). And both two-file leaf functions take the source header from R1 alone, so a mixed pair "processed" would mean a BAM with a header from one arbitrary side. Rejection is not merely the likelier-intent call — the accepting path does not exist in a correct form. This is also the argument for Step 5's belt-and-braces checks. Worded carefully so no future reader takes A6 as "it nearly works already".
- **A7.** `inputs.len() == formats.len()`. This is what the `chunks(2)` indexing actually relies on, and it was unstated in v1. `debug_assert_eq!` in the helper.

---

## 9. Validation

Run from the crate root; `cargo build --release` first — `cargo clippy` leaves the binary stale (handoff §5). Note `CARGO_BIN_EXE_trim_galore` (`integration_clump_only_ubam.rs:16`) means `cargo test` always builds the binary it exercises, so V9/V10 cannot pass against a stale build; the manual V1–V8 runs are what need the explicit build.

| # | Verify | How | Expected |
|---|---|---|---|
| V1 | Mixed pair, both orders, gets the mixed message | `--paired <fq> <bam>` and `--paired <bam> <fq>` | Non-zero; "both inputs of a pair"/"same format"; names both files with their formats |
| V2 | Two distinct BAMs keep the two-BAM message | `--paired <bam1> <bam2>` | Non-zero; "two BAM files is not supported" + interleaved remediation + the `samtools merge -n` hint |
| V3 | **No side effects on rejection** — the §1 outcome most easily lost in a refactor | Automated: pass `-o <fresh_tmpdir>/nested`, assert `!nested.exists()`. **Gotcha:** `fresh_tmpdir` (`integration_ubam_out.rs:28-32`) does `create_dir_all`, so passing it directly guarantees existence and the check would be vacuous | Directory absent; stderr shows no adapter-detection, poly-G, or trimming-banner lines |
| V4 | uBAM-output path fixed too | V1 and V2 with `--output-format ubam` | Same messages, mode-appropriate remediation |
| V5 | clump-only **uBAM-output** messages did not regress after deleting `:536-558` | `--clump_only --paired --output-format ubam` with mixed, then two BAMs | Mixed → mixed message; two-BAM → two-BAM message; remediation names `--clump_only --paired --output-format ubam` |
| V6 | Legal invocations unaffected | Paired FASTQ pair; `--paired one.bam`; multi-pair FASTQ (N=4, distinct paths) | All succeed, output unchanged |
| V7 | Specialty modes untouched | `--hardtrim5 10 --paired <fq> <bam>`; `--clock --paired <fq> <bam>` | Byte-identical to `dev` — hardtrim writes output and exits 0; clock still reports its read-count error |
| V8 | Multi-pair offence names the right pair | `--paired R1.fq.gz R2.fq.gz R2.fq.gz <bam>`. The property that matters is **distinct output paths**, not distinct inputs: v1's `<fq> <fq> <fq> <bam>` collided on output and tripped the collision pre-flight first, so post-fix it would have flipped to the guard and been recorded as a pass while actually testing precedence, not indexing. Reusing `R2.fq.gz` as pair 1's R2 and pair 2's R1 is fine — the four outputs differ and within-pair R1≠R2 holds | Rejected naming **pair 2 of 2** |
| V9 | Existing tests pass **unedited** | `cargo test rejects_two_bam_paired rejects_mixed_format_paired ubam_out_two_bam_pair_rejected` | 3 passed, zero diff in `tests/` for these three |
| V10 | Full suite + lint | `cargo test`; `cargo fmt --all -- --check`; `cargo clippy --all-targets --release -- -D warnings` | 421 + new tests pass; zero warnings |
| V11 | Validation matrix untouched | md5 the PE (`BS-seq_10K_R{1,2}_val_{1,2}.fq.gz`), SE, hardtrim5 and clock outputs **before** touching code; re-compare after | Identical md5s. The guard only ever rejects, so this is belt-and-braces — but V11 is on the do-not-skip list, so make it a measurement, not an inspection |
| **V12** | **Plain+gz FASTQ pair still accepted** — the §3.3 regression a naive `formats[0] != formats[1]` would cause | `--paired plainR1.fastq R2.fastq.gz` | Exit 0; `plainR1_val_1.fq` + `BS-seq_10K_R2_val_2.fq` written (verified working on `dev`) |
| **V13** | **`--clump_only` FASTQ-output paired keeps today's diagnosis** — closes §2.5 gap 2 and is the test that would have failed under v1 | `--clump_only --paired <fq> <bam>` and `--clump_only --paired <bam1> <bam2>`, both **without** `--output-format ubam` | Non-zero; stderr contains `--output-format ubam` and `drop aux tags`; must **not** suggest `--paired interleaved.bam` |
| **V14** | **The absent substring is pinned** — #363's substance | Mixed-pair test asserts `!stderr.contains("uBAM paired mode expects")` | Passes. Without this, the two messages can silently re-converge in a later refactor |
| **V15** | **No partial multi-pair output** | Run V8's command with `-o <fresh>/nested`; assert no `_val_*` and no `*_trimming_report.txt` for pair 1 | None present (they *are* present on `dev` — verified) |
| **V16** | Format labels are human-readable | Any mixed-pair stderr | Contains `uBAM`; does **not** contain `UnalignedBam` or `FastqGz` |
| **V17** | Site `:1268` covered — closes §2.5 gap 1 | `--paired <bam1> <bam2>` **without** `--output-format ubam`, as an integration test | Non-zero with the two-BAM message |

V3, V9, V11, V13, V15 are the ones that would quietly stop holding. Do not mark any of them "doesn't apply" — if a check appears not to apply, run it down. That discipline is what surfaced FastQC-Rust#6 and, this round, the silently-destructive `samtools` hint.

---

## 10. Questions or ambiguities

**Resolved since v1:**

- **R1 — fix shape for the Critical.** Reviewer A's mode enum adopted over Reviewer B's `fmt_flag` derivation (user-confirmed). B's version produces a *valid* command but still leads with "your pair is heterogeneous" when the user's actual problem is BAM-on-a-FASTQ-path — a point B's own report concedes while recommending otherwise. The enum preserves the aux-tag *reason*, which is the actionable part.
- **R2 — helper placement (v1 Q2).** `src/format.rs`. `main.rs` has zero `#[cfg(test)]` modules (verified), so a first bin-target unit test would be off-convention; `format.rs` owns `InputFormat` and has an existing test module.
- **R3 — `{fmt}` rendering.** Promote `input_format_label` (`clump_only.rs:685`) to `format.rs`. `{:?}` would leak enum variant names.
- **R4 — the `samtools` hint.** Replaced with `samtools merge -n`, verified end-to-end on samtools 1.21. v1's `collate` form silently drops half the library (§3.4).
- **R5 — dead code at `main.rs:498-512`.** Retired in this pass (user-approved).

**Open (assumption taken, no blocker):**

1. **Whether the mixed message distinguishes plain from gzipped FASTQ.** Taken: label them distinctly (`FASTQ (plain)` / `FASTQ (gzip)`) since `input_format_label` already does, but attribute the rejection to BAM-vs-FASTQ in the sentence so a plain+gz pair can never read as an error. V12 + V16 guard both halves.
2. **Pair index on single-pair runs.** Taken: omit when `n == 1`. House style at `main.rs:772-780` does print it, so either choice is defensible; "Pair 1 of 1" is noise on the common case.
3. **`Cli::validate_formats(&self, &[InputFormat])` as the eventual home** for all five format-dependent guards now inline in `main.rs` (`:211`, `:233`, `:251`, `:257`, and this one). `main.rs` is ~2 400 lines and validation logic is what a maintainer looks for in `cli.rs`. Out of scope here; recorded as the direction of travel.

**Critical:** none outstanding.

---

## 11. Self-Review

**Efficiency.** No hot path touched. Strictly less work on rejected runs, less I/O on all runs, and a correctness gain from single-source format classification (§6).

**Logic.** Traced against every paired entry point: trim FASTQ (`:1268`), trim uBAM output (`:1750`), clump-only uBAM (`:536-558`), clump-only FASTQ (`:426` N=1 + `clump_only.rs:399-410` N≥2), and the four specialty modes (exempt). Non-applicable paths checked rather than assumed: `--demux` is rejected with `--paired` at CLI level; `--passthrough` is rejected against any BAM at `:251-256`, i.e. earlier; `--retain_unpaired` has no format interaction. `run_paired` and `run_ubam_output` each have exactly one caller.

**Edge cases.** N=1 interleaved uBAM (exempt, A3); N=1 non-BAM (rejected at `:257`); odd N (rejected at `cli.rs:506`); R1==R2 (rejected at `cli.rs:513` — which is why `rejects_two_bam_paired` copies its fixture, per its own comment); **plain+gz pair (V12 — the regression a naive predicate would cause)**; multi-pair offence in a later pair (V8, V15); mixed in both argument orders (V1); `bgzip`-framed FASTQ (A1); empty `cli.input` (unreachable — `:176` indexes `[0]` first); empty BAM as R1 (keeps the empty-BAM error, correct precedence, §2.3).

**Integration.** Five user-visible changes enumerated in §7, up from three in v1; the multi-pair partial-output change is the largest and was missing.

**Remaining risks.**

- *Low:* wording drift breaks an existing test. Mitigated by V9, §3.5's line-break rule, and the fix-the-wording-not-the-test rule.
- *Low:* a future addition to §3.1's exemption list reopens the silent-wrong-output path at the two R1-header leaf functions. Mitigated by Step 5's belt-and-braces checks — which is why they are not optional.
- *Very low:* an unenumerated paired entry point still carries a stale guard. Mitigated by V5–V7 and by `grep -n 'two BAM files' src/` returning exactly one hit afterwards (Step 4's internal-invariant wording keeps that grep meaningful).

---

## 12. Revision history

### v2 — 2026-07-26, after dual independent plan review

**One Critical, from both reviewers independently, against v1's own claim.** v1 §11 recorded having traced the clump-only FASTQ-output path as *"unaffected — confirmed live, case D"*. Only the **N=1** shape is unaffected. For N≥2 the hoisted guard preempts `clump_only.rs:403-409` and replaces a correct, actionable message with a remediation `main.rs:426-432` rejects — the exact failure mode #363 is about, recreated in another mode. v1's live check exercised N=1 and generalised. Fixed by the `PairedShape` enum (§3.2, §5); pinned by V13, which had zero coverage before.

**Corrections to v1's own claims:**

| v1 said | Actually |
|---|---|
| clump-only FASTQ path "unaffected" | Only N=1; N≥2 was preempted (the Critical) |
| `samtools collate -O r1.bam r2.bam` interleaves two BAMs | Exits 0 and emits **half** the records — `collate`'s 2nd positional is the temp prefix. Replaced with verified `samtools merge -n` |
| §7 listed three behaviour changes | Five. The missing one — multi-pair runs no longer write partial output — is the largest |
| "Nothing has been created or read when the guard runs" | `sanity_check_any` and `detect_input_format` both read first; only "created" is load-bearing |
| Helper named for format equality | Predicate is BAM *count*; `InputFormat` has 3 variants and plain+gz pairs are legal (V12) |
| V8's `<fq> <fq> <fq> <bam>` tests pair indexing | Repeated path trips the collision pre-flight; post-fix it would flip to the guard and pass while testing the wrong thing |
| A6: "a mixed pair could be processed" | Only on `--cores > 1`; the default serial path is `FastqReader`-only |
| A4: specialty modes are harmless | True for hardtrim; `--clock`/`--implicon` do pair and are merely out of scope |
| Two stale comments to update | Four (`main.rs:1309-1310` and `:1993-1994` were missed) |
| Q2 (helper placement) open | Resolved: `format.rs`; `main.rs` has no test module |
| Line refs `:543`/`:552`, "~:300", `clump_only.rs:940-944` | `:542`/`:551`, `:321`, `clump_only.rs:942-945` |

**Added:** `PairedShape` enum; promoted `input_format_label`; enforced internal-invariant checks at three callees instead of comments; §2.4 pre-flight precedence; §3.5 Rust line-continuation rule; A7; V12–V17; retirement of dead code at `:498-512`.

### v1 — 2026-07-25

Initial plan. Correct on the core diagnosis (two defective sites, not the three the issue names), the guard position, and `Cli::validate()`'s preconditions — all independently confirmed by both reviewers.

---

## 13. Implementation notes

**Implemented 2026-07-26** on branch `fix/mixed-format-pair-message` off `dev` @ `a4ffd47`. All 12 outline steps done except Step 12 (follow-up issue), which is outward-facing and awaiting the user's go-ahead.

### What was built

| File | Change |
|---|---|
| `src/format.rs` | `PairedShape` enum + `reject_bam_format_mismatch_in_pair` + promoted `pub fn input_format_label`; 9 unit tests |
| `src/main.rs` | Guard wired in after the N=1 check; two defective sites replaced (one by an internal-invariant backstop, one deleted); dead Shape-B check retired; belt-and-braces check at `run_ubam_output_paired_two_files`; 3 stale comments rewritten |
| `src/clump_only.rs` | Private `input_format_label` removed in favour of the promoted one; belt-and-braces check + comment rewrite on Shape A |
| `tests/integration_paired_format_guard.rs` | **New**, 11 tests |
| `CHANGELOG.md` | Appended to the existing `#### Fixes` block |

Net: 441 tests, up from 421.

### Deviations from the plan

1. **Signature — `&[PathBuf]`, not `&[Path]`.** §5 specified `inputs: &[std::path::Path]`, which cannot compile: `Path` is unsized, so `[Path]` is not a valid slice element type. Used `&[std::path::PathBuf]`, which is what `cli.input` already is. Mechanical correction, no design impact.

2. **Pair prefix reads `"Got "` on a single pair, not `""`.** §3.2 said to omit the index when there is one pair. Implemented literally first, and the smoke test showed the result reads as a fragment: *"…`interleaved.bam`. two BAM files: A and B."* The clause needs a subject, so a single pair now yields *"Got two BAM files: A and B"* / *"Got mixed: …"* and multi-pair yields *"Pair 2 of 3 is …"*. Caught by reading the output rather than by a test — worth noting, since no assertion would have failed.

3. **Both defective sites did not get identical treatment.** §4 Step 4/5 anticipated this, but to be explicit: `run_paired`'s site was *replaced* with an enforced internal-invariant `bail!` (the sequential path below it hard-codes `FastqReader::open`, and `--cores 1` is the default), while `run_ubam_output`'s was *deleted* outright — its downstream leaf already receives an equivalent check.

4. **`grep -n 'two BAM files' src/` returns 4 lines, not 1.** §11's completion check expected one hit. All four are inside `format.rs`: one comment, two in the single message literal, one test assertion. The intent — no stale copy of the message outside the shared helper — holds; the literal count was mis-stated in the plan.

### Validation results

| Check | Result |
|---|---|
| V1 mixed message, both orders | Pass — names both inputs with `FASTQ (plain)` / `uBAM` labels |
| V2 two-BAM message + hint | Pass — `samtools merge -n`, both files named |
| V3 no side effects | Pass — output dir absent, no adapter-detection or banner lines. Automated |
| V4 uBAM-output path | Pass |
| V5 clump-only uBAM messages | Pass — remediation echoes `--clump_only --paired --output-format ubam` |
| V6 legal invocations | Pass — PE pair, single interleaved uBAM, multi-pair N=4 (4 outputs) |
| V7 specialty modes | Pass — hardtrim writes output and exits 0; `--clock` still reports its own read-count error |
| V8 multi-pair index | Pass — *"Pair 2 of 2 is mixed"*; three distinct input paths yielding four distinct outputs (see §13 deviation 5) |
| V9 existing tests unedited | Pass — 3/3, and `git diff --stat tests/` is empty |
| V10 suite + fmt + clippy | Pass — 441 tests, `fmt --check` clean, `clippy -D warnings` clean |
| V11 validation matrix | Pass — all 6 baseline outputs (PE ×2, SE, hardtrim5, clock ×2) byte-identical to pre-change md5s |
| V12 plain+gz accepted | Pass — unit + integration |
| V13 clump-only FASTQ diagnosis | Pass — all three shapes keep the aux-tag message; no `interleaved.bam` suggestion |
| V14 absent substring pinned | Pass — `!contains("uBAM paired mode expects")` on the mixed message |
| V15 no partial multi-pair output | Pass — output dir absent after an N=4 rejection |
| V16 human-readable labels | Pass — `uBAM` present, `UnalignedBam`/`FastqGz` absent |
| V17 site `:1268` covered | Pass — `two_bam_pair_keeps_interleave_remediation` |

### Iteration log

**#1 — `samtools` hint verified before shipping (pre-implementation).** The plan flagged v1's `samtools collate -O r1.bam r2.bam` as wrong on Reviewer A's reading that `collate` cannot take two inputs. Tested on samtools 1.21 against two 10 000-record uBAMs: it exits **0** and emits **10 000** records, silently consuming `r2.bam` as the temp prefix — worse than the error A predicted. `samtools merge -n -o interleaved.bam r1.bam r2.bam` yields all 20 000 records, mate-adjacent (flags 69/133), and TrimGalore accepts the result. An earlier attempt at this check was itself invalid: both fixtures were built with `samtools import -1`, so both claimed READ1 and TrimGalore rejected the merge for non-adjacent mates. Rebuilt with `-1`/`-2` before concluding anything.

**#2 — single-pair message read as a fragment.** See deviation 2. Changed the empty prefix to `"Got "`; re-verified all seven shapes by hand.

**#3 — `cargo fmt` reflowed three assertions in `format.rs`.** Cosmetic; clippy was already clean. Re-ran the full suite after formatting — 441 still passing.

### Follow-ups

- **Not done: Step 12**, the follow-up issue for `--hardtrim5/3` silently accepting mixed pairs and `--clock` reporting a misleading read-count error (§2.6). Awaiting go-ahead — filing an issue is outward-facing.
- **Not done:** the correcting comment on #363 itself, noting that `main.rs:542` was the reference implementation rather than a third defect (§2.2). Same reason.
- Nothing committed or pushed.

---

## 14. Verification round — dual code review + coverage audit

Ran 2026-07-26 after implementation: two independent code reviewers plus a plan-manager coverage audit, all in fresh contexts, each told to treat §13 as a claim to verify.

- `CODE_review_reviewer-A.md`, `CODE_review_reviewer-B.md`, `COVERAGE.md`.
- **Coverage verdict: INCOMPLETE — 1 item**, that item being Step 12 (the user-gated follow-up issue). 45 items audited: 42 DONE, 0 PARTIAL, 1 MISSING, 2 DEVIATED (both already documented). Every §9 check was independently re-run and matched §13 in every row.
- **V11 was upgraded from deductive to empirical.** Reviewer A built `dev` @ `a4ffd47` in a throwaway worktree with a separate `--target-dir` and md5-compared 15 outputs across PE, SE, `--hardtrim5`, `--clock`, `--demux` and `--clump_only`: identical. That closes the gap the coverage audit had declared (it could not reach the pre-change baselines).

### Defects found in the implementation, and fixed

| # | Found by | Defect |
|---|---|---|
| 1 | **Both** | **The CHANGELOG asserted a bug that never shipped.** The entry claimed the old two-BAM message "suggested `samtools collate -O r1.bam r2.bam`". It never did — `git show a4ffd47:src/main.rs` shows all three shipped two-BAM messages ending at *"Got two BAM files; one of them is {}"* with **no** `samtools` command. The destructive two-input `collate` form existed only in this plan's v1 draft, corrected before implementation. Shipping the paragraph would have told users TrimGalore once printed a silently-lossy command, and could have cast doubt on the *correct* one-input `collate` hint in `bam.rs:57`. Rewritten to describe what actually changed. Reviewer A also caught a second false claim in the same paragraph — that `--clock`/`--implicon` "accept mixed pairs today", which they do not |
| 2 | **Both** | The internal-invariant message named `run_paired_end`; the enclosing function is `run_paired` (`trimmer::run_paired_end` is a different function named in `CLAUDE.md`). §4 Step 4 specified the wrong name verbatim, so the plan carried it in. Corrected in both plan and code |
| 3 | B | **A fifth stale comment**, missed by Step 6's list of four: `run_ubam_output_paired_two_files` still carried *"Two-file paired-BAM rejection lives up-front in `run_ubam_output`"* — pointing at the loop this change deletes, three lines above the new backstop comment that says the opposite. Deleted |
| 4 | B | `multi_pair_names_the_offending_pair_and_writes_nothing` asserted `contains("Pair 2 of 2")`, which the per-pair progress banner `=== Pair 2 of 2 ===` also satisfies — so a regression where the guard stopped firing and the run *proceeded* would still have matched. Tightened to `"Pair 2 of 2 is mixed"`, which only the guard can produce |

### Recommendations adopted

| # | Found by | Change |
|---|---|---|
| 5 | **Both** | `debug_assert_eq!` → `anyhow::ensure!` + `chunks_exact(2)`. Release builds set no `debug-assertions`, so the assert compiled out of every shipped binary, and Reviewer A identified the worse of the two failure modes: a short `formats` makes `zip` truncate, so later pairs go unexamined and the function returns `Ok(())` **without having checked** — a validation guard that silently passes, precisely the failure class it exists to prevent. Now also checks evenness. Deviates from §5/A7, which specified `debug_assert` |
| 6 | **Both** | Both leaf backstops widened from *mismatch* to *any BAM*. The code they replaced rejected any BAM; narrowing to a mismatch check half-met the stated intent, since a two-BAM pair reaching either leaf would silently discard R2's `@HD`/`@PG` chain and tag dictionary — the same "header from R1 only" hazard the backstops exist for. `n_bam == 0` is the invariant main()'s guard actually establishes for those shapes |
| 7 | A | The negative `samtools` assertion pinned the full `collate -O <fixture> <fixture>` string — a spelling that never existed and that a fixture rename would have turned into a tautology. Now asserts the property: `!contains("samtools collate")` |
| 8 | B | The CHANGELOG's three consequences were stated as applying to all four shapes. Only `--output_dir` does. Adapter detection ran before the old rejection **only** on the ordinary trim path, and the no-partial-output change applies to trim and `--clump_only` FASTQ but not the two uBAM-out paths, which already validated every pair up-front. Scoped per bullet. **Same generalisation error §2.1 records as one of v1's, recurring in a different file** |
| 9 | A | `anyhow::bail!` → `bail!` in `clump_only.rs`, which imports `bail` and uses the bare form in the adjacent Shape B backstop |

### Plan-text corrections from this round

`{fmt_hint}` removed from §3.4 (a placeholder the plan never defined; the implementation renders it empty — flagged by both reviewers *and* the coverage audit); §2.1/§4/§11/§13 `run_paired_end` → `run_paired`; §9 V8's "four distinct paths" restated as the property that actually matters, distinct **output** paths (the test passes three distinct inputs, and Reviewer A independently confirmed with four genuinely distinct paths that the index is still correct, so no test change was needed); §2.3/§6's "three deleted re-detections" corrected to two of three, since `:1269`'s was deliberately retained as the backstop.

### Deviations added to the §13 list

5. **`{fmt_hint}` dropped** from the mixed message (§3.4). Harmless — both legal alternatives are already named for every shape that emits the message — but it was an undocumented divergence from the plan text.
6. **`anyhow::ensure!` instead of `debug_assert`** (§5, A7). See item 5 above; the plan's choice would have compiled out of release builds.
7. **Leaf backstops enforce `n_bam == 0`, not mismatch** (§4 Step 5). See item 6 above.

### Not adopted, surfaced to the maintainer instead

- **`docs/src/content/docs/quickstart.md:57` states *"Paired reads may come as two BAM files"*, which is false** — the binary rejects it, and this change makes the rejection more emphatic. Pre-existing on `dev` and outside this diff, so Reviewer A deliberately did not touch it. #363 is arguably the right issue to close it under.
- **Multi-pair `ClumpOnlyFastqOut` omits the pair index** (B Low-3): `pair_prefix` is computed but that arm bails before using it. A direct consequence of §3.4's decision to preserve the message verbatim, which is what makes V13 a pure regression test. Left as-is deliberately; changing it would weaken that property.
- Editorial nits: `Got mixed:` → `Got a mixed pair:` (A), `Pair 2 of 2 is two BAM files` → `contains two BAM files` (B). Both grammatical as they stand; the shared prefix makes the second awkward to apply.
- **`Cli::validate_formats(&self, &[InputFormat])`** as the eventual home for all five format-dependent guards now inline in `main.rs` (§10 open item 3). Still out of scope.

### Process note

Both reviewers edited the same working tree concurrently; two of B's `Edit` calls collided with A's, and one fix was applied by both. Reviewer B correctly warned that neither report is a complete account of the tree. Every fix above was therefore re-verified directly against the working tree rather than taken from either report, and the full suite, `fmt`, `clippy` and the V11 md5 comparison were re-run after the final edit: **441 tests passing, fmt and clippy clean, all 6 baseline outputs still byte-identical.** For future rounds, dual reviewers should either review read-only or work on separate worktrees.
