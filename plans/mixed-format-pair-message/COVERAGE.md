# Plan Coverage Report

**Mode:** B (code vs. plan §4/§3/§9)
**Plan(s):** plans/mixed-format-pair-message/PLAN.md (v2)
**Date:** 2026-07-26
**Verdict:** INCOMPLETE — 1 item unresolved (§4 Step 12, deliberately user-gated)

## Summary

- Total items: 45
- DONE: 42
- PARTIAL: 0
- MISSING: 1
- DEVIATED: 2 (both documented in §13)

Audit basis: branch `fix/mixed-format-pair-message`, uncommitted working tree, based on `dev` @ `a4ffd47`. Every §9 check was re-run independently rather than read from §13. All 441 tests pass; `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings` clean **from a cold scratch target dir** (not a cached result — exit 0, zero warning/error lines).

## Coverage ledger

### §4 Implementation outline (12 steps)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | Shape enum + pure decision helper added to `src/format.rs` (not `main.rs`), `anyhow::Result` signature | §4 Step 1 | DEVIATED | Both exist: `PairedShape` at `src/format.rs:54`, `reject_bam_format_mismatch_in_pair` at `:108`. Signature takes `&[std::path::PathBuf]`, not §5's `&[std::path::Path]`. Documented as deviation 1 in §13, and mechanically forced — `Path` is unsized, so `[Path]` is not a valid slice element type. No design impact; `cli.input` is already `Vec<PathBuf>`. `debug_assert_eq!` for A7 present at `:113` |
| 2 | Promote `input_format_label` to `format.rs`; update `clump_only.rs` use site | §4 Step 2 | DONE | `pub fn input_format_label` at `src/format.rs:41`, body byte-identical to the removed private version (same three arms, same strings). Private copy deleted from `clump_only.rs`; import updated at `clump_only.rs:51`; both call sites (`:783`, `:1016`) now resolve to the promoted fn. Verified live: a `--clump_only` run still renders the label in `*_clumping_report.txt` |
| 3 | Call the helper in `main()` after the N=1 check, before the `gzip` binding; gate per §3.1; resolve shape per §3.2; comment in the `--phred64` register | §4 Step 3 | DONE | Guard at `src/main.rs:286-301`, comment block at `:267-285`. Ordering confirmed by line number: `sanity_check_any` `:178` → `input_formats` `:184` → `--phred64` → `--preserve-tags` `:235` → `--passthrough` `:253` → N=1 `:259` → **guard `:286`** → `gzip` `:311` → `ensure_output_dir` `:315`. Comment covers both "why here" and "why specialty modes are exempt" |
| 4 | Replace defective site `:1268-1278` with an *enforced* internal-invariant check using internal-invariant wording | §4 Step 4 | DONE | `src/main.rs:1278-1296`. Still a `bail!`, not a comment. Wording: `"internal error: uBAM input reached run_paired (…); the paired format guard in main() should have rejected it. Please report this at …"` — names the real function (`run_paired`, not §4's illustrative `run_paired_end`) and cannot be mistaken for user-facing text |
| 5 | Delete site `:1750-1762` and the redundant clump-only branch `:536-558`, each replaced by a pointer; keep belt-and-braces at the two R1-header leaf functions | §4 Step 5 | DONE | Both deletions done, both replaced by pointer comments (`main.rs:1767-1771`, `main.rs:564-568`). Belt-and-braces present at **both** required leaves: `run_ubam_output_paired_two_files` (`main.rs:1999-2016`) and `clump_only_paired_to_bam_one_pair` Shape A (`clump_only.rs:937-955`). Each `bail!`s on `is_bam(r1) != is_bam(r2)` |
| 6 | Update **all four** stale comments | §4 Step 6 | DONE | All four rewritten: `main.rs:1278-1287` (was "handled in main() before this fn is called"), `main.rs:1327-1333` (was "rejected above"), `main.rs:1999-2010` (was "At this point both inputs are FASTQ"), `clump_only.rs:933-941` (was "Format-guards in main.rs::dispatch"). §13 says "3 stale comments rewritten", which counts only the `main.rs` ones; the fourth is in the `clump_only.rs` row of the same table. All four are done |
| 7 | Retire the dead check at `main.rs:498-512` | §4 Step 7 | DONE | Replaced by an explanatory note at `main.rs:536-542`. Unreachability re-verified live, not assumed: `--clump_only --paired --output-format ubam <fastq>` now reports the N=1 guard's message (`"--paired with a single input file is only legal if that file is a uBAM"`), confirming the retired branch's condition was a strict subset |
| 8 | Unit tests in `format.rs`'s test module — 8 enumerated behaviours | §4 Step 8 | DONE | 9 tests, all 8 behaviours covered. See test table. All-FASTQ pass (across all 4 shapes), plain+gz pass, mixed both orders, two-BAM, `ClumpOnlyFastqOut` on `n_bam == 1` *and* `== 2`, multi-pair index, single-pair index omitted. Plus a 9th (`pair_guard_remediation_echoes_the_users_mode`) not required by the plan |
| 9 | Integration tests per §9, covering both §2.5 gaps | §4 Step 9 | DONE | `tests/integration_paired_format_guard.rs` (new, untracked), 11 tests. Gap 1 (plain `--paired <bam> <bam>`, no `--output-format ubam`) → `two_bam_pair_keeps_interleave_remediation`. Gap 2 (`--clump_only --paired` FASTQ output, N≥2) → `clump_only_fastq_output_keeps_aux_tag_diagnosis`, all three shapes |
| 10 | Verify the three existing tests pass **unedited** | §4 Step 10 | DONE | Verified two ways, not one: `git diff tests/` is **empty** (no tracked test file touched — the only entry under `tests/` in `git status` is the new untracked file), and all three pass when run individually. See V9 |
| 11 | CHANGELOG appended to the **existing** `#### Fixes` block under `### Unreleased`; cite #363; cover 4 points | §4 Step 11 | DONE | Entry at `CHANGELOG.md:130`, inside the single `#### Fixes` block at `:84`. Confirmed there is exactly one `####` Fixes heading between `### Unreleased` (`:4`) and `### Version 2.3.0` (`:165`) — no second block added. All four required points present: wording change, earlier firing, no output dir created, and multi-pair no-partial-output (called out in bold). #363 linked |
| 12 | File the follow-up issue for §2.6 (specialty modes on mixed pairs), cross-referencing #363 | §4 Step 12 | **MISSING** | No such issue exists. `gh issue list --repo FelixKrueger/TrimGalore --state open` returns exactly one open issue: #363 itself. Self-reported as not done in §13 "Follow-ups", awaiting the user's go-ahead because filing is outward-facing |

### §3 Behavior

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 13 | Guard activation: `cli.paired` ∧ `len > 1` ∧ no specialty mode | §3.1 | DONE | `main.rs:286-292`, exactly the five conjuncts specified (`hardtrim5.is_none()`, `hardtrim3.is_none()`, `!clock`, `implicon.is_none()`). `clump_only` correctly **not** in the exemption list |
| 14 | Shape determination (4 shapes) + per-pair decision; `ClumpOnlyFastqOut` predicate is `n_bam >= 1`; others 0/1/2; first offending pair wins; index only when >1 pair | §3.2 | DEVIATED | Shape resolution at `main.rs:294-299` matches the §3.2 table exactly. `ClumpOnlyFastqOut` is checked at `format.rs:142` **before** the `n_bam == 2` branch, so it fires on 1 *or* 2 BAMs — verified live for all three shapes. Loop `bail!`s on the first pair with `n_bam > 0`. **Deviation:** the single-pair prefix is `"Got "`, not `""`. Documented as deviation 2 in §13, with the reason (the clause needs a subject; the literal reading produced *"…`interleaved.bam`. two BAM files: A and B."*). §3.2's actual requirement — omit "Pair 1 of 1" — is met, and pinned by `pair_guard_omits_index_for_a_single_pair` |
| 15 | Predicate is BAM-count, not format equality | §3.3 | DONE | `n_bam = fmts.iter().filter(is_bam).count()` at `format.rs:123`. No `formats[0] != formats[1]` anywhere. Confirmed live: plain+gz pair exits 0 (V12), bgzip-framed FASTQ + gz pair exits 0 (A1) |
| 16 | Four message variants: two-BAM, mixed, `ClumpOnlyFastqOut` verbatim, all ✚ substrings retained | §3.4 | DONE | All four rendered and inspected live (V1/V2/V4/V5/V13). Two-BAM retains `two BAM files is not supported`, `uBAM paired mode expects`, `single interleaved`, names **both** files, and carries the verified `samtools merge -n -o interleaved.bam <r1> <r2>`. Mixed retains `same format` and `mixed`, names both files with labels, and has **no** bare interleaved-file imperative. `ClumpOnlyFastqOut` is verbatim — checked against `git show HEAD:src/clump_only.rs`, which bailed on the first BAM of `[r1, r2]` with the same three clauses; the new `first_bam` selection reproduces that choice. *One discrepancy in the plan, not the code:* §3.4's mixed template ends `"…a single interleaved uBAM{fmt_hint}."` but `{fmt_hint}` is never defined anywhere in the plan; the implementation renders it as empty. Nothing to implement, and nothing testable was dropped |
| 17 | Rust line-continuation rule: no pinned substring straddles a source line break | §3.5 | DONE | Source-verified at `format.rs:159-180`: each of the four phrases sits wholly on one line, and every break carries its space *before* the `\`. Also proven at the stronger level — the rendered stderr reads `"…is not supported. uBAM paired mode expects a single interleaved file: …"` with correct spacing, so no word-joining occurred. The rule is additionally recorded as a maintenance comment at `format.rs:153-158` |
| 18 | No side effects on rejection; non-zero exit | §3.6 | DONE | No reader opened, no adapter detection, no poly-G, no banner, no output dir, no pair processed — all confirmed live and by an automated test whose absence-assertion is non-vacuous (see V3). Non-zero exit on every rejection variant |

### §9 Validation (V1–V17) — each re-run independently

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 19 | V1 — mixed pair, both orders, gets the mixed message naming both files with formats | §9 V1 | DONE | Ran both orders. `--paired requires both inputs of a pair to be the same format. Got mixed: test_files/phred64_test.fastq is FASTQ (plain) and test_files/ubam_test.bam is uBAM. …` and the reverse. Non-zero both times |
| 20 | V2 — two distinct BAMs keep the two-BAM message + interleaved remediation + `samtools merge -n` hint | §9 V2 | DONE | Full message reproduced, including `samtools merge -n -o interleaved.bam test_files/ubam_test.bam test_files/ubam_paired_test.bam` |
| 21 | V3 — **no side effects**, checked non-vacuously | §9 V3 | DONE | **Helper checked, per instruction.** `tempdir()` (`integration_paired_format_guard.rs:37`) does `create_dir_all`, but the `-o` argument comes from `nonexistent_out()` (`:46`), which returns `dir.join("nested_out")` and **asserts the precondition that it does not exist** before use. The `!out.exists()` assertion is therefore real. Independently reproduced on the binary: `-o <tmp>/nested` left the path absent, stderr free of adapter-detection/poly-G/banner lines |
| 22 | V4 — uBAM-output path fixed too, mode-appropriate remediation | §9 V4 | DONE | Mixed + `--output-format ubam` → mixed message. Two-BAM + `--output-format ubam` → remediation reads `trim_galore --paired --output-format ubam interleaved.bam` |
| 23 | V5 — clump-only uBAM-output messages did not regress after deleting `:536-558` | §9 V5 | DONE | Mixed → mixed message; two-BAM → two-BAM message with remediation `trim_galore --clump_only --paired --output-format ubam interleaved.bam` |
| 24 | V6 — legal invocations unaffected | §9 V6 | DONE | Paired FASTQ pair → exit 0, both `_val_{1,2}.fq.gz` + both reports. `--paired <interleaved.bam>` → exit 0, `_val_1.fq`/`_val_2.fq`. Multi-pair N=4 with four distinct paths → exit 0, four `_val_*` files |
| 25 | V7 — specialty modes untouched | §9 V7 | DONE | Tested **all four**, not just the two in the plan's command column. `--hardtrim5 10` → exit 0, both outputs written (including from the BAM input). `--hardtrim3 5` → exit 0, both outputs. `--clock` → still `Paired-end files have different numbers of reads!`. `--implicon` → same pre-existing error. Guard does not fire for any of them |
| 26 | V8 — multi-pair offence names the right pair, four distinct paths | §9 V8 | DONE | `Pair 2 of 2 is mixed: …p2_R1.fastq.gz is FASTQ (gzip) and test_files/ubam_test.bam is uBAM.` Four distinct paths used, so this tests indexing rather than collision precedence |
| 27 | V9 — three existing tests pass **unedited** | §9 V9 | DONE | **Verified by diff, not only by passing.** `git diff tests/` → empty. `rejects_two_bam_paired` ok, `rejects_mixed_format_paired` ok, `ubam_out_two_bam_pair_rejected` ok (1 passed each, run individually) |
| 28 | V10 — full suite + fmt + clippy | §9 V10 | DONE | `cargo test --release` → 441 passed, 0 failed across 8 test binaries (369 lib + 12 + 15 + 1 + 11 + 2 + 8 + 23). `cargo fmt --all -- --check` clean. `cargo clippy --all-targets --release -- -D warnings` clean — re-run **cold in a scratch `--target-dir`** so the result cannot be a stale cache: exit 0, zero warning/error lines |
| 29 | V11 — Perl byte-identity paths unmoved | §9 V11 | DONE | See "Verification limitations" — confirmed by exhaustive diff classification plus a clean run of all five CI matrix paths; the plan's claimed pre-change md5 baselines were not accessible to me |
| 30 | V12 — plain+gz FASTQ pair still accepted | §9 V12 | DONE | Exit 0; `plain_R1_val_1.fq` + `BS-seq_10K_R2_val_2.fq` written. Also covered by a unit test and an integration test |
| 31 | V13 — `--clump_only` FASTQ-output paired keeps today's diagnosis | §9 V13 | DONE | All three shapes (fq+bam, bam+fq, bam+bam) produce `uBAM input under --clump_only requires --output-format ubam (using the FASTQ output path with uBAM input would drop aux tags). Input: test_files/ubam_test.bam`. Contains `--output-format ubam` and `drop aux tags`; contains **no** `interleaved.bam` |
| 32 | V14 — the absent substring is pinned | §9 V14 | DONE | Pinned in two places: `format.rs` unit test (`!msg.contains("uBAM paired mode expects")`) and `integration_paired_format_guard.rs:93` (plus a second assertion that the mixed message contains no `two BAM files` at all). Confirmed live |
| 33 | V15 — no partial multi-pair output | §9 V15 | DONE | Ran V8's command with `-o <fresh>/nested`: directory absent, zero files. Non-vacuous for the same reason as V3 |
| 34 | V16 — format labels are human-readable | §9 V16 | DONE | Live stderr contains `uBAM`, `FASTQ (plain)`, `FASTQ (gzip)`; contains no `UnalignedBam`, `FastqGz`, or `FastqPlain`. Asserted in both the unit and integration tests |
| 35 | V17 — site `:1268` covered by an integration test | §9 V17 | DONE | `two_bam_pair_keeps_interleave_remediation` runs `--paired <bam1> <bam2>` with no `--output-format ubam` and asserts the two-BAM message. Also reproduced manually |

### §11 enumerated edge cases

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 36 | N=1 interleaved uBAM — exempt (A3) | §11 | DONE | `--paired test_files/ubam_paired_test.bam` → exit 0, `_val_{1,2}.fq`. Guard's `len > 1` conjunct keeps it out |
| 37 | N=1 non-BAM — rejected by the pre-existing guard | §11 | DONE | Reports `--paired with a single input file is only legal if that file is a uBAM.` (`main.rs:259`), not a guard message |
| 38 | Odd N — rejected in `Cli::validate()` | §11 | DONE | `Paired-end mode requires an even number of input files (R1/R2 pairs), got 3` — fires before the guard, so `chunks(2)` never sees a short final chunk |
| 39 | R1 == R2 — rejected in `Cli::validate()` | §11 | DONE | `Read 1 and Read 2 appear to be the same file: test_files/ubam_test.bam` — precedes the guard even when both are BAM |
| 40 | plain+gz pair | §11 | DONE | Same as V12 |
| 41 | Multi-pair offence in a later pair | §11 | DONE | Same as V8/V15 |
| 42 | Mixed in both argument orders | §11 | DONE | Same as V1 |
| 43 | bgzip-framed FASTQ (A1) | §11 | DONE | Tested end-to-end with real `bgzip`, both directions: bgzip-FASTQ + gz-FASTQ → **accepted**, exit 0 with `bg_R1_val_1.fq.gz`; bgzip-FASTQ + real uBAM → **mixed** message labelling the bgzip file `FASTQ (gzip)`. The `BAM\1`-in-payload discriminator holds; the guard adds no format logic of its own |
| 44 | Empty `cli.input` — unreachable | §11 | DONE | `sanity_check_any(&cli.input[0])` at `main.rs:178` indexes `[0]` before the guard, and clap requires ≥1 input. Guard is not reachable with an empty vector |
| 45 | Empty BAM as R1 in a mixed pair — empty-BAM error must win | §11 | DONE | `Error: Input file '…/empty.bam' is empty`, from `sanity_check_any` at `:178`. Correct precedence per §2.3 — the guard did not preempt it |

## Gaps (detail)

### Item 12 — §4 Step 12: file the follow-up issue for §2.6

**Expected:** A GitHub issue filed against `FelixKrueger/TrimGalore` covering the two out-of-scope specialty-mode defects — `--hardtrim5/3 --paired <fastq> <bam>` silently accepting the mixed pair and writing output, and `--clock`/`--implicon --paired <fastq> <bam>` failing with the misleading `Paired-end files have different numbers of reads!` (which would have *proceeded* had the record counts matched) — cross-referencing #363.

**Found:** Nothing. `gh issue list --repo FelixKrueger/TrimGalore --state open --limit 15` returns exactly one open issue, #363. No follow-up issue exists in any state that the listing surfaces.

**Gap:** File the issue. Both defects were re-confirmed live during this audit (item 25), so the content is ready: `--hardtrim5 10 --paired phred64_test.fastq ubam_test.bam` exits 0 and writes `phred64_test.10bp_5prime.fq` + `ubam_test.10bp_5prime.fq`; `--hardtrim3 5` behaves the same; `--clock` and `--implicon` both emit the read-count error. Note that §8 A4 records this issue as *"doing real work and not optional"* — the `--clock`/`--implicon` exemption is "out of scope", not "correct as-is".

This is the only unresolved plan item. §13 parks it explicitly as awaiting the user's go-ahead because filing is outward-facing, which is a reasonable hold — but the plan lists it, so it is MISSING until filed.

### Item 1 — §5 signature deviation (documented)

**Expected:** `inputs: &[std::path::Path]`.
**Found:** `inputs: &[std::path::PathBuf]` (`src/format.rs:109`).
**Gap:** None. The plan's signature cannot compile — `Path` is unsized, so `[Path]` is not a valid slice element type. Documented as deviation 1 in §13. `PathBuf` is what `cli.input` already holds, so no conversion is introduced at the call site.

### Item 14 — §3.2 single-pair prefix deviation (documented)

**Expected:** §3.2 — "The message names the 1-based pair index **when there is more than one pair**", read literally as an empty prefix on a single pair.
**Found:** `"Got "` on a single pair (`src/format.rs:134`); `"Pair {i} of {n} is "` on multi-pair (`:132`).
**Gap:** None. Documented as deviation 2 in §13, with the reason: the literal reading produced a sentence fragment (*"…`interleaved.bam`. two BAM files: A and B."*). The substantive requirement — no `"Pair 1 of 1"` noise — is met and is pinned by `pair_guard_omits_index_for_a_single_pair`. Worth noting that §13 is right that no assertion would have caught this; it was found by reading the output.

## Test verification

| Test name | File | Status |
|-----------|------|--------|
| `pair_guard_accepts_all_fastq` | `src/format.rs` | PASS |
| `pair_guard_accepts_plain_plus_gzip_fastq` | `src/format.rs` | PASS |
| `pair_guard_rejects_mixed_in_both_orders` | `src/format.rs` | PASS |
| `pair_guard_rejects_two_bam_with_interleave_remediation` | `src/format.rs` | PASS |
| `pair_guard_clump_only_fastq_out_keeps_aux_tag_diagnosis` | `src/format.rs` | PASS |
| `pair_guard_clump_only_fastq_out_names_the_bam_side` | `src/format.rs` | PASS |
| `pair_guard_reports_the_offending_pair_index` | `src/format.rs` | PASS |
| `pair_guard_omits_index_for_a_single_pair` | `src/format.rs` | PASS |
| `pair_guard_remediation_echoes_the_users_mode` | `src/format.rs` | PASS |
| `mixed_pair_reports_format_mismatch_not_two_bams` | `tests/integration_paired_format_guard.rs` | PASS |
| `mixed_pair_rejection_has_no_side_effects` | `tests/integration_paired_format_guard.rs` | PASS |
| `two_bam_pair_keeps_interleave_remediation` | `tests/integration_paired_format_guard.rs` | PASS |
| `two_bam_pair_rejected_with_ubam_output_and_echoes_mode` | `tests/integration_paired_format_guard.rs` | PASS |
| `mixed_pair_rejected_with_ubam_output` | `tests/integration_paired_format_guard.rs` | PASS |
| `clump_only_fastq_output_keeps_aux_tag_diagnosis` | `tests/integration_paired_format_guard.rs` | PASS |
| `clump_only_ubam_output_distinguishes_mixed_from_two_bam` | `tests/integration_paired_format_guard.rs` | PASS |
| `multi_pair_names_the_offending_pair_and_writes_nothing` | `tests/integration_paired_format_guard.rs` | PASS |
| `plain_plus_gzip_fastq_pair_is_still_accepted` | `tests/integration_paired_format_guard.rs` | PASS |
| `single_interleaved_ubam_still_accepted` | `tests/integration_paired_format_guard.rs` | PASS |
| `hardtrim_still_accepts_a_mixed_pair` | `tests/integration_paired_format_guard.rs` | PASS |
| `rejects_two_bam_paired` (pre-existing, unedited) | `tests/integration_clump_only_ubam.rs` | PASS |
| `rejects_mixed_format_paired` (pre-existing, unedited) | `tests/integration_clump_only_ubam.rs` | PASS |
| `ubam_out_two_bam_pair_rejected` (pre-existing, unedited) | `tests/integration_ubam_out.rs` | PASS |
| Full suite | 8 test binaries + doctests | **441 passed, 0 failed, 0 ignored** |

`git status` confirms `tests/integration_paired_format_guard.rs` is the only new file under `tests/` and that no tracked test file was modified (`git diff tests/` is empty).

## Verification limitations

Three things I could not confirm to the same standard as the rest, stated so they are not read as verified:

1. **V11 (Perl byte-identity) — confirmed deductively and by clean runs, not against the plan's baselines.** The pre-change md5s §13 reports were not available to me, and a permission prompt blocked building a pre-change binary out-of-tree. What I did instead:

   - **Classified every hunk of the diff for accepted-path effect.** `src/format.rs` is pure addition — no existing function altered. Every piece of new code in `main.rs`/`clump_only.rs` is `bail!`-or-fall-through: the hoisted guard, the `run_paired` internal-invariant check, and both belt-and-braces checks contribute nothing on an accepted run. Every deletion (`main.rs:498-512`, `:536-558`, `:1750-1762`) removed `bail!`-only code plus some redundant `detect_input_format` reads. The one moved function, `input_format_label`, has a byte-identical body. The one behavioural refinement — `run_ubam_output_paired_two_files` now reads `source_header` from a cached `fmt_r1` rather than re-calling `detect_input_format` — yields the same value. No accepted-path output byte can change.
   - **Checked the guard cannot newly reject anything.** Every paired N≥2 shape containing a BAM was *already* rejected before this change (trim `:1268`, uBAM-out `:1750`, clump-only uBAM-out `:536-558`, clump-only FASTQ-out `clump_only.rs:401`, confirmed against `git show HEAD:src/clump_only.rs`). The CI matrix paths themselves contain no BAM input, and SE/demux are not `--paired` while hardtrim5/clock are exempt.
   - **Ran all five CI matrix paths clean** and recorded the current md5s of the gzip-decompressed outputs so a future comparison has a written anchor: SE `illumina_10K_trimmed` `ffc7cc762bbca5317074e3c3d2c3dce0`; PE `BS-seq_10K_R1_val_1` `12fde5c81a6780ecfa57e54b7953894a`, `_R2_val_2` `6c3857913950a39a6e5ecef10a637129`; hardtrim5 `illumina_10K.30bp_5prime` `cf1221cacdc8c0e5e101e6fe5bb7d9e5`; clock `clock_10K_R1.clock_UMI.R1` `0d2730066e76adba9d08540a27538135`, `…R2` `92418fb4b14a0f39f759ffaccbd32dd3`; demux 8 files incl. `demux_test_trimmed` `50d1c4ffc8ca7b23b6e17ae055f356bd`.
   - **Noted supporting empirical coverage.** The committed pre-change reference fixtures (`ubam_out_se_REFERENCE.bam`, `ubam_out_pe_REFERENCE.bam`, `ubam_paired_test_val_{1,2}_REFERENCE.fq`) are compared byte-wise by tests that pass, and `ubam_out_pe_two_fastq_interleaved` exercises `run_ubam_output_paired_two_files` — the one function that gained a real code change on an accepted path — and passes.

   The authoritative check remains the CI `validation` job against Perl 0.6.11.

2. **§3.4's `{fmt_hint}`.** The plan's mixed-message template contains a `{fmt_hint}` placeholder that the plan never defines. The implementation renders it empty. Nothing verifiable was dropped, but the plan text and the code do not literally agree, and a future reader of §3.4 may look for it.

3. **§11's completion grep.** `grep -n 'two BAM files' src/` returns 4 lines, not the 1 §11 predicted — §13 deviation 4 reports this accurately. I confirmed all four are inside `src/format.rs` (one maintenance comment at `:154`, two inside the single message literal at `:161`/`:164`, one test assertion at `:431`) and that `grep -rln` matches `src/format.rs` and nothing else. The check's *intent* — no stale copy of the message outside the shared helper — holds.

Separately, the §13 follow-up noting that a correcting comment should be posted on **#363** (recording that `main.rs:542` was the reference implementation, not a third defect, per §2.2) has **not** been done: `gh issue view 363` shows `state=OPEN, comment_count=0`. This is not a §4 step, so it is not counted in the ledger, but it is outstanding. Per instruction I did not file or comment on anything.

## Verdict

**INCOMPLETE — 1 item unresolved.**

All code, test, and documentation work in the plan is complete and independently verified. Every one of §4's steps 1–11, all six §3 behaviour requirements, all seventeen §9 validation checks, and all ten §11 edge cases hold. The two deviations are both documented in §13, both immaterial, and one of them (the `&[Path]` signature) was forced by the language. §13's self-reported validation table matched what I measured in every row.

One item remains:

- **§4 Step 12 — file the follow-up GitHub issue for §2.6.** Content: `--hardtrim5/3 --paired <fastq> <bam>` silently accepts the mixed pair and writes output (re-confirmed live: exit 0, both outputs written, for both `--hardtrim5 10` and `--hardtrim3 5`); `--clock --paired <fastq> <bam>` and `--implicon --paired <fastq> <bam>` fail with `Paired-end files have different numbers of reads!`, misleading in the same way as #363 and would have proceeded had the record counts matched. Cross-reference #363. §8 A4 records this as required work rather than optional, because the `--clock`/`--implicon` exemption is "out of scope", not "correct as-is".

Also outstanding, though not a numbered plan step: the correcting comment on #363 itself (§2.2 — `main.rs:542` was the reference implementation, not a third cause site). #363 currently has zero comments.

Both are outward-facing actions that §13 deliberately parked pending the user's go-ahead. Nothing in the code needs to change.
