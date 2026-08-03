# PLAN review — Reviewer A

**Plan:** `plans/mixed-format-pair-message/PLAN.md` (issue #363)
**Repo state reviewed:** `dev` @ `a4ffd47`, release binary rebuilt from that tree
**Method:** every load-bearing `file:line` claim re-read in the source; nine invocations run against the freshly built binary; scratch output written to `$TMPDIR`, no repo file touched except this report.

Verdict: **the diagnosis and the guard placement are correct, and most of the plan's claims check out exactly.** One finding is Critical — the plan's own §11 "unaffected" claim about the clump-only FASTQ-output path is wrong, and the change as specified would hand those users a remediation command that is guaranteed to fail. That is the same defect class #363 is about, re-introduced in a path with zero test coverage.

---

## 0. What I verified as sound (stated briefly, not re-litigated)

| Plan claim | How verified | Result |
|---|---|---|
| Two defective sites at `main.rs:1268-1278` and `main.rs:1750-1762`, both bail on the *first* BAM | read both; ran cases A/B/C | **Confirmed.** Both loop and emit the two-BAM message for a pair containing one BAM, in either argument order |
| Those are the *only* such sites | `grep -rn "two BAM files" src/ tests/` → `main.rs:543`, `:1271`, `:1754` only | **Confirmed.** No fourth site |
| `main.rs:537-557` is already correct, not a cause site | read `:536-558`; ran case F | **Confirmed.** Mixed pair yields *"requires both input files to be the same format. Got mixed: A and B."* The issue mis-attributes it |
| Guard placement `:264` precedes every side effect | traced `:264` → `:272` gzip → `:273` output_dir → `:277 ensure_output_dir` → `~:321` specialty → `:409` clump-only → `:654`/`:661`/`:675` trim | **Confirmed.** Nothing is created or opened between `:264` and any dispatch |
| A2: `Cli::validate()` guarantees even N (`cli.rs:506`), R1≠R2 (`:513`), N=1 carve-out (`:503`) | read; confirmed `--paired` always routes through `validate_paired_input("Paired-end")` at `cli.rs:593-597` | **Confirmed.** `chunks(2)` cannot see a short chunk |
| §3.1 field types (`hardtrim5/3: Option<usize>`, `clock: bool`, `implicon: Option<usize>`, `clump_only: bool`) | `grep -n "pub hardtrim5\|…" src/cli.rs` → `:372 :377 :383 :391 :208` | **Confirmed.** Expression is type-correct |
| Specialty exemption cannot accidentally skip clump-only | `cli.rs:805-814` makes `--clump_only` mutually exclusive with all four specialty modes | **Confirmed**, and dispatch order (`:321` … `:409`) makes it doubly safe |
| Existing three tests' pinned substrings survive §3.3's wording | read `integration_clump_only_ubam.rs:400-403`, `:618-621`; `integration_ubam_out.rs:461-465` | **Confirmed** — subject to I-5 below |
| §2.5 specialty behaviour (hardtrim silently accepts; clock fails with a read-count error) | ran cases G/H | **Confirmed verbatim** |
| §2.4 / §3.3 line references, and the 421-test count in V10 | read each; `grep -rc "#\[test\]"` → 421 | **All accurate** |
| No paired entry point missed | `--demux` + `--paired` rejected at `cli.rs:923-925`; `--passthrough` + any BAM rejected at `main.rs:251-255`; `--retain_unpaired` rides `run_paired`; `run_paired`/`run_ubam_output` each have exactly one caller (`:761`, `:654`) | **Confirmed.** §11's enumeration is complete |
| Deleting all three guards is reachability-safe | activation set of each deleted guard is a strict subset of the hoisted guard's (`run_paired` needs `paired ∧ N>1 ∧ Fastq-out ∧ ¬specialty ∧ ¬clump_only`; `run_ubam_output` paired needs the same with `UBam`; clump-only Shape A needs `clump_only ∧ paired ∧ N>1 ∧ UBam`) | **Confirmed** |

Two smaller correctness points I checked because they could have bitten silently and did not: `input_formats` at `:182-186` really is computed for *every* input (so the helper needs no I/O, and in fact **removes** the two redundant `detect_input_format` calls per pair that the deleted sites made — the plan understates its own win); and flipping the guard ahead of the collision pre-flight does not disturb `pe_bam_collision_preflight_case_folded`, whose inputs are all `.fq.gz` copies, so `n_bam = 0` and the collision is still what gets reported.

---

## 1. Logic review

### C-1 (Critical) — the hoisted guard preempts `--clump_only` with FASTQ output and replaces a correct diagnosis with a command that cannot work

§3.1 states `cli.clump_only` is **not** excluded from the guard. §3.3 keys `fmt_flag` on whether uBAM output was requested. Put together, `--clump_only --paired` with **default FASTQ output** now falls into the new messages with `fmt_flag = ""`.

What that path does **today** (ran both; `dev` @ `a4ffd47`):

```console
$ trim_galore --clump_only --paired BS-seq_10K_R1.fastq.gz ubam_test.bam -o D
Error: processing --clump_only pair 1 of 1 (…)
Caused by:
    uBAM input under --clump_only requires --output-format ubam (using the FASTQ
    output path with uBAM input would drop aux tags). Input: …/ubam_test.bam

$ trim_galore --clump_only --paired ubam_test.bam ubam_copy.bam -o E
    … same message …
```

That message is correct, specific, and immediately actionable — it names the flag to add. It comes from `clump_only.rs:401-410`.

What the plan's guard would say instead, for the mixed case:

> `--clump_only --paired` requires both inputs of a pair to be the same format. Pair 1 of 1 is mixed: … . **Pass two FASTQ files, or a single interleaved uBAM.** If you meant two FASTQ files, check for a mis-typed filename.

and for the two-BAM case:

> … uBAM paired mode expects a single interleaved file: **`trim_galore --clump_only --paired interleaved.bam`** …

Both remediations are wrong, and I verified the second one is *hard* wrong:

```console
$ trim_galore --clump_only --paired ubam_paired_test.bam -o remed2
Error: --clump_only --paired requires two FASTQ input files. Single-file
interleaved uBAM input needs --output-format ubam (add `--output-format ubam`
to the command line).
```

That is `main.rs:426-432`. So the plan's message hands a clump-only user a pasteable command whose only effect is a second error, and it *drops* the one piece of information they needed (`--output-format ubam`). This is precisely the failure mode #363 is about — the remediation pointing the wrong way — recreated in a different mode.

Two aggravating factors:

- **No test would catch it.** `tests/integration_clump_only.rs:303 rejects_ubam_input_without_ubam_output` is **single-end** (`--clump_only <bam>`, no `--paired`), so it exercises `clump_only.rs:265`, not `:402`. `grep -rn "requires --output-format ubam" tests/` returns nothing else. The paired FASTQ-output rejection has zero coverage, and the plan's V-table (V5 covers only `--clump_only --paired --output-format ubam`) does not add any.
- **§11 asserts the opposite.** It records having traced "clump-only FASTQ (`:426`, which has its own N=1 rejection and its own 'uBAM input requires --output-format ubam' error — confirmed live, case D, and unaffected)". Only the **N=1** half is unaffected; for N≥2 the hoisted guard fires first and `clump_only.rs:401-410` becomes unreachable. The plan's live check appears to have exercised the N=1 shape and generalised.

**Recommended fix.** Replace the two opaque `&str` parameters with a mode enum so the remediation is exact per path, e.g.

```rust
enum PairedShape { Trim, TrimUbamOut, ClumpOnlyFastqOut, ClumpOnlyUbamOut }
```

`Trim` and `TrimUbamOut` keep §3.3's wording (I confirmed `trim_galore --paired interleaved.bam` genuinely works — it routes through `run_paired_ubam_single_file` and produced `ubam_paired_test_val_{1,2}.fq` plus both reports). `ClumpOnlyUbamOut` keeps §3.3's wording with the flag appended. `ClumpOnlyFastqOut` must **preserve today's diagnosis** — "uBAM input under `--clump_only` requires `--output-format ubam` (the FASTQ output path would drop aux tags)" — because for that mode the fix is a flag, not a re-shaped input. The minimum viable variant of this fix is to derive `fmt_flag` as `" --output-format ubam"` whenever `cli.clump_only`, but that still loses the aux-tag reason, which is the *why* a user needs.

### I-2 (Important) — the deletions convert three independent guards into a single point of failure in front of a silent-wrong-output path

The plan's risk register rates "an unenumerated paired entry point still carries a stale guard" as *very low*. The inverse risk — the hoisted guard being weakened later — is not registered, and its failure mode is not an error but silently wrong output. Two leaf functions read the source header from **R1 only** and open each side by per-file detection:

- `main.rs:1991-2001` — `run_ubam_output_paired_two_files`: `detect_input_format(input_r1)` decides the `@PG`/`@HD` lineage, then `open_sync_reader` is called on each path independently. Its comment at `:1993-1994` literally says *"At this point both inputs are FASTQ"* — an invariant that, after this change, is maintained by exactly one site 1,700 lines away.
- `clump_only.rs:946-958` — Shape A: same shape, `fmt` from R1, `open_sync_reader` per side.

A mixed pair reaching either would emit a BAM mixing FASTQ-derived records (no aux, no source header) with BAM-derived ones, with no error. Any future addition to §3.1's exemption list — and §3.1 already has four exemptions plus a follow-up issue explicitly contemplating revisiting them — reopens that door.

Cheapest mitigation that keeps the plan's line-count win: at each leaf, one `debug_assert!` (or a `bail!` guarded by a comment naming the single upstream guard) asserting R1's and R2's BAM-ness agree. `run_ubam_output_paired_two_files` already calls `detect_input_format(input_r1)`; adding R2 costs one 24-byte peek per pair on an I/O-bound path.

At an absolute minimum, **both stale comments must be updated**, not just the two the plan lists. §4 step 3/4 names `main.rs:1266-1267` and `clump_only.rs:940-944`. It misses `main.rs:1309-1310` (*"paired-BAM rejected above, so both are FastqReader"* — "above" ceases to exist) and `main.rs:1993-1994`.

### I-3 (Important) — the `samtools collate` hint is technically wrong for the two-BAM case

§3.3 / Q3 justify the hint on the grounds that "the phrasing exists" in `bam.rs`. It does — at `bam.rs:52-58`, `GROUPED_INPUT_ERR`:

```
Re-interleave the input first: `samtools collate -O input.bam tmp > interleaved.bam`
```

but that fires on **one** file whose mates are non-adjacent, where `collate` is exactly right. The two-BAM branch has **two** files. `samtools collate` takes a single input; it cannot interleave `r1.bam` and `r2.bam`. The user needs a merge or concatenate step first, e.g.

```
samtools merge -n -o - r1.bam r2.bam | samtools collate -O - tmp > interleaved.bam
```

Copying the hint across without adjusting for the different situation is the exact class of error this plan is fixing elsewhere. Either give the two-step command or drop the hint and say "combine the two BAMs into one mate-adjacent file".

### I-4 (Important) — an un-enumerated user-visible behaviour change: multi-pair runs stop producing partial output

§7 lists three visible changes. There is a fourth, and it is the largest. Today the trim FASTQ path's guard is *inside* the per-pair loop, so earlier pairs complete before the offending pair fails. Verified (N=4, offence in pair 2):

```console
$ trim_galore --paired BS-seq_10K_R1.fastq.gz BS-seq_10K_R2.fastq.gz \
              phred64_test.fastq ubam_test.bam -o I
Error: processing pair 2 of 2 …
$ ls I
BS-seq_10K_R1_val_1.fq.gz   BS-seq_10K_R1.fastq.gz_trimming_report.{txt,json}
BS-seq_10K_R2_val_2.fq.gz   BS-seq_10K_R2.fastq.gz_trimming_report.{txt,json}
```

After the hoist, that run produces **nothing**. I think fail-fast is the right call — failing in a second beats failing after hours on pair 97 of 100 — but it changes what a pipeline sees, it is more consequential than the loss of the `processing pair N of M` wrapper that §7 does list, and it belongs in the CHANGELOG. It also means **V8 is under-specified**: it checks the pair index in the message but not the delta that matters. V8 should additionally assert that pair 1's `*_val_1.fq*` and `*_trimming_report.txt` are absent.

### I-1 (Important) — the helper's name and doc contradict §3.2, and the obvious wrong implementation breaks a working invocation

§3.2 is right: decide on `n_bam ∈ {0,1,2}`. But the name `reject_mismatched_pair_formats` and the doc's opening line — *"Reject paired inputs whose two members do not share an input format"* — describe a different predicate, and the natural transcription of that sentence is `formats[0] != formats[1]`. `InputFormat` has **three** variants (`format.rs:25-34`: `FastqPlain`, `FastqGz`, `UnalignedBam`), so that predicate rejects a plain+gzipped FASTQ pair. Verified that pair works today:

```console
$ gunzip -c BS-seq_10K_R1.fastq.gz > plainR1.fastq
$ trim_galore --paired plainR1.fastq BS-seq_10K_R2.fastq.gz -o mixcomp
exit=0   →  plainR1_val_1.fq  BS-seq_10K_R2_val_2.fq  + both reports
```

Rename to something predicate-accurate (`reject_bam_format_mismatch_in_pair`, or `check_pair_bam_shape`), restate the doc in terms of BAM count, and add a positive integration test pinning plain+gz acceptance. Note also that `InputFormat` has no `Display` impl (`format.rs:24` derives only `Debug, Clone, Copy, PartialEq, Eq`), so §3.3's `{fmt1}`/`{fmt2}` needs a local label function — `clump_only.rs:689` already does `UnalignedBam => "uBAM"` and is worth matching for consistency. Whatever it renders, plain and gzipped FASTQ must not read as two different things in a message that is about them being different.

### I-5 (Important) — a concrete way V9 fails: `\`-continuation splitting a pinned substring

`rejects_two_bam_paired` pins `"single interleaved"`. §3.3 renders the message with its line break falling **exactly between "single" and "interleaved"**. Transcribed into a Rust multi-line literal, `"… expects a single\` + newline + indented `interleaved file"` yields `"a singleinterleaved file"`, because `\`+newline consumes the newline *and* all leading whitespace on the next line. The existing code gets this right by ending each line with a space before the backslash (`main.rs:1271-1274`), but the plan's rendering invites the wrong break. Reflow §3.3 so none of the four pinned substrings — `two BAM files is not supported`, `single interleaved`, `uBAM paired mode expects`, `same format` — straddles a source-line break, and keep V9 as the backstop.

### I-6 (Important) — validation gaps

V3, V9, V11 are well chosen and the plan is right to flag them as the ones that quietly stop holding. But four things the plan cares about are checked only by eye:

1. **Output-dir absence (V3) is manual.** It is the §1 outcome most easily lost in a refactor and there is no test for it. Make it an integration assertion. Gotcha for whoever writes it: the test helper `fresh_tmpdir` (`integration_ubam_out.rs:28-32`) does `create_dir_all`, so passing it to `-o` guarantees the directory exists regardless. The test must pass `-o <fresh_tmpdir>/nested` and assert `!nested.exists()`.
2. **The substance of #363 has no negative assertion.** V1 says "**no** bare interleaved-file imperative" but nothing pins it. Add `assert!(!stderr.contains("uBAM paired mode expects"))` to the mixed-pair test — that is the one-line regression guard that stops the two messages from re-converging in a future refactor.
3. **C-1's blind spot:** `--clump_only --paired` with FASTQ output, mixed and two-BAM. Both currently untested and both change behaviour under this plan.
4. **I-1's regression:** plain+gz FASTQ pair still accepted.

---

## 2. Assumptions

- **A1** — confirmed against `format.rs:36-83` and its `detect_bgzipped_fastq_is_fastq_not_bam` test at `:176`. Sound.
- **A2** — confirmed exactly, line for line. One thing worth adding to the helper: the doc states the even-count and R1≠R2 preconditions but not `inputs.len() == formats.len()`, which is what its own indexing actually depends on. A `debug_assert_eq!` costs nothing.
- **A3** — confirmed (`main.rs:257-263`).
- **A4** — the *reasoning* is fair (hardtrim genuinely treats each file independently, and I confirmed `--hardtrim5 10 --paired <fq> <bam>` writes both outputs and exits 0). But the assumption is stated as covering "specialty modes" as a category, and case H shows `--clock` is *not* harmless: it fails with `Paired-end files have different numbers of reads!` and, as the plan itself notes in §2.5, would have **proceeded** had the counts matched. So the exemption is defensible for `--hardtrim5/3` and indefensible-but-deferred for `--clock`/`--implicon`. That is a scope decision the user already took, and I am not reopening it — but A4's justification should say so plainly rather than resting on "harmless in hardtrim", because two of the four exempted modes are not harmless. The follow-up issue (§4 step 9) is doing real work here and should not be optional.
- **A5** — confirmed by reading all three assertions; see I-5 for the way it can still fail.
- **A6** — sound, and stronger than the plan claims. `main.rs:1991-2001` and `clump_only.rs:946-958` both take the source header from R1 alone, so "a mixed pair *could* be processed" would in practice mean a BAM with a header from one arbitrary side. Rejection is not merely the likelier-intent call; the accepting path does not exist in a correct form. Worth saying, because it is also the argument for I-2.

Two implicit assumptions the plan does not surface:

- **That `sanity_check_any(&cli.input[0])` running at `:176` — before the guard — is acceptable.** It is (a read, not a write), but it means a mixed pair whose R1 is an *empty* BAM still reports the empty-BAM error rather than the format error. Fine; just not the "guard fires first" story §3.4 tells.
- **That "no banner" means no *trimming* banner.** §3.4 says "no banner"; the version banner at `main.rs:172-174` prints unconditionally on every run, rejected or not. V3's wording is correct; §3.4's is loose.

---

## 3. Efficiency

Nothing to dispute — and the plan is slightly pessimistic about itself.

- One pass over an in-memory `Vec<InputFormat>` of ≤ N elements. Immeasurable.
- On rejected runs it removes the full adapter-detection scan (up to 1 M reads, `adapter.rs`) plus poly-G detection plus `ensure_output_dir`. Verified from case A's stderr that all three currently run before the error.
- Not claimed but true: the deleted sites each re-called `detect_input_format` (`main.rs:1269`, `:1752`, `:537-538`) — 2 extra file opens and up to 2 BGZF-block decompresses **per pair**. The helper takes pre-computed formats, so those disappear on the accepted path too, not just the rejected one.
- Scalability of the guard itself: O(N/2) with an early exit. The only scale-relevant consequence is I-4 (all-or-nothing across pairs), which is a semantics change, not a cost one.

---

## 4. Alternatives

1. **Put the helper in `src/format.rs`, not `main.rs`.** Q2 has a definite answer I can supply: `src/main.rs` contains **no** `#[cfg(test)]` module (`grep -n "cfg(test)" src/main.rs` → nothing). Both of the plan's options work — `Cargo.toml` does not set `test = false` on the `[[bin]]`, so a test module added inside `main.rs` would run under `cargo test` — but `format.rs` is the better home on the merits: it owns `InputFormat`, it already has a six-test module including the load-bearing bgzip-vs-BAM case, and a pure decision function over `&[InputFormat]` is exactly its remit. It also removes the "library placement as a fallback" branch from Step 5.
2. **`Cli::validate_formats(&self, formats: &[InputFormat])` in `cli.rs`.** The plan's rationale for `main.rs` — "validate() runs before format detection" — is about *timing*, not *location*. A second, format-aware validate method called from `main.rs:264` satisfies both, and would collect what are now four inline format-dependent guards (`--phred64` `:211`, `--preserve-tags` `:233`, `--passthrough` `:251`, `--paired` N=1 `:257`, plus this one) into the module a maintainer looks in for validation. `main.rs` is 2,400 lines. Not required for this fix; worth a line in §10 as the direction of travel.
3. **Mode enum instead of two `&str`s** — see C-1. This is not merely tidier; it is what makes the clump-only-FASTQ-output remediation expressible at all. It also removes the risk of a caller passing a `mode` string that does not match the actual flags.
4. **Fold in `main.rs:498-512`.** That clump-only Shape B "must be BAM" check is *already* unreachable today: its condition (`clump_only ∧ UBam ∧ paired ∧ N=1 ∧ ¬BAM`) is a strict subset of `main.rs:257-263`'s. Same family of redundancy the plan is retiring, one branch away from the code it is already editing. Optional, but cheap while the file is open.

---

## 5. Action items

### Critical

- **C-1.** Fix the clump-only-with-FASTQ-output remediation. `--clump_only --paired <fq> <bam>` and `<bam> <bam>` (no `--output-format ubam`) must keep today's *"uBAM input under `--clump_only` requires `--output-format ubam` (the FASTQ output path would drop aux tags)"* diagnosis, not be told to pass a single interleaved uBAM — a command `main.rs:426-432` rejects. Adopt a mode enum (§4 alt 3) rather than `mode`/`fmt_flag` strings. Correct §11's "unaffected" claim: only the N=1 clump-only-FASTQ shape is unaffected. Add integration coverage for both shapes (currently zero).

### Important

- **I-1.** Rename the helper and restate its doc in terms of BAM count, not format equality; add a test pinning that a `FastqPlain` + `FastqGz` pair is still accepted (verified working today). Add a format-label helper — `InputFormat` has no `Display`.
- **I-2.** Update `main.rs:1309-1310` and `main.rs:1993-1994` as well as the two comments §4 already lists, and keep a one-line `debug_assert!`/`bail!` at `run_ubam_output_paired_two_files` and `clump_only_paired_to_bam_one_pair` Shape A. Both take the source header from R1 alone, so a mixed pair slipping through is silent wrong output, not an error.
- **I-3.** Fix the two-BAM `samtools` hint: `collate` alone cannot interleave two separate files. Give the merge-then-collate pipeline, or drop the command and describe the goal.
- **I-4.** Add the fourth behaviour change to §7 and the CHANGELOG: multi-pair runs no longer produce output for pairs preceding the offending one (verified). Strengthen V8 to assert pair 1's outputs and reports are absent.
- **I-5.** Reflow §3.3 so no pinned substring straddles a source-line break (`"a single" / "interleaved file"` is the live hazard, via Rust `\`-continuation whitespace stripping).
- **I-6.** Promote V3 to an automated assertion (`-o <existing>/nested`, then `!nested.exists()` — `fresh_tmpdir` pre-creates its own directory). Add `assert!(!stderr.contains("uBAM paired mode expects"))` to the mixed-pair test as the permanent #363 guard.

### Optional

- **O-1.** §2.1 generalises one site's stderr to both. On `run_ubam_output` (`:1750`) the error fires *before* adapter detection — only `ensure_output_dir` and the uBAM NOTEs precede it (verified, case C). Also, §3.4's "no banner" should read "no *trimming* banner"; `main.rs:172-174` always prints.
- **O-2.** Prefer `src/format.rs` for the helper and close Q2: `main.rs` has no test module, and `format.rs` owns `InputFormat` plus an existing test module.
- **O-3.** `main.rs:498-512` is already redundant with `:257-263`; retire it in the same pass.
- **O-4.** §4 step 8: `### Unreleased` already has a `#### Fixes` subsection (`CHANGELOG.md:84`) — append to it rather than creating a second one.
- **O-5.** "Pair 1 of 1" is noise in the single-pair case (the overwhelmingly common one). Consider omitting the index when `n == 1`; house style at `main.rs:772-780` does print it, so either choice is defensible.
- **O-6.** Add `debug_assert_eq!(inputs.len(), formats.len())` and state that precondition in the doc — it is what the indexing actually relies on.
- **O-7.** Consider `Cli::validate_formats(&self, &[InputFormat])` as the eventual home for all five format-dependent guards (§4 alt 2). Out of scope here; worth a §10 note.

---

## 6. Evidence log

Built `cargo build --release` from `/Users/fkrueger/Github/TrimGalore` at `a4ffd47`; all output to `$TMPDIR`.

| Case | Command | Observed |
|---|---|---|
| A | `--paired phred64_test.fastq ubam_test.bam` | adapter + poly-G detection ran, banner printed, dir created, then two-BAM message. §2.1 reproduced verbatim |
| B | `--paired ubam_test.bam phred64_test.fastq` | same, names the BAM |
| C | `--paired --output-format ubam <fq> <bam>` | two-BAM message **before** adapter detection; dir created (O-1) |
| D | `--clump_only --paired <fq> <bam>` (FASTQ out) | *"uBAM input under --clump_only requires --output-format ubam"* — correct today (C-1) |
| E | `--clump_only --paired <bam> <bam-copy>` (FASTQ out) | same message (C-1) |
| F | `--clump_only --paired --output-format ubam <fq> <bam>` | *"requires both input files to be the same format. Got mixed: …"* — §2.2 confirmed |
| G | `--hardtrim5 10 --paired <fq> <bam>` | exit 0, both outputs written — §2.5 confirmed |
| H | `--clock --paired clock_R1.fq.gz ubam_test.bam` | `Paired-end files have different numbers of reads!` — §2.5 confirmed |
| I | `--paired <fq> <fq> <fq> <bam>` (N=4) | pair 2 rejected **after** pair 1's outputs + reports were written (I-4) |
| remed-1 | `--paired ubam_paired_test.bam` | succeeds — the trim-path remediation is valid |
| remed-2 | `--clump_only --paired ubam_paired_test.bam` | **rejected** — the clump-only remediation is not (C-1) |
| mixcomp | `--paired plainR1.fastq BS-seq_10K_R2.fastq.gz` | exit 0 — plain+gz pair works today (I-1) |

Distinguishing confirmation from suspicion: C-1, I-1, I-3, I-4, I-5, O-1, O-2, O-3, O-4 are **confirmed** by reading the cited code and/or running the cited command. I-2 is a **confirmed structural fact** (both leaf functions do read the header from R1 only, and both comments will be stale) with an **inferred** likelihood of future breakage. I-6 is a confirmed absence of coverage. Nothing in this report is speculative about current behaviour.
