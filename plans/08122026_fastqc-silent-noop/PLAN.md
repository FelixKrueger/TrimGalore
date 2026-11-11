# PLAN — #421: make `--fastqc` either work or say why not

**Issue:** [#421](https://github.com/FelixKrueger/TrimGalore/issues/421)
**Branch base:** `dev` @ `3ffcef4`
**Type:** behaviour change — a currently-succeeding flag combination starts failing, plus one
silently-ignored flag starts working
**Decision taken before drafting (Felix, 2026-08-12):** **refuse the combination at the CLI** on
the arms that cannot run FastQC, rather than wiring FastQC into them or warning and continuing.
**Decisions taken after dual plan review (2026-08-13):** the **sixth arm** both reviewers found
is fixed in this PR rather than deferred to a new issue (§2.3); the guard block moves **late** in
`validate` (§3); arm 5's message names the input *shape* rather than asserting a format (§5.1).

---

## 1. Goal

`--fastqc` is accepted and silently ignored on five dispatch arms: the run exits 0, writes its
normal output, produces no QC report, and never mentions FastQC. Dual review found a **sixth**
case of the same defect from the other direction — `--fastqc_args` without `--fastqc` is silently
ignored on all four `--clump_only` drivers, which are otherwise capable.

Refuse the five, wire the sixth, so the flag either works or says why not.

---

## 2. Context

### 2.1 Where FastQC is reachable, and under which gate

Twelve `fastqc::run` call sites. **The gate is not uniform** — that is the finding that reshaped
this plan:

| Site | Arm | Gate |
|---|---|---|
| `main.rs:1459` | SE trim | `cli.fastqc \|\| cli.fastqc_args.is_some()` |
| `main.rs:1778`, `:1784`, `:1796` | PE trim, incl. `--passthrough` carrier | same disjunction |
| `main.rs:2242` | SE uBAM output | same disjunction |
| `main.rs:2399` | PE uBAM output, FASTQ pair in | same disjunction |
| `main.rs:2522` | PE uBAM output, **single interleaved BAM in** | same disjunction |
| `clump_only.rs:374` | `clump_only_single` | **bare `if fastqc`** |
| `clump_only.rs:554`, `:555` | `clump_only_paired` (R1 and R2) | **bare `if fastqc`** |
| `clump_only.rs:847` | `clump_only_single_to_bam` | **bare `if fastqc`** |
| `clump_only.rs:1095` | `clump_only_paired_to_bam_one_pair` | **bare `if fastqc`** |

So the disjunction covers 7 of 12 calls (5 gate expressions in `main.rs`); the other 5 read a
`fastqc: bool` parameter that `main.rs` feeds with bare `cli.fastqc` at `:636`, `:671`, `:715`,
`:778`, `:819`. Verified by reading all five sites.

`src/specialty.rs` contains zero `fastqc` references, and `run_paired_ubam_single_file`
(`main.rs:1810-2003`) contains none. Those absences are arms 1–5.

### 2.2 The five arms to refuse

`main()`'s dispatch early-returns in a fixed order; both reviewers independently derived the same
table. Arms 1–4 are guarded by flag presence, so arm 5's predicate needs only the terms that arms
1–4 do not already catch:

| | Arm | Dispatch | CLI-visible condition | FastQC |
|---|---|---|---|---|
| 1 | `--hardtrim5` | `main.rs:464` | `hardtrim5.is_some()` | no |
| 2 | `--hardtrim3` | `:495` | `hardtrim3.is_some()` | no |
| 3 | `--clock` | `:526` | `clock` | no |
| 4 | `--implicon` | `:548` | `implicon.is_some()` | no |
| — | `--clump_only` | `:580` | — | yes (see §2.3) |
| — | uBAM output | `:845` | — | yes ×3 |
| 5 | paired **FASTQ** from one interleaved uBAM | `:853` | `paired && input.len() == 1 && output_format != UBam && !clump_only` | no |
| — | PE trim | `:867` | — | yes |
| — | SE trim | loop | — | yes |

**Arm 5's `output_format` term is load-bearing.** `main.rs:845` sits *above* `:853`, so
`--paired --output-format ubam <one.bam>` reaches `run_ubam_output_paired_single_file`, which runs
FastQC at `:2522`. Dropping the term refuses a working feature. It is guarded locally as well as in
CI: `tests/integration_ubam_out.rs:572` runs that exact invocation under plain `cargo test` and
asserts one `_fastqc.zip`, and `.github/workflows/ci.yml:1173` asserts the same in CI.

The `!clump_only` term is about message precedence, not capability: `main.rs:597-602` already
refuses `--clump_only --paired` with a single input and points at `--output-format ubam`, which is
a better diagnosis than a FastQC complaint.

Both reviewers hunted for an accepted-today invocation these five conditions would refuse and
found none, and for a silently-ignoring invocation in the `paired && len == 1` family that they
would miss, and found none.

### 2.3 The sixth arm — the same defect, from the other side

`--clump_only --fastqc_args "--quiet" a.fastq` exits 0, produces no report and says nothing.
Reviewer A reproduced it; the mechanism is §2.1's non-uniform gate. This is exactly the failure
this plan's own reasoning predicts — "a refusal keyed on `cli.fastqc` alone would leave
`--fastqc_args` silently ignored, the same defect one flag over" — applied to code the first draft
declared safe on the strength of a grep for `fastqc::run` rather than for its guard.

**Decision: fix it here, not in a new issue.** The arm is *capable*; it needs one flag wired
through, so refusing it would remove a working feature. The fix is `Cli::fastqc_requested()` plus
five one-word substitutions, and it closes #421's stated goal on all six arms rather than shipping
a `--fastqc_args`-shaped hole in the same PR that refuses `--fastqc_args` elsewhere. It also keeps
two in-repo promises honest: `cli.rs:333` ("`--fastqc_args` … Implies --fastqc") and
`docs/.../modes/clump-only.md:77` ("`--fastqc` — runs on the reordered output").

### 2.4 The house pattern, and where the block belongs

`cli.rs:586-634` (§3.4a) is the model for the guards: a flat sequence of
`if cond { anyhow::bail!(…) }`, each naming flag, mode and remedy.

**But the block must be sited late in `validate`, not after §3.4a.** Reviewer B found that
several per-flag checks live below line 634 and would be pre-empted:

| Invocation | Message today | If sited at 634 |
|---|---|---|
| `--hardtrim5 0 --fastqc a.fq` | `--hardtrim5 must be between 1 and 999, got 0` (`cli.rs:996`) | FastQC refusal |
| `--hardtrim5 20 --hardtrim3 20 --fastqc a.fq` | cannot combine the two | FastQC refusal |
| `--clock --fastqc r1 r2 r3` | `--clock mode requires an even number of input files` (`cli.rs:1015`) | FastQC refusal |
| `--implicon --fastqc r1.fq r1.fq` | `Read 1 and Read 2 appear to be the same file` | FastQC refusal |

**Both of those rows said `--paired` in v2 and cited `:1004`/`:1007`, and both were wrong.** With
`--paired` set, `validate_paired_input("Paired-end")` at `cli.rs:650` fires first, so `--clock`'s and
`--implicon`'s own calls at `:1015`/`:1018` are reachable **only without** `--paired`. Reviewer B
proved it by mutation: the v2 test used `--paired` and stayed green when the whole guard block was
relocated above `:997`. Corrected here and in the test.

Verified: the hardtrim range checks are at `cli.rs:996/1001`, `validate_paired_input` runs for
`--clock` at `:1004` and `--implicon` at `:1007`, and `check_restartable_input` at `:1020`. A
malformed *value* or a structurally invalid input list should be reported before an unsupported
*combination*, which is the argument the codebase already makes at `main.rs:344-351` ("Sited here
… so those report their own defect first"). **Site the block after the last fallible check in
`validate`** — which is the `--adapter2` parse, not the `check_restartable_input` loop; v2 said
the loop and that was one check short (see Implementation notes, iteration #5). Cost: zero.

### 2.5 #422's tripwire becomes the second half of this change

`tests/integration_preflight_tripwire.rs` carries a `fastqc_capable` field, a
`.fastqc_incapable()` builder marking eight cases, and a match arm whose message reads
*"#421 may have been fixed — drop the marker"*.

Both reviewers verified the deletion is safe and that it is safe for a stronger reason than v1
claimed: the `(true, false)` branch is not merely unfired but **unexercisable**, because every
direct `compare()` caller passes `fastqc_capable: true`. The eight marked cases map onto exactly
arms 1–5, and collapsing the match to `if active && allowed.is_empty()` *widens* the surviving
rule to all 53 cases rather than weakening it.

Two corrections to v1's plan for it:

- **The deletion inventory was short.** `fastqc_capable` is also `compare()`'s 6th positional
  parameter (`:348`), passed at `Tripwire::run` (`:311`), with nine literal `true` arguments in the
  `comparison_*` and helper callers. Under `-D warnings` a half-done deletion fails on an unused
  binding, so the list has to be complete.
- **Deleting a guard without leaving one behind loses the only in-repo record that these arms
  cannot run FastQC.** Replace it with five *active* cases — non-zero exit, the message, and
  nothing written — following `a_refused_run_writes_nothing` (`:1350`), the file's existing idiom.
  That converts a guard that never fired into one that fires on every run.

### 2.6 Nothing in the repo currently consumes the refused combinations

Both reviewers swept `tests/`, `.github/workflows/`, `docs/`, `scripts/`, `examples/` and
`README.md`: no test, CI step or doc pairs a FastQC flag with any of arms 1–5. Recorded here so it
is not re-derived.

Blast radius outside the repo is small: nf-core's `trimgalore` module takes `--fastqc` only from
`task.ext.args` (it is not in the module command) and nf-core/rnaseq runs a standalone `FASTQC`
module rather than the flag, so the flagship consumer does not break. The residual risk is a
wrapper that sets `--fastqc` globally while exposing specialty modes — that wrapper goes from
silently QC-less to hard-failing.

---

## 3. Behavior

1. Add `Cli::fastqc_requested(&self) -> bool` returning `self.fastqc || self.fastqc_args.is_some()`.
2. In `Cli::validate`, **after the `--adapter2` parse** — the last fallible check (§2.4) — if
   `self.fastqc_requested()`, refuse for each of arms 1–5 (§2.2), one `bail!` per arm.
3. At `main.rs:636`, `:671`, `:715`, `:778`, `:819`, pass `cli.fastqc_requested()` instead of
   `cli.fastqc`, so `--fastqc_args` alone activates FastQC on the `--clump_only` drivers.
4. Nothing else changes. No arm gains or loses a `fastqc::run` call site.

### Edge cases

| Input | Behaviour |
|---|---|
| `--hardtrim5 20 --fastqc a.fq` | error naming `--hardtrim5` |
| `--hardtrim5 20 --fastqc_args "--quiet" a.fq` (no `--fastqc`) | error — the disjunction |
| `--clock --paired --fastqc r1 r2` | error naming `--clock` |
| `--implicon --paired --fastqc r1 r2` | error naming `--implicon` |
| `--paired --fastqc one.bam` | error naming the `--paired`-with-one-input shape |
| **`--paired --output-format ubam --fastqc one.bam`** | **accepted** — capable, pinned locally and in CI |
| `--paired --output-format ubam --fastqc r1.fq r2.fq` | accepted — capable (`main.rs:2399`) |
| `--output-format ubam --fastqc one.bam` (SE) | accepted — capable (`main.rs:2242`) |
| `--clump_only --fastqc a.fq` | accepted, report produced — unchanged |
| **`--clump_only --fastqc_args "--quiet" a.fq`** | **accepted, and now produces a report** — §2.3. Red today |
| `--clump_only --paired --fastqc one.bam` | unchanged: `main.rs:597`'s message still wins |
| `--hardtrim5 0 --fastqc a.fq` | **range error still wins** — §2.4 |
| `--clock --paired --fastqc r1 r2 r3` | **odd-count error still wins** — §2.4 |
| `--paired --fastqc single.fastq` | FastQC refusal, whose wording must also be true for this user — §5.1 |
| `--hardtrim5 20 a.fq` (no FastQC flag) | unchanged |

---

## 4. Signature

```rust
/// Whether the run asked for FastQC. `--fastqc_args`' "Implies --fastqc" is doc text, not a
/// clap `requires`, so `fastqc` alone is false on an args-only run.
pub fn fastqc_requested(&self) -> bool {
    self.fastqc || self.fastqc_args.is_some()
}
```

One accessor, on `Cli`. It replaces a predicate currently spelled out at five `main.rs` gate
expressions and derived a second time in the tripwire (`fastqc_active`, `:223`, whose own comment
argues for one source of truth). No message helper: §3.4a's flat `if cond { bail!() }` sequence is
the file's pattern eight times over, and five arms do not justify departing from it.

---

## 5. Implementation outline

1. **`src/cli.rs`** — add `fastqc_requested()`; add the five guards after the `--adapter2`
   parse, the last fallible check in `validate` (§2.4). Message shape:

   > `--fastqc is not supported with --hardtrim5, which writes no QC report; run FastQC
   > separately on the hardtrimmed output, or drop --fastqc.`

   **Arm 5's message must not assert the input is a uBAM.** The reachable-but-wrong case is a user
   who forgot R2, which is likely the *more* common mistake; v1's wording told them their FASTQ was
   an interleaved uBAM and offered `--output-format ubam`, which still fails for them at
   `main.rs:334`. Name the shape instead:

   > `--fastqc is not supported with --paired and a single input file: that shape is the
   > interleaved-uBAM paired path, which writes FASTQ output with no QC report. Use
   > --output-format ubam (which does produce one), pass R1 and R2 as two files, or run FastQC
   > separately.`

2. **`src/main.rs`** — five one-word substitutions: `cli.fastqc` → `cli.fastqc_requested()` at
   `:636`, `:671`, `:715`, `:778`, `:819`. No `clump_only.rs` signature changes.

3. **`src/cli.rs` unit tests** — one per refused arm, following the
   `output_format_ubam_plus_*_rejected` naming. Plus the acceptance twins that catch a
   mis-specified condition: `--paired --output-format ubam --fastqc one.bam`,
   `--paired --output-format ubam --fastqc r1 r2`, `--clump_only --fastqc a.fq`. Plus one
   `--fastqc_args`-only refusal, and two precedence tests (`--hardtrim5 0 --fastqc` → range error;
   `--clock --paired --fastqc r1 r2 r3` → odd-count error).

4. **`tests/integration_preflight_tripwire.rs`** — delete `fastqc_capable`, `.fastqc_incapable()`
   and its eight call sites, `compare()`'s 6th parameter (`:348`), the argument at `:311`, and the
   nine literal `true`s in the `comparison_*` callers; collapse the match to
   `if active && allowed.is_empty()`. Rewrite the module-doc paragraph (`:25-36`). **Add five
   refusal cases** — non-zero exit, message, nothing written — per `a_refused_run_writes_nothing`
   (`:1350`).

5. **Integration tests** — in `tests/integration_paired_format_guard.rs`, **not**
   `integration_no_args_help.rs` (which is narrowly about the bare-invocation help path). That file
   already owns `--paired` single-file refusals, has a `run(args) -> (bool, String)` helper and a
   `nonexistent_out()` no-side-effects idiom, and its module doc already covers `--clump_only`
   message precedence. Add: arm 5's refusal on `test_files/ubam_paired_test.bam` (the only
   non-trivial condition, and the only one whose *accept* twin has a binary-level test while its
   refuse case would not), the `--clump_only --paired --fastqc one.bam` precedence case, and —
   added in round 2 — arm 5's *acceptance* twin under `--output-format ubam`. Three, not two.

6. **`--clump_only --fastqc_args` regression test** — `--clump_only --fastqc_args "--quiet"` now
   produces `*_clumped_fastqc.{html,zip}`. Red today; it is what makes step 2 self-proving.

7. **Help text** — `--fastqc` (`cli.rs:327-331`) gains one clause naming the restriction.
   `--fastqc_args` (`:333-338`) keeps "Implies --fastqc", which step 2 makes true everywhere.

8. **Docs** — four lines, not one:
   - `guide/outputs.md:143` — "Works on both FASTQ output (default) and uBAM output" is a blanket
     promise with no exception.
   - `guide/outputs.md:145` — "PE FASTQ produces two reports (one per mate)" is false for the
     interleaved-uBAM→FASTQ arm, and separately false under `--retain_unpaired` (four output files,
     two reports).
   - `quickstart.md:53` — "Works on both FASTQ and uBAM output paths."
   - `guide/flags.md` "Cross-flag interactions" — the canonical home for refusals; has no FastQC
     entry.

9. **`CHANGELOG.md`** — one sentence with the `**Breaking:**` prefix, since previously-succeeding
   invocations now exit non-zero.

---

## 6. Efficiency

Five boolean tests plus one `Option::is_some()`, once, at argument-validation time before any I/O;
messages allocate only on the error path. `fastqc_requested()` is an inlined bool disjunction. No
memory or scalability dimension.

---

## 7. Integration

- **Reads/writes:** `src/cli.rs` (accessor, guards, tests, help), `src/main.rs` (five one-word
  substitutions), `tests/integration_preflight_tripwire.rs` (deletions + five new cases),
  `tests/integration_paired_format_guard.rs` (two cases), one `--clump_only` regression test, four
  docs lines, one CHANGELOG line.
- **`src/clump_only.rs`:** four `fastqc: bool` parameters renamed `fastqc_requested` (round 2).
  The sixth arm existed *because* a parameter named `fastqc` invites `cli.fastqc`; the rename is the
  only change that stops it recurring. No signature shape change, no call-site change.
- **Not touched:** `src/specialty.rs`, `src/fastqc.rs`, and every `fastqc::run` call site.
- **CI:** the four `--fastqc` invocations (`ci.yml:584`, `:1146`, `:1173`, `:1191`) are all on
  capable arms and must keep passing. Three have local `cargo test` equivalents (§9 V10). The Perl
  validation matrix passes no `--fastqc`, so byte-identity is unaffected.
- **#422's tripwire:** loses its FastQC allowance scoping — a workaround for this defect, as its
  own comment says — and gains two tests: one looping over all five refused arms, and (round 3)
  `clump_only_se_fastqc_args_only`, which puts the sixth arm under the file's own
  `active && allowed.is_empty()` rule. Verified: reverting the five pass-throughs makes that case
  fail with "a FastQC flag is set but no FastQC artefact was written", so the mechanism that missed
  the original defect now covers it.
- **Follow-on:** none filed. Two adjacent "which outputs get QC'd?" questions are recorded in §11
  rather than opened as issues, per the filing bar set on 2026-08-13.

---

## 8. Assumptions

| | Assumption | Status |
|---|---|---|
| A1 | Arms 1–5 are the complete set from which `fastqc::run` is unreachable | **Verified** by both reviewers, tracing all nine dispatch arms; four confirmed empirically. Multi-pair specialty input adds no arm (the guards key on the flag); `--demux`, `--passthrough`, `--clumpify` and `--basename` are all inside capable arms |
| A2 | `paired && len == 1 && !ubam && !clump_only` identifies arm 5 exactly | **Verified** independently by both, from dispatch order (`main.rs:845` above `:853`) rather than from the plan's prose |
| A3 | `--paired --output-format ubam <one.bam>` is capable and must stay accepted | **Verified** — `main.rs:2522`, `tests/integration_ubam_out.rs:572` (local), `ci.yml:1173` (CI) |
| A4 | ~~The disjunction is the gate at all 12 call sites.~~ **FALSE** | **Falsified by both reviewers.** It is the gate at 5 `main.rs` expressions (7 calls); the 5 `clump_only.rs` calls read a bare `fastqc: bool`. Corrected in §2.1, fixed by step 2 |
| A5 | ~~All four `--clump_only` drivers are capable, so no refusal is needed.~~ **Half true** | Capable for `--fastqc`, **not** for `--fastqc_args` alone. No refusal is still right — the fix is to wire the flag (§2.3) |
| A6 | Deleting the tripwire's `fastqc_capable` machinery breaks no case | **Verified**, and more strongly than v1 claimed: the branch is *unexercisable* because every `compare()` caller passes `true`. Deletion scope corrected in §5 step 4 |
| A7 | Refusing (rather than wiring) is right for arms 1–5 | **Decided by Felix**, 2026-08-12. Not re-litigated; no evidence against it surfaced, and §2.6's nf-core check argues the blast radius is small |
| A8 | Siting the guards late in `validate` preserves every existing message | **Verified** — the pre-empted checks sit at `cli.rs:996/1001/1004/1007/1020`, all below §3.4a. Pinned by two precedence tests |

---

## 9. Validation

| | Verify | How | Expected |
|---|---|---|---|
| V1 | The defect exists before the fix | `--hardtrim5 20 --fastqc a.fq` on `dev` | exit 0, no `*_fastqc.*`, "fastqc" absent from stderr |
| V1b | The **sixth** arm's defect exists before the fix | `--clump_only --fastqc_args "--quiet" a.fq` on `dev` | exit 0, no `*_fastqc.*` — reproduced by Reviewer A |
| V2 | Each of arms 1–5 is refused | Unit tests over `validate`, one per arm | non-zero, message names the arm |
| V3 | **Arm 5 did not over-reach** | Unit: `--paired --output-format ubam --fastqc one.bam` | **accepted** — catches a dropped `output_format` term |
| V4 | The refusal keys on the disjunction | Unit: `--fastqc_args "--quiet"` + `--hardtrim5`, no `--fastqc` | refused |
| V5 | **The sixth arm is fixed** | Integration: `--clump_only --fastqc_args "--quiet" a.fq` | `*_clumped_fastqc.{html,zip}` appear. **Red today** |
| V6 | Capable arms keep working | Unit + existing tests: SE trim, PE FASTQ, both other uBAM-out arms, `--clump_only --fastqc` | all accepted |
| V7 | Existing messages keep precedence | Unit: `--hardtrim5 0 --fastqc a.fq`; `--clock --fastqc r1 r2 r3` (no `--paired`, so it reaches `:1015`) | range error; `--clock mode requires` |
| V7b | A malformed `--adapter2` still reports first | Unit: `--paired --fastqc --adapter2 "NOT!A!SEQ" one.bam` | `only DNA characters`. Added by iteration #5; spell the flag long, `parse_from` skips `rewrite_perl_short_flags` |
| V7c | **Plain paired FASTQ + `--fastqc` stays accepted** | Unit: `--paired --fastqc R1 R2`, and `--fastqc R1` | accepted. Added in round 2 — dropping `input.len() == 1` refuses the commonest FastQC run and left all 87 unit tests green |
| V8 | Arm 5 refuses through the binary | `integration_paired_format_guard.rs`, `--paired --fastqc ubam_paired_test.bam` | non-zero, message, nothing written |
| V9 | `--clump_only --paired one.bam --fastqc` keeps its better message | Integration | `main.rs:597`'s wording |
| V10 | The tripwire's five new refusal cases fire | `cargo test --test integration_preflight_tripwire` | 53 + 5 cases pass; each asserts non-zero + nothing written |
| V11 | No regression | `cargo test`, `--release`, `cargo nextest run`, clippy `-D warnings`, `cargo fmt --check`, docs build | 663 + new tests, zero warnings |
| V12 | The capable `--fastqc` shapes are still exercised end-to-end | `cargo test --test integration_ubam_out --test integration_clump_only_ubam` | green. Covers `ci.yml:584/1146/1173`; only `:1191` is CI-only |

**V12 replaces v1's V10, which was not runnable.** Its four CI steps write to `/tmp` and use
`rm -rf` (both unavailable here), need a `--release` build, and `ci.yml:584` is preceded by
`! command -v fastqc`, which fails on any machine with FastQC installed. Three of the four already
have local equivalents, so `cargo test` is the honest gate; only `ci.yml:1191` is CI-only and CI on
the PR covers it.

---

## 10. Questions or ambiguities

The design decision (refuse, not wire or warn) was taken before drafting; the sixth arm, guard
placement and arm 5's wording were settled by dual review.

**Resolved at review:**

- **O1 — message precedence for `--paired --fastqc single.fastq`.** v1 accepted a wart. Both
  reviewers found the wording was not merely mis-ordered but **false** for that user, and its
  remedy did not work for them. Resolved by rewording to name the shape (§5.1), not by relocating
  the check to `main.rs` — which both reviewers also rejected, as it splits the family across two
  files for a message-ordering nicety.
- **O2 — should messages name a remedy?** Yes, one clause; arm 5 names three, because it has both
  a working alternative and a likely-different user.
- **O3 — `--help` scope.** One clause on `--fastqc`. `--fastqc_args` needs no edit because step 2
  makes its existing "Implies --fastqc" true.

**Open:** none.

---

## 11. Self-Review

**Logic.** Arms 1–5 are a CLI-level refusal; the sixth is a one-word wiring fix. No dispatch,
capability or output path moves, and no `fastqc::run` call site changes.

**Edge cases addressed.** The `--fastqc_args`-only path in both directions (refused on arms 1–5,
now working on `--clump_only`); the capable look-alike `--paired --output-format ubam <one.bam>`;
`--clump_only --paired` with one input, whose better message is preserved by an explicit exclusion;
four message-precedence cases that guard placement would otherwise have inverted; and the
FASTQ-single-file user, for whom arm 5's message must be true as well as ordered.

**Efficiency.** Nothing to consider.

**Integration.** The real couplings are #422's tripwire (a deletion plus five stronger cases) and
`main.rs`'s five pass-through sites, which v1 had listed as "not touched".

**What dual review changed.** v1's highest-risk claim — arm 5's four-term condition — was
**correct**, and both reviewers confirmed it independently from the dispatch. What was wrong was
the claim v1 marked "Verified" most confidently:

1. **A4 was false.** v1 asserted the disjunction was the gate at all 12 call sites, on the strength
   of grepping for `fastqc::run`. Grep finds *calls*, not *guards*; five of them read a different
   predicate, and that hid a sixth instance of the very defect being fixed. The generalisable check
   is "for each call site, what is the guard expression?" — twelve sites, two answers.
2. **Guard placement was wrong**, and would have silently pre-empted four existing messages,
   including a malformed-value error. The two reviewers disagreed on severity — one called the
   shifts benign, one called it a regression — and the second is right: reporting an unsupported
   combination ahead of `--hardtrim5 must be between 1 and 999` inverts the useful order, and the
   codebase argues that way in as many words at `main.rs:344-351`.
3. **Arm 5's message was false, not just early.** It told a user who passed one FASTQ that their
   input was an interleaved uBAM, and offered a remedy that fails for them.
4. **Remaining-risk #1 was wrong on the facts.** v1 said CI was the only end-to-end guard for three
   of the four capable shapes; `tests/integration_ubam_out.rs:538-572` is a "--fastqc regression
   guards" block whose own comment says "caught locally, not only in CI".
5. Smaller: the deletion inventory was short by seven sites; the docs scope by three lines; the
   integration-test home was the wrong file; §2.1's table listed 11 of 12 sites.

**Remaining risks.**
1. **Over-rejection is still the only way this can do harm.** V3, V6, V7 and V12 are the guards;
   both reviewers hunted for a broken working invocation and neither found one.
2. Users passing `--fastqc` unconditionally on a specialty arm get a hard failure where they got
   silence. That is the point; the CHANGELOG carries `**Breaking:**`.
3. The tripwire deletion touches a file whose purpose is anti-decay. Mitigated by replacing a
   never-fired guard with five that fire every run.

**Recorded rather than filed, per the 2026-08-13 filing bar.** Two adjacent "which outputs get
QC'd?" questions surfaced in review and are *not* #421's defect, because a report is produced in
both cases: `--fastqc` with `--demux` reports on the pre-demux trimmed file (`main.rs:1458` runs
before `:1468`), and PE `--retain_unpaired` produces no report for the `*_unpaired_*` files. Both
are undocumented rather than broken. If either ever warrants an issue, the user-visible symptom to
lead with is "the report describes a different file than you think".

---

## Implementation notes (2026-08-13)

Implemented on `dev` at `3ffcef4`. The release binary was built **before** any edit, so the
pre-fix matrix below is measured on genuinely unmodified code — worth stating because the build
and the first edit raced, and I checked the binary was still pre-fix (`--hardtrim5 --fastqc`
exited 0) rather than assuming the ordering held.

**Pre-fix matrix, measured. `artifacts` counts `*_fastqc.*` in the output dir.**

| Invocation | exit | artifacts | "fastqc" in log |
|---|---|---|---|
| `--hardtrim5 20 --fastqc` | 0 | 0 | 0 |
| `--hardtrim3 20 --fastqc` | 0 | 0 | 0 |
| `--clock --paired --fastqc` | 0 | 0 | 0 |
| `--implicon --paired --fastqc` | 0 | 0 | 0 |
| `--paired --fastqc <one.bam>` | 0 | 0 | 0 |
| **`--clump_only --fastqc_args "--quiet"`** | 0 | **0** | 0 |
| `--fastqc` (SE trim, control) | 0 | 2 | 4 |
| `--clump_only --fastqc` (control) | 0 | 2 | 1 |
| `--paired --output-format ubam --fastqc <one.bam>` (control) | 0 | 2 | 1 |

**Post-fix, same matrix:** all five arms exit 1 with a message naming the arm;
`--hardtrim5 --fastqc_args` alone also exits 1; **the sixth arm now produces 2 artifacts**; all
three controls still produce 2. `--hardtrim5 0 --fastqc` reports `must be between 1 and 999` and
an odd input count reports the even-count error — the two precedence cases the guard placement
was moved for.

**What was built.** `Cli::fastqc_requested()` (`cli.rs:585-593`) plus five guards sited after the `--adapter2`
parse, the last fallible check in `validate`; `main.rs` now calls the accessor at **ten** sites — the five
`--clump_only` pass-throughs *and* the five gate expressions, which previously spelled the
predicate out. 12 unit tests, 5 integration tests, one combined tripwire case covering five arms,
four docs edits, one help clause, one CHANGELOG line.

**Deviations from the plan.**

1. **The accessor replaced the five gate expressions too**, not just the five pass-throughs. §4
   said the accessor "replaces a predicate currently spelled out at five `main.rs` gate
   expressions", so this is the plan's intent; the outline only spelled out the pass-throughs.
   `grep "fastqc || .*fastqc_args"` now matches exactly once, in `cli.rs`.
2. **The five tripwire refusal cases are one `#[test]` looping over five arms**, not five
   separate tests. The file's convention is one `#[test]` per *dispatch arm* for the tripwire
   cases proper, but these are refusals, and the loop keeps the shared fixture setup in one
   place. Case count 53 → 54 rather than → 58.
3. **`--fastqc_args`' help text was left alone** (§5 step 7 predicted this): its "Implies
   --fastqc" is now true everywhere, so there is nothing to qualify.

**Iteration log.**

- **#1** — the background release build raced my first `cli.rs` edit. Verified the binary was
  still pre-fix by running `--hardtrim5 --fastqc` and checking for exit 0 before trusting any
  pre-fix measurement.
- **#2** — the `--implicon` probe exited 1 for an unrelated reason: `--implicon` is
  `require_equals = true`, so `--implicon 8` parsed `8` as a third input file and tripped the
  even-count check. Re-probed with bare `--implicon`; the silent no-op was then confirmed. A
  probe that fails for the wrong reason is worse than no probe.
- **#3** — the first tripwire case used `std::fs::copy("test_files/…")`; the file documents a
  `copy_fixture` helper resolving from `CARGO_MANIFEST_DIR` precisely because cases run with the
  case root as CWD. Switched to the helper.
- **#4** — `tempdir()` in `integration_paired_format_guard.rs` returns `PathBuf`, not `TempDir`;
  two `dir.path()` calls failed to compile. Fixed.

- **#5** — my own enumeration of `validate` (prompted by the code reviewers dying before they
  could do it) found the guard block still masked one message. `Cli::validate` has ~60 fallible
  checks before the block and exactly one after it: `parse_adapter_specs_quiet(&self.adapter2)`.
  So `--paired --fastqc --adapter2 "NOT!A!SEQ" one.bam` reported the FastQC refusal instead of
  `Adapter sequence must contain only DNA characters`. That is the same
  malformed-value-loses-to-unsupported-combination shape the block was moved late to avoid, one
  check short. Moved below the `-a2` parse; §2.4 and §5 step 1 corrected. Pinned by
  `malformed_adapter2_precedes_the_fastqc_refusal`.
- **#6** — that test first failed with `Input file not found: NOT!A!SEQ`. **`Cli::parse_from`
  bypasses `rewrite_perl_short_flags()`**, so `-a2` was never rewritten and its value parsed as a
  positional. Spelled `--adapter2` instead. Worth knowing generally: a `parse_from` test that
  asserts a *refusal* using a Perl-era short flag can pass for the wrong reason, since the run
  errors either way.

### Round 3 — dual code review + coverage audit (2026-08-13)

`CODE_REVIEW_A_gen2.md`, `CODE_REVIEW_B_gen2.md`, `COVERAGE_gen2.md`. Three of the first four agents
died on infrastructure errors having produced nothing; the gen-2 set completed. **Coverage verdict:
COMPLETE** — 48 items, 45 DONE, 0 MISSING, 3 documented deviations, every claimed number
re-measured. **Both code reviewers: no over-rejection**, no Critical, no High.

**Both reviewers independently found the same defect I introduced:** the new tripwire `#[test]` was
inserted *between* `a_refused_run_writes_nothing` and its `#409` doc comment, re-attaching that
comment to the FastQC test and orphaning the contract record — in the file whose whole purpose is
anti-decay.

**The reviewers' deaths produced one finding of their own.** With three agents dead I did the
read-only enumeration myself, which is how iteration #5 (the `-a2` masking) was found. Reviewer A
then re-did that enumeration independently and confirmed the *new* position is right: `validate` ends
at `cli.rs:1131`, and the only statements below the block are five `eprintln!` deprecation warnings,
so the block is genuinely the last fallible check. That was the one question I could not judge,
having both made and diagnosed the change.

**Fixed in round 3:**

- `--help` still asserted "given a single interleaved uBAM" — the exact claim round 2 removed from
  the runtime message (A's M-3). **The corrective-copy pattern, third occurrence in this issue
  alone**: fix the message, leave the identical claim in the help text.
- `integration_clump_only.rs`'s doc comment: historical narrative *and* stale after the parameter
  rename (it still said `fastqc: bool`). Cut to one line.
- `integration_paired_format_guard.rs`'s module doc opened "Three properties are pinned here" over a
  closed numbered list; the diff added a fourth. Extended.
- Two comments trimmed (L-1, L-2); the surviving line is the invariant that stops the block being
  moved back up — the mistake iteration #5 made.
- Unit coverage for arm 5's `!clump_only` term (L-5): mutation M7 was caught only by an integration
  test, so `cargo test --lib` left it unguarded.
- The tripwire case for the sixth arm (L-6), mutation-verified.

**One of A's suggestions was wrong, and finding out was worth more than the suggestion.** L-4 asked
for `!out.exists()` on `clump_only_paired_single_file_keeps_its_own_message`, reasoning that a refused
run writes nothing. It **fails**: `Cli::validate` runs at `main.rs:231`, before `ensure_output_dir` at
`:420`, but `main.rs:599`'s `--clump_only` refusal runs *after* it. Validate-layer refusals create
nothing; main-layer refusals create the output directory. That asymmetry is pre-existing and is now
pinned in the direction it actually holds, with the contrast to its sibling documented.

**Still open by decision:** nothing. The two adjacent "which outputs get QC'd?" questions stay
recorded in §11 rather than filed, per the 2026-08-13 filing bar.

---

## Revision History

- **v1** (2026-08-12) — first draft. Arm conditions read off the dispatch and CI workflow rather
  than the issue text, which identified arm 5's `output_format` term and the `--clump_only`
  message-precedence exclusion.
- **v2** (2026-08-13) — dual plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`). Arm 5's
  condition **confirmed** by both, independently. **A4 falsified:** the gate is not uniform, and a
  sixth silent no-op (`--clump_only --fastqc_args`) is live — folded into this PR as a wiring fix
  rather than deferred to a new issue. Guard block moved late in `validate` to stop it pre-empting
  four existing messages. Arm 5's message reworded to name the shape rather than assert a format.
  Tripwire deletion inventory completed and paired with five replacement refusal cases. Docs scope
  widened to four lines; integration-test home corrected; v1's unrunnable V10 replaced by V12;
  Remaining-risk #1 corrected.
