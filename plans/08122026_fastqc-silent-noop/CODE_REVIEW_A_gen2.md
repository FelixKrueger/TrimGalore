# CODE REVIEW A (gen2) — #421 `--fastqc` silent no-op

**Reviewer:** A (fresh context, no shared state with Reviewer B)
**Date:** 2026-08-13
**Base:** `dev` @ `3ffcef4`, working diff uncommitted
**Method:** all builds, tests, mutations and binary probes run in an isolated worktree
(`$TMPDIR/cr-A2-421`), never in the main tree.

---

## Important: the tree changed mid-review

I began against a snapshot with 9 modified files (`src/cli.rs` +213). Partway through, a **round-2
edit landed** in the main tree: 10 files, `src/cli.rs` +236, and `src/clump_only.rs` newly modified.
I detected this while cross-checking a `clump_only.rs` parameter name against `3ffcef4`, re-synced my
worktree, and **re-ran every check against the current state.** Everything below is against the
current `git diff`, not the earlier snapshot.

Round 2 independently fixed eight things I had recorded against snapshot 1 — I list them under
*Already fixed in round 2* so they are not re-raised, and because the fact that a second reviewer
found them matters more than which of us wrote them down first. Three of my findings survive round 2,
plus eight lower-priority ones.

---

## Summary

**The change is correct.** Over-rejection — the only way this can do harm — does not occur: I
characterised the newly-refused set exactly from the guard predicate, verified the two "no
`fastqc::run` is reachable" premises by direct inspection rather than by grep, and swept 16
invocations through a purpose-built binary. Guard placement is right, and I confirmed it
independently rather than taking the plan's word: `Cli::validate` ends at `cli.rs:1131`, and the only
statements below the #421 block are five `eprintln!` deprecation warnings — the block *is* the last
fallible check, and nothing below it should precede it.

All nine mutations I designed go red, including the four the team lead predicted. The tripwire
deletion is safe and the surviving rule is strictly wider than what it replaced.

No Critical or High findings. The three Medium findings are all **stale or inaccurate prose in
user-facing or contract-bearing text** — the same failure mode as the last five PRs in this repo
(`feedback_corrective_copy_defect`): round 2 corrected the runtime message that named the input
format wrongly, and left the identical claim standing in `--help`.

**Verification reproduced independently in my worktree:** 682 tests pass (`cargo test --release
--no-fail-fast`, 15 binaries), `cargo clippy --all-targets --release -- -D warnings` clean,
`cargo fmt --all -- --check` clean.

---

## 1. Over-rejection (the load-bearing question)

### 1.1 The newly-refused set, characterised exactly

A run is newly refused iff `fastqc_requested() && (hardtrim5 | hardtrim3 | clock | implicon | arm5)`,
where arm5 = `paired && input.len() == 1 && output_format != UBam && !clump_only`. So over-rejection
is possible only if one of those arms *can* reach `fastqc::run`. Both premises hold:

- **`src/specialty.rs` contains zero occurrences of `fastqc`** (`grep -c` → 0), covering arms 1–4
  including their `--output-format ubam` variants (`specialty::hardtrim5_to_bam` etc.).
- **`run_paired_ubam_single_file` (`main.rs:1803–2003`) contains no `fastqc::run`.** All 12 sites are
  at `main.rs:1459/1778/1784/1796/2242/2399/2522` and `clump_only.rs:374/554/555/847/1095`. Plan
  §2.1's table is accurate, all 12 accounted for.
- **Dispatch order matches the plan's table**, read off `main.rs` directly: `:464` hardtrim5 → `:495`
  hardtrim3 → `:526` clock → `:548` implicon → `:580` clump_only → `:845` uBAM output → `:853`
  paired+1 → `:867` paired. So `:845` sitting above `:853` is what makes the `output_format` term
  load-bearing, and the four specialty terms in arm 5 are genuinely redundant because those guards
  bail earlier in `validate`.
- **`fastqc_args` has no clap `default_value`** (`cli.rs:340`, plain `Option<String>`), so
  `fastqc_requested()` cannot be true on a run where the user passed neither flag.

### 1.2 Empirical sweep (release binary built from the current diff)

| Invocation | rc | reports | verdict |
|---|---|---|---|
| `--fastqc` SE | 0 | 2 | accepted |
| `--paired --fastqc r1 r2` | 0 | 4 | accepted |
| `--output-format ubam --fastqc one.bam` | 0 | 2 | accepted |
| `--paired --output-format ubam --fastqc one.bam` | 0 | 2 | accepted |
| `--clump_only --fastqc` | 0 | 2 | accepted |
| `--clump_only --fastqc_args --quiet` | 0 | **2** | **sixth arm fixed** |
| `--clump_only --paired --fastqc_args --quiet r1 r2` | 0 | 4 | driver 2 wired |
| `--clump_only --output-format ubam --fastqc_args --quiet` | 0 | 2 | driver 3 wired |
| `--clump_only --paired --output-format ubam --fastqc_args --quiet r1 r2` | 0 | 2 | driver 4 wired |
| `--demux … --fastqc` | 0 | 2 | accepted |
| `--paired --passthrough … --fastqc` | 0 | 3 | accepted |
| `--hardtrim5 20 --fastqc` | 1 | 0 | intended refusal |
| `--hardtrim5 20 --output-format ubam --fastqc` | 1 | 0 | intended refusal |
| `--hardtrim5 20 --fastqc` ×2 pairs (multi-pair) | 1 | 0 | intended refusal |
| `--clock --paired --output-format ubam --fastqc` | 1 | 0 | §3.4a's clock+ubam message wins ✓ |
| `--paired --fastqc single.fastq` | 1 | 0 | see §1.3 |

All four `--clump_only` drivers now honour `--fastqc_args` alone — that is the sixth arm proven end
to end, not just in the one regression test.

**Harness caveat worth recording:** my first sweep reported `clump_only --fastqc_args` → **0
reports**, which would have been a serious finding. It was a harness artefact: the mutation driver
restores sources but does not rebuild, so `target/release/trim_galore` was still the M8 (mutated)
binary. Rebuilt and re-ran; 2 reports. A negative from a binary you did not just build is not a
negative.

### 1.3 One residual message-ordering wart, and it is acceptable

`--paired --fastqc single.fastq` (user forgot Read 2) previously got `main.rs:334`'s message:

> `--paired with a single input file is only legal if that file is a uBAM. Provide R1 and R2 as separate files for FASTQ paired-end mode.`

It now gets the FastQC refusal instead. Round 2 reordered that message so the two-file remedy leads
and the uBAM remedy is qualified ("for a genuine interleaved uBAM"), which is the right fix — the
correct remedy is now first, and `fastqc_plus_paired_single_file_rejected` pins the ordering with a
real `find()` comparison rather than the vacuous negative assertion snapshot 1 had. `main.rs:334` is
still more direct for that user, but the wart is now cosmetic. Not worth relocating the check.

**Conclusion: no invocation that works on `3ffcef4` newly fails, other than the five arms this PR
exists to refuse.**

---

## 2. Guard placement — enumerated independently

`validate` runs `cli.rs:596–1131`. I enumerated its checks rather than trusting §2.4.

- **Last fallible check before the block:** `crate::adapter::parse_adapter_specs_quiet(&self.adapter2)`
  at `cli.rs:1054`. The block is at `1061–1101`.
- **Everything below the block:** five `eprintln!` deprecation warnings (`--gzip`,
  `--path_to_cutadapt`, `--cutadapt_args`, `--suppress_warn`, `--keep`, `cli.rs:1103–1128`), then
  `Ok(())`. **No fallible check follows the block**, so nothing is masked.
- **Nothing below it should precede it.** The only behavioural consequence of the placement is that a
  refused run no longer prints those deprecation warnings — irrelevant on a path that exits non-zero.

The upstream precedence chain is intact and I verified two links empirically: `--clock --paired
--output-format ubam --fastqc` reports §3.4a's clock+ubam message (top of `validate`), and
`--hardtrim5 20 --fastqc missing.fq` reports the input-not-found error from the
`check_restartable_input` loop at `:1030`. The order "malformed value / structurally invalid input →
unsupported combination" is preserved throughout, matching the argument the codebase already makes at
`main.rs:344–351`.

The `-a2` interaction is worth one note: `--hardtrim5 20 --fastqc --adapter2 GCAT a.fq` prints the
`-a2`-unused **warning** (`cli.rs:1050`) and then the FastQC **error**. A warning about an ignored
flag on a run that is about to fail is mild noise, not a defect.

---

## 3. The ten `fastqc_requested()` call sites

Confirmed 5 + 5, and that the split is exactly what the plan describes:

- **Five pass-throughs** (`main.rs:636, 671, 715, 778, 819`) feed the four `clump_only.rs` drivers.
  These are the semantic change: they previously passed bare `cli.fastqc`.
- **Five gate expressions** (`main.rs:1458, 1777, 2241, 2398, 2521`) previously spelled
  `cli.fastqc || cli.fastqc_args.is_some()`. Textually identical substitution — no behaviour change,
  and `grep -rn 'cli\.fastqc'` outside `main.rs`/`cli.rs` returns nothing, so there is no other
  reader.

**No site wants the narrower `cli.fastqc`.** `--fastqc_args`' help text (`cli.rs:336`) has promised
"Implies --fastqc" all along, and all four `clump_only` drivers are capable — so the narrow reading
was never the intent anywhere, which is exactly why the defect was silent.

**A detail worth putting on the record.** At `3ffcef4` the parameter was `fastqc: bool`; round 2
renamed it `fastqc_requested` at all four drivers. That rename is the most durable part of the fix:
the sixth arm existed *because* a parameter named `fastqc` invites `cli.fastqc` at the call site, and
a parameter whose name states the disjunction makes the wrong feed read as wrong. Plan §7 documents
this and correctly removes `clump_only.rs` from "Not touched".

---

## 4. Mutation testing — 9/9 caught

Applied in the worktree, `cargo test --release --no-fail-fast` over `--lib` plus
`integration_{preflight_tripwire, paired_format_guard, clump_only, ubam_out}`.

| # | Mutation | Tests that go RED |
|---|---|---|
| M1 | neuter `hardtrim5` guard | `fastqc_plus_hardtrim5_rejected`, `fastqc_args_alone_also_rejected`, `refuses_fastqc_on_every_arm_that_cannot_run_it` |
| M2 | neuter `hardtrim3` guard | `fastqc_plus_hardtrim3_rejected`, `refuses_fastqc_…` |
| M3 | neuter `clock` guard | `fastqc_plus_clock_rejected`, `refuses_fastqc_…` |
| M4 | neuter `implicon` guard | `fastqc_plus_implicon_rejected`, `refuses_fastqc_…` |
| M5 | neuter arm-5 guard | `fastqc_plus_paired_single_file_rejected`, `fastqc_refused_on_paired_interleaved_fastq_output`, `refuses_fastqc_…` |
| M6 | drop arm 5's `output_format != UBam` term | `fastqc_plus_paired_ubam_output_single_file_accepted`, `fastqc_accepted_on_paired_interleaved_ubam_output`, `ubam_out_pe_fastqc_produces_exactly_one_report` |
| M7 | drop arm 5's `!clump_only` term | `clump_only_paired_single_file_keeps_its_own_message` |
| M8 | revert the 5 pass-throughs to `cli.fastqc` | `clump_only_fastqc_args_alone_produces_a_report` |
| M9 | move the guard block above the `-a2` parse | `malformed_adapter2_precedes_the_fastqc_refusal` |

All four predictions in the brief confirmed. Two observations:

- **My first mutation run under-reported.** With cargo's default fail-fast, `--lib` failing meant the
  integration binaries never ran, so M1–M5 looked like they were caught by unit tests only. Re-run
  with `--no-fail-fast`, the tripwire fires on all five. Worth noting because "which layers catch
  this" was the actual question.
- **M7 and M8 are each caught by exactly one test** — the two thinnest points in the suite. See L5
  and L6.

---

## 5. The tripwire deletion

- **The surviving rule is load-bearing and strictly wider.** `compare()` collapses to `if active &&
  allowed.is_empty()` and now applies to all 54 cases instead of being scoped away on eight. The
  removed `(true, false)` arm was unexercisable (every `compare()` caller passed `true`), so nothing
  is lost. The deletion inventory is complete: field, builder, 8 call sites, the 6th parameter, the
  `Tripwire::run` argument and nine literal `true`s — clippy `-D warnings` clean confirms no unused
  binding survived.
- **`refuses_fastqc_on_every_arm_that_cannot_run_it` (`:1340`) asserts all three things**: `!ok`,
  `stderr.contains("--fastqc is not supported")`, and an empty `snapshot(&root).difference(&before)`.
  `case_root` (`:87–92`) does `remove_dir_all` then `create_dir_all`, so each of the five iterations
  gets a clean root and the diff is meaningful. The "nothing written" property is structurally sound
  (the guard is in `validate()`, before any I/O) rather than mutation-proven — I could not construct a
  cheap mutation that writes-then-refuses.
- **Round 2 fixed the placement bug I found in snapshot 1**: the new test had been inserted *between*
  `a_refused_run_writes_nothing` and its `#409` doc comment, silently re-attaching that comment to the
  FastQC test. It is now appended after, and `:1312–1314` documents `a_refused_run_writes_nothing`
  again.

---

## 6. `parse_from` and Perl-era short flags — no existing defect

The trap is real (`Cli::parse_from` bypasses `rewrite_perl_short_flags`, so `-a2` parses as a
positional and a refusal test passes for the wrong reason). I swept for it: `cli.rs` contains exactly
six raw `-a2`/`-r1`/`-r2` occurrences in tests (`:1446, :1454, :1518, :1520, :1533, :1564/:1566`) and
**all six pass through the `rewrite()` helper at `:1439` first**, including both
`test_rewrite_*_end_to_end_via_parse_from` tests. No test in `src/` or `tests/` has the defect.

Related, and a genuine catch by the implementer in round 2:
`odd_input_count_precedes_the_fastqc_refusal` **was** passing for the wrong reason. With `--paired`
set, `validate_paired_input("Paired-end")` fires at `~:650` and its message also contains "even
number of input files", so the assertion could not distinguish it from `--clock`'s check at `:1015` —
the check the test exists to pin. Round 2 dropped `--paired` and asserts `"--clock mode requires"`.
I verified the binary now emits `--clock mode requires an even number of input files (R1/R2 pairs),
got 3`. This is the same class of defect as the `parse_from` trap and the reason to look for it.

---

## 7. Docs and CHANGELOG accuracy

Both claims the brief asked about, checked against code:

- **`--retain_unpaired` produces no report for `*_unpaired_*`** — **true.** `run_paired`
  (`main.rs:1777–1798`) calls `fastqc::run` on `output_r1`, `output_r2` and the `--passthrough`
  carrier only; the unpaired writers are never passed.
- **"Four paths write no report"** — **was wrong** (it enumerated five: four specialty modes plus the
  paired-interleaved shape). Round 2 rewrote it to "Two situations … the four specialty modes …, and
  `--paired` given …", which is right.
- **`--passthrough`** — I verified empirically that `--paired --passthrough … --fastqc` writes
  **three** reports (`…_I1_passthrough_fastqc.zip`, `…_R1_val_1_fastqc.zip`, `…_R2_val_2_fastqc.zip`).
  Snapshot 1's rewrite of that sentence enumerated the `--retain_unpaired` exception but not this one;
  round 2 added it. Both are now accurate.
- **CHANGELOG** — one semicolon-joined sentence with the `**Breaking:**` prefix, matching the shape of
  the neighbouring #415 entry. Within house style.

---

## 8. Findings

### Medium

**M-1 — `tests/integration_clump_only.rs:402–404`: historical narrative, and now factually stale.**

```rust
/// #421 — the four `--clump_only` drivers gate FastQC on a `fastqc: bool` parameter,
/// which `main.rs` feeds from `Cli::fastqc_requested()`. Passing bare `cli.fastqc` made
/// `--fastqc_args` alone a silent no-op on this path.
```

Two defects in three lines. The last sentence describes the state being fixed, which
`~/.claude/CLAUDE.md` bans outright ("A comment is read by someone who has never seen the old code").
And `fastqc: bool` is **stale as of round 2** — the parameter is now `fastqc_requested: bool`, so the
comment misdescribes the code it points at. Suggested replacement, one line:

```rust
/// #421 — `--fastqc_args` alone must activate FastQC on the `--clump_only` drivers.
```

**M-2 — `tests/integration_paired_format_guard.rs:10`: the module doc is an enumerated contract, and
the diff adds a fourth item to a list of three.**

Line 10 reads "Three properties are pinned here that unit tests cannot reach:" followed by a numbered
1./2./3. list. The diff adds three tests (`:445`, `:474`, `:497`) pinning a fourth property — the
#421 refusal and its acceptance twin on the `--paired`-single-input family. Because the list is
numbered and closed ("Three properties"), leaving it makes the file's own contract wrong. Add a
fourth item naming the FastQC refusal.

**M-3 — `src/cli.rs:330–333`: the `--fastqc` help text still asserts the input format that round 2
removed from the runtime message.**

```
Rejected with the specialty modes (--hardtrim5/3, --clock, --implicon) and with
--paired given a single interleaved uBAM and FASTQ output, none of which write a
report. For the last, --output-format ubam does produce one; otherwise run FastQC
separately.
```

Round 2 correctly fixed the remedy scoping here (snapshot 1 offered `--output-format ubam` for all
five arms, where it helps only one — and I confirmed `--hardtrim5 20 --output-format ubam --fastqc`
is still refused, so the old wording sent users into a second failure). But **"given a single
interleaved uBAM" is the exact claim §5.1 established is wrong**: the guard fires on *any* single
input under `--paired`, and the likelier user forgot Read 2 and has no uBAM. The runtime message was
reworded to name the shape; the help text carrying the same false claim was not. Suggest "and with
`--paired` given a single input file (the interleaved-uBAM path) writing FASTQ output".

### Low

- **L-1 — `src/cli.rs:1087–1088`.** "Names the input shape rather than the format: this also catches a
  user who meant to pass Read 2, whose file is not a uBAM at all." Two lines justifying a
  message-wording choice — that argument belongs in the commit message. It also sits slightly at odds
  with the message it guards, which *does* name the format ("the interleaved-uBAM paired path").
- **L-2 — `src/cli.rs:1059–1060`.** "…so the flag would write no report and say nothing" is a
  counterfactual justification for the guard's existence. The second sentence — "Last of the fallible
  checks: malformed values report first" — is the part that earns its place: it is the invariant that
  stops the next person moving the block up, which is precisely the mistake iteration #5 made. Keep
  that line, drop the first.
- **L-3 — `tests/integration_paired_format_guard.rs:474` is redundant.**
  `fastqc_accepted_on_paired_interleaved_ubam_output` and
  `integration_ubam_out.rs::ubam_out_pe_fastqc_produces_exactly_one_report` went red together on M6,
  so it adds no discriminating power, and it costs a full FastQC run. Keep only if the local contrast
  with the refusal case two tests above is worth the runtime.
- **L-4 — `tests/integration_paired_format_guard.rs:505`.**
  `clump_only_paired_single_file_keeps_its_own_message` passes `nonexistent_out(&dir)` to `-o` but
  never asserts `!exists()` afterwards. The helper's documented purpose (`:44–45`) is to make that
  absence assertion non-vacuous; using it without the assertion silently drops the property. The
  sibling at `:445` does assert it.
- **L-5 — arm 5's `!clump_only` term has no unit coverage.** M7 is caught only by an integration test,
  so `cargo test --lib` leaves it unguarded. One `validate()` unit test asserting `--clump_only
  --paired --fastqc one.bam` does *not* produce the FastQC message would close it cheaply.
- **L-6 — the tripwire has no `--clump_only --fastqc_args` case.** The sixth-arm fix is guarded only by
  `integration_clump_only.rs:406`. A tripwire case would put it under the file's own `active &&
  allowed.is_empty()` rule — which is the mechanism that would have caught the original defect, and
  the reason it did not is that no case paired `--clump_only` with `--fastqc_args` alone. Adding one
  closes the loop the plan's §11 point 1 identifies.
- **L-7 — tripwire `:36` says "The `refuses_fastqc_*` cases below"**, plural, for a glob matching one
  `#[test]`. Defensible (five loop iterations) but reads as several functions. Deviation #2 chose one
  test deliberately; the doc predates that choice.
- **L-8 — minor plan/code drift** (not blocking, but the plan is the audit trail): §5 step 5 says two
  cases in `integration_paired_format_guard.rs` where three landed; §7 says the tripwire "gains five
  active refusal cases" where it gained one test with five iterations (documented in deviation #2, but
  §7 was not updated); §5 step 2's "No `clump_only.rs` signature changes" is now only true in the
  narrow sense that the *shape* is unchanged, which §7 already explains better.

### Already fixed in round 2 (recorded, not re-raised)

Found against snapshot 1; all verified fixed in the current tree. Listed so the coverage of the
review is auditable, and because several are the "corrective copy contains the defect" pattern.

1. Tripwire: the new test was inserted between `a_refused_run_writes_nothing` and its `#409` doc
   comment, re-attaching that comment to the wrong function.
2. Tripwire `:219`: "Production's gate is `cli.fastqc || cli.fastqc_args.is_some()` (`main.rs:1458` and
   four siblings)" — stale after the substitution, in the one file whose purpose is anti-decay.
3. Tripwire test body: "Replaces the `fastqc_incapable` marker … That marker was latent …" —
   historical narrative about deleted code.
4. `outputs.md`: "Four paths write no report" enumerating five.
5. `outputs.md`: the report-count sentence omitted `--passthrough`'s third report.
6. `--fastqc` help text offered `--output-format ubam` as a remedy for all five arms.
7. `fastqc_plus_paired_single_file_rejected` asserted `!err.contains("your input is")` — a negative
   against a string no candidate message ever contained, so it could not fail.
8. `odd_input_count_precedes_the_fastqc_refusal` passed via the generic `"Paired-end"` even-count
   message rather than `--clock`'s (§6).

---

## 9. Recommendations, prioritised

**Critical / High:** none.

**Medium — worth fixing before merge:**
1. M-3 — `--fastqc` help text: stop asserting the input is a uBAM (`cli.rs:330–333`). This is
   user-facing and is the *same* claim round 2 just removed from the runtime message.
2. M-1 — rewrite `integration_clump_only.rs:402–404` to one line; it is both historical narrative and
   stale on the parameter type.
3. M-2 — extend `integration_paired_format_guard.rs:10`'s numbered list to four.

**Low — batch or defer:**
4. L-5 and L-6 — the two single-test coverage points (arm 5's `!clump_only` term; the sixth arm under
   the tripwire's own rule). Both are one small test each and would make the suite proof against the
   thinnest mutations.
5. L-1, L-2 — trim two comments to one line each.
6. L-4 — add the missing `!out.exists()` assertion.
7. L-3, L-7, L-8 — cosmetic.

**Merge verdict:** approve after the three Medium prose fixes. The behaviour, the guard placement, the
ten call sites and the tripwire rewrite are all correct and adequately pinned; nothing I found touches
runtime behaviour.

---

## Appendix — what I ran

- Worktree: `git worktree add --detach "$TMPDIR/cr-A2-421" 3ffcef4`, then `git apply` of the working
  diff. Re-synced (`git checkout -- . && git apply`) after detecting the round-2 edit. Main tree never
  built in.
- `cargo build --release`; `cargo test --release --no-fail-fast` → **682 passed, 0 failed** across 15
  binaries; `cargo clippy --all-targets --release -- -D warnings` → clean; `cargo fmt --all -- --check`
  → clean.
- Mutation driver (9 mutations × full re-test):
  `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/d03b1224-8e75-4004-b69d-2e4f77c29bdd/scratchpad/mutate.py`
- 16-invocation binary sweep for over-rejection (§1.2), plus targeted probes for message precedence
  (`--clock` odd count, `--clock`+ubam, arm 5 remedy ordering) and the `--passthrough` report count.
