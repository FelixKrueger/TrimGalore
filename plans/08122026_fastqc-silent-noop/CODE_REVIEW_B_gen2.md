# CODE REVIEW B (gen2) — #421 `--fastqc` silent no-op

**Reviewer:** B (independent, fresh context; gen2 — the gen1 instance died before writing)
**Base:** `dev` @ `3ffcef4`, work uncommitted
**Method:** read-only in the main tree; all execution in a throwaway worktree
`$TMPDIR/cr-B2-421` (detached @ `3ffcef4` + `git apply` of the submitted diff). The main tree was
never built, tested or edited — `git diff --stat` is unchanged, and the worktree was byte-compared
back to the submitted patch after every mutation.
**Verdict:** the behaviour change is correct. I found **no over-rejection** — no invocation that
works on `3ffcef4` and now fails. Ten defects, none in the guard logic: one mis-attached doc
comment, one wrong remedy in `--help`, one diagnostic regression for the likeliest user error, and
one test that does not test what it claims.

> **The working tree changed while this review was running.** Findings below are against the diff I
> was handed (snapshotted at start). Between then and writing, `src/cli.rs` was edited to move the
> guard block below the `#369 -a2` block, shorten its comment to two lines, and add
> `malformed_adapter2_precedes_the_fastqc_refusal`. That independently resolves **LOW-8** and
> **LOW-9 item 1** — both are marked RESOLVED in place. I re-verified the new state in my worktree:
> `cli::tests` 88/88 (was 87), tripwire 54/54, PFG 14/14, `fmt` clean, and mutation **M8** (block
> back above `-a2`) turns the new test red, so it is a real placement pin. **Every other finding
> below still stands against the current tree** — the tripwire doc comment, the `--help` remedy, the
> arm-5 message order, the missing acceptance twins, the `--paired` precedence test and the docs
> count are all untouched.

---

## 1. Summary

| Question asked | Answer |
|---|---|
| Over-rejection? | **None found.** Arms 1–4 cannot over-reject (dispatch is unconditional on the flag); arm 5's four terms exactly reproduce the `main.rs:853` dispatch |
| Guard placement | **Correct, and protects far more than the plan claims.** One malformed-value check (`-a2`) is still pre-empted, contrary to the block's own comment |
| Ten `fastqc_requested()` sites | **All ten semantically identical.** No bare `cli.fastqc` read remains in `src/` |
| Can the new tests fail? | **Yes — 7 mutations, all caught.** Except one: the mutation that breaks `--paired --fastqc R1 R2` leaves all 87 `cli::tests` green |
| Tripwire deletion | **Safe, and the surviving rule is live** — proved by mutation, fires across 8 cases |
| Docs / CHANGELOG accuracy | One wrong count, one stale claim; the `--retain_unpaired` claim is **true** (verified in code) |
| Comment style | 3 comments over the two-line max; 2 narrate the state being fixed, one of them naming a now-deleted identifier |

Independent green run in the worktree: `cli::tests` 87/87, `integration_preflight_tripwire` **54/54**,
`integration_clump_only` 13/13, `integration_paired_format_guard` 14/14, `cargo fmt --check` clean.

---

## 2. Issues

### HIGH-1 (Structure) — the new tripwire test stole `a_refused_run_writes_nothing`'s doc comment

`tests/integration_preflight_tripwire.rs:1312-1316`. The new `#[test]` was inserted **between** an
existing doc comment and the item it described. Rust attaches `///` to the *following* item, so:

```rust
/// The tripwire tests the pre-flight's *coverage*; this pins its *purpose*. Two inputs
/// collapsing onto one output under `-o` must be refused before anything is written —
/// the contract #409 violated.
#[test]
fn refuses_fastqc_on_every_arm_that_cannot_run_it() {   // ← inherits a comment about `-o` collisions
...
#[test]
fn a_refused_run_writes_nothing() {                      // ← the `-o` collision test, now undocumented
```

The new FastQC test now carries a description of output-path collisions, and the #409 purpose record
is orphaned from the test that holds it. In a file whose entire reason for existing is anti-decay,
and against a house rule that comments must describe the current state of the thing they sit on, this
is the kind of defect that survives for years.

**Fix:** move `refuses_fastqc_on_every_arm_that_cannot_run_it` *below*
`a_refused_run_writes_nothing`, or give it its own `///` and restore the original comment directly
above the collision test. Zero risk.

### HIGH-2 (Errors) — `--help` offers a remedy that does not work on four of the five refused arms

`src/cli.rs:330-332`, rendered by `--help` as one paragraph:

> Rejected with the specialty modes (--hardtrim5/3, --clock, --implicon) and with --paired given a
> single interleaved uBAM and FASTQ output, none of which write a report; **add --output-format ubam**
> or run FastQC separately

The remedy reads as applying to the whole list. It does not. Verified in the worktree:

```
$ trim_galore --hardtrim5 20 --output-format ubam --fastqc BS-seq_10K_R1.fastq.gz
Error: --fastqc is not supported with --hardtrim5, which writes no QC report; …
```

The guard keys on `hardtrim5.is_some()` alone, so `--output-format ubam` changes nothing for four of
the five arms. A user reading `--help` will try it and get the same refusal.

**Fix:** attach the remedy to its case — e.g. `…none of which write a report; for the last,
--output-format ubam does produce one, otherwise run FastQC separately`.

### MEDIUM-3 (Logic) — the forgot-R2 user now gets a worse message than before

Not over-rejection (the run fails either way), but a diagnostic regression on what the plan itself
calls the likelier mistake. `cli.validate()` runs at `main.rs:231`, *above* `main.rs:334`. Verified
in the worktree:

| Invocation | `3ffcef4` | with the fix |
|---|---|---|
| `--paired r1.fastq.gz` | `--paired with a single input file is only legal if that file is a uBAM. Provide R1 and R2 as separate files for FASTQ paired-end mode.` | unchanged |
| `--paired --fastqc r1.fastq.gz` | **same, correct message** | the arm-5 refusal |
| `--paired --output-format ubam --fastqc r1.fastq.gz` (the remedy the new message offers **first**) | — | `--paired with a single input file is only legal if that file is a uBAM…` |

So the user who forgot Read 2 loses a message that named their exact error, and the first of the
three remedies they are handed fails with precisely the message they used to get. They have to read
to the second clause. §2.4/O1 anticipated this user and resolved it by rewording, but the reworded
message still leads with the wrong remedy for them, and its central claim ("that shape is the
interleaved-uBAM paired path") is a statement about dispatch, not about their file.

**Fix (cheap):** reorder the remedies — `pass Read 1 and Read 2 as two files` first, `--output-format
ubam` second, qualified as being for a genuine interleaved uBAM.
**Fix (correct, declined by the plan):** site arm 5's guard in `main.rs` immediately after `:334`,
where `input_formats` is known, so the FASTQ user keeps `:334`'s message and the uBAM user gets a
refusal that can truthfully assert the format. Splitting the family across two files is a real cost;
if that call stands, the wording fix above is required.

### MEDIUM-4 (Logic/tests) — the worst mutation is invisible to the unit suite

Mutation M7: delete `&& self.input.len() == 1` from arm 5. The guard becomes
`paired && !ubam && !clump_only`, which refuses **`--paired --fastqc R1 R2`** — plain paired-end
FASTQ with FastQC, the most common FastQC invocation there is.

Result: **all 87 `cli::tests` pass.** All 14 `integration_paired_format_guard` tests pass. Only two
*pre-existing* tripwire cases go red (`paired_trim_fastqc`, `paired_trim_passthrough_fastqc`).

The change is therefore guarded, but not where the plan says its guards are: §5 step 3's acceptance
twins are all uBAM-output or `--clump_only` shapes, and both uBAM twins survive M7 because the
`output_format` term still fires. The plan calls over-rejection "the only way this change can do
harm" — the two cheapest, most-run tests for that risk are missing.

**Fix:** two unit twins, four lines total, next to the existing ones:

```rust
#[test] fn fastqc_plus_paired_fastq_output_accepted() {
    let cli = Cli::parse_from(["trim_galore", "--paired", "--fastqc", R1, R2]);
    assert!(cli.validate().is_ok(), "plain PE FASTQ + --fastqc must stay accepted");
}
#[test] fn fastqc_plus_single_end_accepted() {
    let cli = Cli::parse_from(["trim_galore", "--fastqc", R1]);
    assert!(cli.validate().is_ok(), "SE trim + --fastqc must stay accepted");
}
```

### MEDIUM-5 (Logic/tests) — `odd_input_count_precedes_the_fastqc_refusal` pins the wrong thing, and §2.4's evidence is misattributed

§2.4 row 3 says `--clock --paired --fastqc r1 r2 r3` reports "`--clock mode requires an even number
of input files` (`cli.rs:1004`)". It does not. Verified in the worktree:

```
$ trim_galore --clock --paired --fastqc r1 r2 r1
Error: Paired-end mode requires an even number of input files (R1/R2 pairs), got 3   ← cli.rs:650
$ trim_galore --clock --fastqc r1 r2 r1          # no --paired
Error: --clock mode requires an even number of input files (R1/R2 pairs), got 3      ← cli.rs:1015
```

`validate_paired_input("Paired-end")` at `cli.rs:650` fires first whenever `--paired` is set, so the
`--clock`/`--implicon` calls at `:1015`/`:1018` are reachable **only without `--paired`** — and the
test uses `--paired`. Consequence, confirmed by mutation M5: relocating the whole guard block above
`#386` (i.e. above `:997`, `:1004`, `:1015`, `:1018`) leaves this test **green**. It distinguishes
only "at the plan's old 634 site" from "anywhere below `:650`"; it says nothing about the range and
paired-input checks the placement was moved for.

`hardtrim5_range_error_precedes_the_fastqc_refusal` **does** go red under M5, so the late placement
is not unguarded. But the second precedence test is worth one word:

**Fix:** drop `--paired` from `odd_input_count_precedes_the_fastqc_refusal` and assert
`"--clock mode requires"`. That case genuinely exercises `:1015` and goes red if the block moves up.
The same applies to §2.4 row 4 (`--implicon --paired … r1.fq r1.fq`): its R1==R2 message also comes
from `:650`, not `:1018`.

### LOW-6 (Structure/docs) — "Four paths" enumerates five

`docs/…/guide/outputs.md:145`: *"**Four** paths write no report and reject the flag …: the specialty
modes (`--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`), and `--paired` given a single
interleaved uBAM with FASTQ output."* That is four specialty modes **plus** the interleaved arm =
five. Neither reading gives four. Probably a slip from "the four specialty modes", which
`CHANGELOG.md`, `flags.md` and `quickstart.md` all phrase correctly.

**Fix:** "Five paths…", or "Two situations…: the four specialty modes …, and …".

### LOW-7 (docs) — the reports-per-run sentence still omits `--passthrough`

`outputs.md:147`, as edited: *"FastQC is invoked once per trimmed primary output — so PE FASTQ
produces two reports (one per mate)…"*. Under `--passthrough` it produces **three**
(`main.rs:1795-1797`), and `docs/…/modes/passthrough.md:65` says so in as many words: "FastQC runs
on all three outputs, including the passthrough." The edit fixed two of the three exceptions the
sentence had (§5 step 8 named the interleaved arm and `--retain_unpaired`) and left this one, while
"once per trimmed primary output" arguably makes the carrier harder to account for than the old
"once per output file" did.

**Fix:** add `--passthrough`'s carrier to the same clause as the `*_unpaired_*` note.

### LOW-8 — RESOLVED mid-review (Logic) — one malformed-value check was still pre-empted, contrary to the block's own comment

> **Fixed in the working tree after my snapshot.** The block now sits below the `-a2` block
> (`cli.rs:1059`), and `malformed_adapter2_precedes_the_fastqc_refusal` pins it — verified red under
> mutation M8. The original finding is kept below for the record.


The block's comment claims it is "sited below the range and paired-input checks so a malformed value
or input list reports its own defect first". One malformed-value check remains **below** it:
`crate::adapter::parse_adapter_specs_quiet(&self.adapter2)` at `cli.rs:1098`. So
`--paired --fastqc -a2 <malformed> one.bam` now reports the FastQC refusal instead of the adapter
parse error. (For arms 1–4 the `-a2` branch only warns, and the run now fails anyway — immaterial.)
Everything else below the block is deprecation `eprintln!`s.

**Fix (optional, ~25 lines):** move the block below the `#369 -a2` block, above the deprecation
warnings. Then the stated rule holds without exception.

Related, for the plan rather than the code: §6's "before any I/O" is not accurate — the block sits
after the `check_restartable_input` loop, which opens every input (and blocks on a writer-less FIFO).
That is the right trade for message ordering, but it is not a pre-I/O check.

### LOW-9 (Structure) — three comments over the two-line max, two of them previous-state narrative

Against `~/.claude/CLAUDE.md`'s "one line default, two maximum" and "current state only — never
reference the state you are fixing":

1. ~~`src/cli.rs:1034-1036` — three lines.~~ **RESOLVED mid-review** — now two lines at
   `cli.rs:1059-1060`, current-state phrasing, no narrative.
2. `tests/integration_preflight_tripwire.rs:1317-1318` — three lines, **and** pure previous-state
   narrative: *"Replaces the `fastqc_incapable` marker (#421). That marker was latent — no case
   paired a FastQC flag with a marked arm, so it never fired. These five run."* `fastqc_incapable`
   now appears **nowhere else in the tree** (verified by grep over `src tests docs README.md
   CHANGELOG.md` — this comment is the only hit), so a future reader cannot look up what it refers
   to. That belongs in the commit message.
3. `tests/integration_clump_only.rs` doc comment, last sentence: *"Passing bare `cli.fastqc` made
   `--fastqc_args` alone a silent no-op on this path."* Previous-state narrative. Weakest of the
   three — this file's existing comments do carry issue history ("the contract #409 violated") — but
   the first two sentences already say everything a reader needs.

### LOW-10 (tests) — two assertions that cannot fail / are looser than they read

- `cli.rs`, `fastqc_plus_paired_single_file_rejected`: `!err.contains("your input is")` is vacuous —
  no message anywhere in the codebase contains that string, so this half of the assertion can never
  go red. The paired `err.contains("two files")` is the real one. Drop the dead half, or replace it
  with an assertion about remedy *order* (which is what MEDIUM-3 is about).
- `refuses_fastqc_on_every_arm_that_cannot_run_it` asserts only the shared prefix
  `"--fastqc is not supported"`, so a guard that bailed with the *wrong arm's* message passes all
  five cases. The `cli::tests` twins do check arm names, so combined coverage is adequate; closing it
  needs a per-case expected fragment (the labels don't map to flag spellings — `paired_interleaved`).

### LOW-11 (Structure) — `clump_only.rs`'s `fastqc: bool` parameter name is what caused the bug

Four signatures take `fastqc: bool` (`clump_only.rs:258`, `:395`, `:742`, `:892`), now fed
`cli.fastqc_requested()`. The sixth arm existed precisely because a parameter named `fastqc` looks
like it should be fed `cli.fastqc`. Renaming the parameter `fastqc_requested` in those four
signatures makes the next reader unable to make the same mistake. The plan deliberately declined
`clump_only.rs` signature changes (§7), so this is a judgement call, not a defect — but it is the one
change that prevents recurrence rather than fixing this instance.

### LOW-12 (docs/process) — CHANGELOG bundles two changes under one `**Breaking:**` entry

The entry pairs the refusal (breaking) with the `--fastqc_args`-under-`--clump_only` fix
(non-breaking) via a semicolon. It is one sentence and matches the neighbouring entries' shape, so it
is not a style violation — but a user scanning for "did `--fastqc_args` start working?" will not look
under a Breaking heading. Two one-sentence entries would be more discoverable. Also: the four
specialty modes counted correctly here, which is what makes LOW-6 look like a copy slip.

### LOW-13 (housekeeping, not this change)

The working tree carries untracked cruft that `git add -A` would sweep into this PR:
`10K_150bp.fastq.gz`, `10K_150bp_trimmed.fq.gz`, `examples/fastqc_only.rs`, `optimus_prime/`,
`target_repro_b/`, `SESSION_HANDOFF.md`. `examples/fastqc_only.rs` in particular is a cargo example
that `cargo clippy --all-targets` compiles locally but CI would never see.

---

## 3. What I verified rather than assumed

### 3.1 Over-rejection hunt (the only way this can do harm)

**Arms 1–4 cannot over-reject.** The four specialty dispatches are unconditional on their flag —
`main.rs:464` `if let Some(n) = cli.hardtrim5`, `:495` `if let Some(n) = cli.hardtrim3`, `:526`
`if cli.clock`, `:548` `if let Some(umi_len) = cli.implicon` — and all four sit *above*
`--clump_only` (`:580`) and the uBAM-output dispatch (`:845`). So whenever one of those flags is set,
control reaches an arm with no `fastqc::run` in it, whatever else is on the command line
(`--output-format ubam`, `--clump_only`, `--clumpify`, `--passthrough`, multi-pair input). Confirmed
`src/specialty.rs` has zero `fastqc` references.

**Arm 5's four terms exactly reproduce `:853`.** `:845` (uBAM out) and `:580` (`--clump_only`) both
return above it, which is what the `output_format` and `clump_only` terms encode; `paired && len == 1`
is `:853` verbatim. Confirmed `run_paired_ubam_single_file` (`:1803`+) contains no `fastqc::run`.

**Ten `fastqc_requested()` sites.** `grep` for `cli\.fastqc\b` over `src/` after the change returns
nothing — every read now goes through the accessor. The five gate expressions were
`cli.fastqc || cli.fastqc_args.is_some()`, which is the accessor's body, so those five are identity
substitutions. The five `--clump_only` pass-throughs feed `if fastqc { fastqc::run(...) }`, and the
disjunction is the established semantics for that predicate everywhere else (pinned pre-change by the
tripwire's `se_trim_fastqc_args_only`), so none of them wanted the narrower `cli.fastqc`.

### 3.2 Guard placement — the full enumeration asked for

**Between the plan's old candidate site (~`:634`, just after §3.4a) and the new one (`:1034`)**, in
order — every one of these would have been pre-empted at 634, and the plan names only four:

| Lines | Check |
|---|---|
| 647-650 | `validate_paired_input("Paired-end")` — even count, R1==R2, duplicate pairs |
| 653, 680 | `--basename` conflicts (SE multi-input, PE >2 inputs) |
| 662-679 | #383 duplicate-input detection (SE) |
| 686, 693, 704, 713, 722, 726 | **malformed-value checks** — `--error_rate`, `--stringency`, `--nextseq`, `--consider_already_trimmed`, `--times`, `--cores` |
| 697 | `--nextseq` + `-q` conflict |
| 730-760 | `--clumpify` exclusion matrix (incl. clock/implicon/hardtrim5/hardtrim3/demux) |
| 771-912 | `--clump_only` exclusion matrix (incl. "`--clump_only` and `--hardtrim5` are mutually exclusive") |
| 917-993 | `--passthrough` envelope (incl. clock/implicon/hardtrim), + passthrough file existence |
| 997 | hardtrim5+hardtrim3 conflict |
| 1004, 1009 | hardtrim range (`must be between 1 and 999`) |
| 1014-1018 | `validate_paired_input` for `--clock` / `--implicon` (reachable only **without** `--paired` — see MEDIUM-5) |
| 1020-1027 | `--demux` + `--paired`; barcode file existence |
| 1030-1032 | input file existence |

So the late siting is a bigger win than §2.4 claims: `--clump_only --hardtrim5 --fastqc a.fq` keeps
"mutually exclusive", `--clumpify --clock --fastqc` keeps "not yet supported with --clock",
`--cores 0 --hardtrim5 20 --fastqc` keeps "must be at least 1", and so on.

**After the new site:** only the `#369 -a2` block (LOW-8) and five deprecation `eprintln!`s. Nothing
else. No check exists that the block *should* have preceded.

### 3.3 Mutation matrix — can the new tests fail?

Every mutation was applied in the worktree, run, then reverted; the worktree was byte-compared to the
submitted patch at the end (`WORKTREE CLEAN: identical to submitted diff`).

| # | Mutation | Result |
|---|---|---|
| M1 | delete the whole guard block | **6 `cli::tests` red** (5 arms + `fastqc_args_alone_also_rejected`), `refuses_fastqc_on_every_arm…` red, `fastqc_refused_on_paired_interleaved_fastq_output` red. Both precedence tests stay green (by design — they test siting, not existence) |
| M2 | revert the 5 `--clump_only` pass-throughs to `cli.fastqc` | **`clump_only_fastqc_args_alone_produces_a_report` red** — `expected _fastqc.html + _fastqc.zip, got []`. The sixth-arm fix is self-proving |
| M3 | drop `!matches!(output_format, UBam)` from arm 5 | **`fastqc_plus_paired_ubam_output_single_file_accepted` red** *and* `fastqc_accepted_on_paired_interleaved_ubam_output` red. The acceptance twins are real |
| M4 | drop `!self.clump_only` | **`clump_only_paired_single_file_keeps_its_own_message` red.** No unit test covers this term (correctly — it is a `main.rs` message-ordering property) |
| M5 | relocate the block above `#386` (`:997`) | **`hardtrim5_range_error_precedes_the_fastqc_refusal` red**; `odd_input_count_precedes_the_fastqc_refusal` **stays green** → MEDIUM-5 |
| M6 | make the tripwire's FastQC allowance never apply | **`se_trim_fastqc`, `se_trim_fastqc_args_only`, `se_trim_fastqc_args_outdir` red** with "the allowance is masking a broken case" → the collapsed `if active && allowed.is_empty()` rule is live, and now exercised by 8 FastQC-flagged cases instead of being scoped away from 5 arms |
| M7 | drop `input.len() == 1` (refuses plain `--paired --fastqc R1 R2`) | **87/87 `cli::tests` green**, 14/14 PFG green; only pre-existing `paired_trim_fastqc` + `paired_trim_passthrough_fastqc` red → MEDIUM-4 |
| M8 | *(against the post-snapshot tree)* move the block back above the `-a2` block | **`malformed_adapter2_precedes_the_fastqc_refusal` red.** The mid-review relocation is genuinely pinned |

### 3.4 Tripwire deletion

Deletion is complete: `fastqc_capable` (field, builder, 6th parameter, `:311` argument) and all nine
literal `true`s are gone; `grep` finds no surviving reference except the comment in LOW-9. All
**54** tests in the file pass. The new case does assert what it claims — non-zero exit
(`assert!(!ok)`), the message, and nothing written (`snapshot().difference(&before)` empty, with the
fixtures staged *before* the snapshot so they are not counted). It correctly uses the file's
`copy_fixture` helper (iteration-log #3) rather than a bare `fs::copy`. Caveats in LOW-10.

The module-doc rewrite (`:34-37`) is accurate and consistent with the collapsed rule; the older
paragraph at `:25-28` ("only on arms that pass a FastQC flag, and on those arms the allowed set must
be non-empty") also still reads true.

### 3.5 Docs claims, checked one by one

| Claim | Verdict |
|---|---|
| `--retain_unpaired`'s `*_unpaired_{1,2}` files get no report | **True** — `main.rs:1777-1798` runs FastQC on `output_r1`, `output_r2` and the `--passthrough` carrier only |
| "Four paths write no report" | **False** — five (LOW-6) |
| "PE FASTQ produces two reports (one per mate)" | **Incomplete** — three under `--passthrough` (LOW-7) |
| `flags.md`: "rejected on the four specialty modes and on `--paired` given a single interleaved uBAM with FASTQ output … adding `--output-format ubam` does produce a report" | **True**, and correctly scopes the remedy to the last case — this is the wording `--help` should copy (HIGH-2) |
| `quickstart.md`: "rejected on the specialty modes and on `--paired` with a single interleaved uBAM writing FASTQ" | **True** |
| `CHANGELOG.md`: "the four specialty modes …; `--fastqc_args` alone activate[s] FastQC under `--clump_only`" | **True** (see LOW-12 for the bundling nit) |
| `outputs.md:143` blanket "Works on both FASTQ output … and uBAM output" left unedited | **Acceptable** — the exceptions paragraph now immediately follows it; §5 step 8's concern is met in substance |
| Any doc promising FastQC on hardtrim/clock/implicon | **None** — `modes/hardtrim.md`, `clock.md`, `implicon.md` contain no `fastqc` reference at all, so §2.6's sweep holds |

Efficiency: nothing to report. Five boolean tests plus one `Option::is_some()` on the
argument-validation path; all messages are `&'static str` formatted only on the error path; the
accessor is a trivial disjunction.

---

## 4. Recommendations, in order

**Before merge**

1. **HIGH-1** — move the new tripwire test below `a_refused_run_writes_nothing` (or re-pair the doc
   comments). One-line edit, zero risk, and it is the sort of thing that never gets fixed later.
2. **HIGH-2** — scope `--help`'s `--output-format ubam` remedy to the arm it works on. Copy
   `flags.md`'s phrasing, which already gets this right.
3. **MEDIUM-3** — reorder arm 5's three remedies so "pass Read 1 and Read 2 as two files" comes
   first. If relocating the guard to `main.rs:334`+ is on the table instead, that is strictly better
   for the user and worse for cohesion; the plan's call to keep it in `validate` stands, but then the
   wording must not lead with a remedy that fails for the likelier user.
4. **MEDIUM-4** — add the two acceptance twins (`--paired --fastqc R1 R2`, `--fastqc R1`). Four
   lines; closes the only mutation that slipped past the unit suite.

**Worth doing in the same PR**

5. **MEDIUM-5** — drop `--paired` from `odd_input_count_precedes_the_fastqc_refusal` so it pins
   `cli.rs:1015`; correct §2.4 rows 3–4, which attribute both messages to the wrong call site.
6. **LOW-6, LOW-7** — the count and the `--passthrough` clause in `outputs.md`.
7. **LOW-9** — cut the previous-state narrative from the tripwire comment, which names an identifier
   that no longer exists anywhere in the tree. (The `cli.rs` half is already done.)
8. **LOW-10** — delete the vacuous `!err.contains("your input is")` assertion.

**Optional / judgement**

9. ~~**LOW-8**~~ — done mid-review.
10. **LOW-11** — rename `clump_only.rs`'s `fastqc: bool` parameter to `fastqc_requested`; this is the
    only change that prevents the sixth arm from coming back.
11. **LOW-12** — split the CHANGELOG entry in two.
12. **LOW-13** — check `git status` before staging; six untracked paths would otherwise ride along.
13. Consider one line on `modes/hardtrim.md` / `clock.md` / `implicon.md` noting the refusal — those
    pages are where a user of those modes looks, and they currently say nothing about FastQC either
    way.

Nothing here blocks the behaviour change: the guards are correct, the sixth-arm fix is self-proving,
and the tripwire is stronger after the deletion than before it.

---

**Report:** `/Users/fkrueger/Github/TrimGalore/plans/08122026_fastqc-silent-noop/CODE_REVIEW_B_gen2.md`
**Worktree:** `$TMPDIR/cr-B2-421` (detached @ `3ffcef4` + patch; restored, byte-identical to the
submitted diff — safe to `git worktree remove`)
