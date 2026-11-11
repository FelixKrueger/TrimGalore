# PLAN v2 — #441 follow-up: refuse `--library` where it cannot apply, and close the test seams

**Feature dir:** `plans/08162026_pr441-library-presets/` (alongside two code reviews and two plan reviews)
**Base:** `dev` after PR #441 merges. Anchors below are quoted at `145432e` and verified there by both
plan reviewers (all 14 resolve correctly). See §10 R-4 for why #439 merging first does not shift them.
**Status:** revised after dual plan review. Awaiting the implementation trigger.

## Revision history

**v2 (2026-08-16)** — revised after `PLAN_REVIEW_reviewer-A_g1.md` and `PLAN_REVIEW_reviewer-B_g1.md`.
Three Critical items from A, all confirmed by B, plus Felix's two decisions:

| Change | Source |
|---|---|
| Refusal wording is now **"does not honour the clip flags"** | Both reviewers; Felix's decision. v1's "does not clip" was **false** — `specialty.rs:107,219` call `record.clip_5prime()` and `:241,:352` document clipping |
| #441's CHANGELOG entry is **edited in place**, no second entry added | A's C-1 + Felix's decision. v1 added a new entry while leaving a falsified sentence standing |
| Stderr log line gains test coverage | A's C-2 / B's #3 — v1 left half its own headline fix silently revertible |
| `--implicon` test shape corrected; assumption 6 deleted | A's C-3 / B's #2 — v1's `--implicon 8` does not parse |
| **Both** override-cardinality tests kept (1 and 4) | B's #5 — v1's step 6 destroyed the only cardinality-1 coverage |
| V5 and V13 replaced; every V-row must name a test that runs | A's I-4/O-4 — both were vacuous |
| The asymmetry argument replaced with the reviewers' stronger one | A's I-2 / B — both independently found it |
| `json["mode"]`, not a `paired` key | A's O-1 — v1 cited a key that does not exist |
| rustfmt note added to step 2 | A's I-7 |
| §11 added: disposition for every omitted finding | A's I-6 |

---

## Goal

Close the two behavioural findings both code reviewers of #441 agreed on, plus the test seams they
found between them, in a follow-up PR of ours rather than a review round on the contributor's fork.

1. `--library` is **refused** on the four specialty modes (`--hardtrim5`, `--hardtrim3`, `--clock`,
   `--implicon`), which return before any clip site is reached.
2. The report's `Clipping in force:` line **omits the R2 values on single-end runs**, where they
   cannot apply.

Plus test hardening so each seam has a test that provably fails without the code under it.

## Context

### Placement and delivery

- **#441 merges first, unmodified.** Spec-conformant on all five agreed clauses, 9/9 CI green at
  `145432e`, `--rrbs` byte-identity preserved (two independent base-commit differential builds). The
  contributor is a `FIRST_TIME_CONTRIBUTOR`; their commit lands as they wrote it.
- **This is a separate PR off `dev`**, branch `fix/440-library-refusals`.
- **No release-timing risk.** The release is parked, so #441 and this follow-up ship together. The
  refusal lands *before* `--library` is ever published, so it is not a breaking change — tightening
  after a release would be.

### Decisions taken (Felix)

| Decision | Choice |
|---|---|
| Delivery | Merge #441, then our follow-up PR |
| Scope | Substantive findings + test hardening. **Editorial excluded** |
| Specialty modes | **Refuse** in `Cli::validate()`, not warn |
| SE reporting line | **Gate** the R2 halves on `paired`, not reword |
| Refusal wording | `"<mode> does not honour the clip flags; --library is not compatible"` |
| CHANGELOG | Edit #441's entry **in place**; add no second entry |

### Existing code this must align with

- **`cli.rs:1071-1097`** — the #421 `--fastqc` refusal on the same four modes. The structural template.
- **`cli.rs:2313-2351`** — that block's four `fastqc_plus_<mode>_rejected` tests (note: `-2351`, not
  `-2350` as v1 said).
- **`cli.rs:2466`** — `hardtrim5_range_error_precedes_the_fastqc_refusal`, which pins validation
  *ordering*.
- **`cli.rs:833-850`** — the `--clump_only` family, whose `<mode> does not X; --library is not
  compatible` shape we keep.
- **`cli.rs:1050-1055`** — the `-a2` block, which already phrases per-mode incapability correctly
  ("`--clock` performs no adapter trimming"). Our wording is the clip-flag analogue.
- **`main.rs:273-277`** — the `--phred64` precedent: uniform refusal over a mode-dependent rule,
  accepting the removal of working invocations.

### The load-bearing detail

`flag_summary()` has **exactly two** call sites — confirmed by both reviewers, three grep hits total
including the definition, and **no test calls it**:

- `main.rs:1121` — the stderr log line (`cli.paired` in scope)
- `report.rs:279` — the text report (`config.paired` in scope)

Gating one without the other desyncs the log from the report for the same run. v1 got this right and
then **tested only one of them**; §V6a now covers the stderr half.

### Established empirically by Reviewer A, before implementation

Applying steps 1-2 verbatim and running `cargo test --no-fail-fast` gives **751 tests, exactly one
failure** — `both_reports_record_the_preset_and_its_expanded_values` — emitting byte-for-byte:

```
Library preset: pbat
Clipping in force: --clip_R1 8 --three_prime_clip_R1 8
```

The blast radius is one test, measured rather than asserted. `cargo clippy -D warnings` was clean on
the patched tree.

## Behavior

### B1 — `--library` with a specialty mode

1. `Cli::validate()` reaches the new block at **line 1111** — the first statement after the
   `fastqc_requested` block closes at `1110`. `1103` is the last `bail!` in the file and everything
   from `1110` to `Ok(())` at `1140` is `eprintln!` warnings, so **our block becomes the last
   fallible check** and every precedence guarantee follows automatically.
2. If `self.library.is_some()` and any of `hardtrim5` / `hardtrim3` / `clock` / `implicon` is set,
   bail naming the mode and `--library`.
3. Exit is non-zero and **no output file is created** — `validate()` precedes `ensure_output_dir` and
   every dispatch arm.

Edge cases, all verified by both reviewers:
- `--hardtrim5 0` reports its `1 and 999` range error first (`cli.rs:1019`).
- `--hardtrim5` + `--hardtrim3` reports #386's bail first (`cli.rs:1009`).
- `--clock` / `--implicon` arity errors report first (`cli.rs:1027`/`:1030`).
- Missing input files report first (`cli.rs:1042-1044`) — a case v1 did not mention.
- **`--clump_only` short-circuits before the new block entirely.** `--clump_only --hardtrim5
  --library emseq` yields `--clump_only does not clip; --library is not compatible` from `cli.rs:850`,
  because the clip family precedes the specialty mutual-exclusions at `:893-903` and the whole
  `clump_only` block at `:780` precedes ours. **Our block is unreachable on clump-only runs.** Correct
  and pre-existing; stated here so B1 does not read as if it always fires.

### B2 — `Clipping in force:` on single-end

1. `flag_summary` takes `paired: bool`.
2. R2 entries are skipped when `!paired`; R1 entries unaffected in both modes.
3. Paired output is byte-identical to today. Single-end becomes `--clip_R1 8 --three_prime_clip_R1 8`
   for `pbat`.

Edge cases:
- **Override lines are deliberately NOT gated.** `--library emseq --clip_R2 5` on a single-end run
  still logs `--clip_R2 5 was given on the command line and overrides the emseq preset value 10`.
  That sentence is true about what the user typed and claims nothing about application.
  `Clipping in force:` does make that claim, which is why it alone is gated.
- **The JSON keeps all four values on single-end.** See §8.3 for the decisive reason.

## Signature

```rust
/// The clipping in force, rendered as the command line that would reproduce it.
/// R2 values are listed for paired runs only.
pub fn flag_summary(&self, paired: bool) -> String
```

## Implementation outline

### 1. `src/library.rs` — gate the R2 halves

Replace `flag_summary` (`:121-134`):

```rust
/// The clipping in force, rendered as the command line that would reproduce it.
/// R2 values are listed for paired runs only.
pub fn flag_summary(&self, paired: bool) -> String {
    let mut parts = Vec::new();
    for (flag, value, read_2) in [
        ("--clip_R1", self.clip_r1, false),
        ("--clip_R2", self.clip_r2, true),
        ("--three_prime_clip_R1", self.three_prime_clip_r1, false),
        ("--three_prime_clip_R2", self.three_prime_clip_r2, true),
    ] {
        if read_2 && !paired {
            continue;
        }
        if let Some(v) = value {
            parts.push(format!("{flag} {v}"));
        }
    }
    parts.join(" ")
}
```

Delete the two justification sentences from the existing doc comment (they restate design rationale;
CLAUDE.md puts that in the commit message). Both reviewers noted this incidentally satisfies part of
B's code-review C-1.

### 2. Both call sites — and reformat

- `src/main.rs:1121` → `clips.flag_summary(cli.paired)`
- `src/report.rs:279` → `library.flag_summary(config.paired)`

**`report.rs:279` must then be reformatted**: the added argument pushes it past rustfmt's width and
rustfmt splits the line. Run `cargo fmt` after this step, or CI's `cargo fmt --all -- --check` fails
on a change nobody would think to re-format (A's I-7).

### 3. `src/cli.rs:1111` — refuse `--library` on the four specialty modes

Insert immediately after the `if self.fastqc_requested() { … }` block closes at `1110`:

```rust
// The specialty arms return before setup_trimming, so no preset reaches a clip site.
if self.library.is_some() {
    if self.hardtrim5.is_some() {
        anyhow::bail!("--hardtrim5 does not honour the clip flags; --library is not compatible");
    }
    if self.hardtrim3.is_some() {
        anyhow::bail!("--hardtrim3 does not honour the clip flags; --library is not compatible");
    }
    if self.clock {
        anyhow::bail!("--clock does not honour the clip flags; --library is not compatible");
    }
    if self.implicon.is_some() {
        anyhow::bail!("--implicon does not honour the clip flags; --library is not compatible");
    }
}
```

**Wording rationale (do not revert to "does not clip").** All four modes *do* perform end clipping:
`specialty.rs:107` and `:219` call `record.clip_5prime(start)`, `:241` documents `--clock` as "then
clips: R1 at position 13 …, R2 at position 15", and `:352` documents `--implicon` as "then clips R2 by
`umi_length` bases". What they do not do is honour `--clip_R1`/`--clip_R2`/`--three_prime_clip_*` —
`specialty.rs` contains zero references to any of them. "Does not clip" would contradict the tool's
own documentation in a user-facing error.

### 4. `src/cli.rs` tests — new section after the #421 block's tests

```
// ── #440: --library refused where no clip site is reachable ──
```

1. `library_plus_hardtrim5_rejected` — `["trim_galore", "--hardtrim5", "20", "--library", "emseq", R1]`
2. `library_plus_hardtrim3_rejected` — same with `--hardtrim3`
3. `library_plus_clock_rejected` — `["trim_galore", "--clock", "--paired", "--library", "pbat", R1, R2]`
4. `library_plus_implicon_rejected` — **`["trim_galore", "--implicon", "--paired", "--library", "pbat", R1, R2]`**

   **Bare `--implicon`, no separate value.** `cli.rs:419-421` sets `require_equals = true` with
   `default_missing_value = "8"`, so `--implicon 8` leaves the `8` in `input` and produces a
   misleading "requires an even number of input files" error. Reviewer B proved this by running it.
   Mirrors the existing `fastqc_plus_implicon_rejected` at `cli.rs:2345`. (`--implicon=8` also works.)

Each asserts the error contains both `--library` and the mode name.

5. `hardtrim5_range_error_precedes_the_library_refusal` —
   `["trim_galore", "--hardtrim5", "0", "--library", "emseq", R1]` must report the `1 and 999` range
   error, **not** the refusal. Mirrors `cli.rs:2466`.

### 5. `src/report.rs` tests

6. **Keep `test_json_report_records_library_block` as it is** — one override, `emseq`, `clip_r1: 12`.
   It is the **only** coverage of override-cardinality 1, which is the branch where `report.rs:894`'s
   `if i + 1 < len` must *not* emit a comma (B's #5). v1 would have deleted it.

7. **Add `test_json_report_records_all_four_overrides`** — `accel` with all four flags overridden to
   **distinct** values:

   ```rust
   config.library = Some(crate::library::resolve(
       Some(crate::library::LibraryPreset::Accel),
       crate::library::UserClips {
           clip_r1: Some(1), clip_r2: Some(2),
           three_prime_clip_r1: Some(3), three_prime_clip_r2: Some(4),
       },
   ));
   ```

   Assert the block's four keys are `1, 2, 3, 4` respectively — four *distinct* values is what kills
   key permutation, since `accel` alone (10/15/10/10) catches only swaps touching `clip_r2` and `pbat`
   (8/8/8/8) catches none. **This is the T-1 fix.** In the same test assert the four `ClipOverride`
   records carry the right flag names and the right displaced preset values (`10, 15, 10, 10`) — the
   T-2 fix — which also exercises the comma path at cardinality 4.

   Reviewer A confirmed the mechanism: `resolve()` evaluates `pick` in struct-literal field order, so
   `overrides` is deterministically `[--clip_R1(10→1), --clip_R2(15→2), --three_prime_clip_R1(10→3),
   --three_prime_clip_R2(10→4)]`. Resolved and displaced values are independently assertable.

8. `test_text_report_omits_r2_clipping_on_single_end` — new. `config.paired = false`, `pbat`, assert
   the text contains `Clipping in force: --clip_R1 8 --three_prime_clip_R1 8` and **not** `--clip_R2`.

9. `test_text_report_records_preset_name_and_expanded_values` already sets `config.paired = true`, so
   it stays green unchanged and becomes the paired half of the pair.

10. `test_json_report_library_is_null_without_a_preset` — add a **presence** assertion beside the
    existing null check (F-5a): `assert!(json["parameters"].as_object().unwrap().contains_key("library"))`.
    `serde_json` yields `Value::Null` for a *missing* key, so `is_null()` alone cannot distinguish
    absent from present-and-null.

11. Add two direct unit tests for `flag_summary` itself (B) — one paired, one single-end. It currently
    has none, and they are the cheapest guard for B2.

### 6. `tests/integration_library.rs`

12. **`both_reports_record_the_preset_and_its_expanded_values` must be updated.** It is single-end
    (`--dont_gzip --library pbat -o dir R1`) and asserts all four values, which the gate falsifies.
    Change the expectation to `--clip_R1 8 --three_prime_clip_R1 8`. Reviewer A proved this is the
    **only** test in 751 that the gate breaks.

13. Add the **F-5b agreement loop** to the same test — for all four keys, assert
    `json["parameters"][k] == json["parameters"]["library"][k]`.

14. **`explicit_clip_flag_overrides_the_preset_and_is_logged` gains a stderr assertion** —
    `assert!(!err.contains("--clip_R2"))`. This is the **C-2 fix** and the most important single
    addition in v2: without it the `main.rs:1121` half of B2 is silently revertible, because the
    existing `err.contains("--clip_R1 12")` is satisfied by the *override* line regardless of whether
    the banner was gated.

15. `a_run_without_a_preset_reports_no_library_block` — add the same `contains_key` presence assertion.

16. New `paired_report_lists_all_four_values` — a `--paired --library accel` run asserting the text
    report carries all four **and** that the stderr banner does too, so an over-eager gate cannot
    silently drop R2 from paired runs.

17. New `preset_with_multiple_overrides_round_trips_through_json` — real binary, `--library accel` with
    all four clip flags overridden; assert the JSON parses and all four `overrides` entries are present.

### 7. `tests/integration_ubam_out.rs` — T-3

18. New test: `--rename --output-format ubam --library emseq` on FASTQ input is refused. **Assert the
    error text, not merely a non-zero exit** (A's I-5), and assert the whole output listing is empty
    following `integration_ubam_out.rs:1032-1038`'s pattern (B). Reviewer A confirmed `main.rs:391` is
    the only `--rename`+uBAM guard, so reverting it makes the run succeed — V10 is not vacuous.

### 8. Docs — factual correctness only

`docs/src/content/docs/guide/library-presets.md`:

19. The override transcript at `:48-54` is single-end and shows all four flags. Replace with the real
    two-flag output.
20. The text-report transcript at `:62-66` is single-end and shows all four. Same fix.
21. **`:60` prose is falsified** — it says "Both trimming reports record the preset name **and the
    four values** that were used". Rewrite, and site the text/JSON explanation **at 60-79**, not only
    in Overrides: that passage shows the text report and the JSON block for the *same* single-end run
    seven lines apart, so a sentence placed elsewhere reads as contradicted by the code block below it
    (A's I-3, B's docs note).
22. §"What a preset does *not* set" gains one bullet: the four specialty modes refuse `--library`,
    because they do not honour the clip flags.
23. **Decide** whether `flags.md:14` gains a refusal clause (A's O-5). It is not falsified by B2 — it
      says `--library` "expands into the four clip flags", which stays true. Recommendation: add a
      half-sentence, since it is the flag-interaction page.

### 9. `CHANGELOG.md` — edit #441's entry in place, add nothing

**Per Felix's decision, and A's C-1: no second entry.** A reader should not have to reconcile two
bullets about a flag that was never released.

Two edits inside #441's existing entry:

- The sentence at `CHANGELOG.md:21-22` — "Both trimming reports record the preset name and **all four
  values** in force" — is **falsified** by B2. Rewrite so the text report's behaviour is stated
  correctly: all four in JSON, and in the text report the values that apply to the run.
- Fold the refusal into the paragraph that already explains the membership criterion (the "A preset is
  only a clip bundle …" sentence), as one clause naming the four modes.

Side benefit: because we **modify** existing lines rather than **insert** at the top of
`#### Changes`, this reduces the conflict surface with #439, which inserts there. See R-4.

## Efficiency

No hot path. `flag_summary` runs once per report and gains one `bool` comparison per entry; the array
becomes a 3-tuple, still four stack-allocated iterations. The `validate()` block is four
`Option::is_some()` checks once per invocation, short-circuited on `self.library.is_some()`. Both
reviewers agreed and moved on.

## Integration

**Reads:** `Cli::library`, `hardtrim5`, `hardtrim3`, `clock`, `implicon`; `ResolvedClips`;
`TrimConfig::paired`.
**Writes:** stderr log line (`main.rs`), text trimming report (`report.rs`). **JSON untouched.**

**Order:** the refusal is in `Cli::validate()`, before `ensure_output_dir` and every dispatch arm, so a
refused run creates nothing.

**All three report sites carry `paired`** (A's claim 9): `main.rs:1449` (`run_single_file`,
`paired: false`), `:2239` (`run_ubam_output_single`, `paired: false`), `:2618` (`write_paired_reports`,
`paired: true`), each setting `library: Some(cli.effective_clips())`. So `config.paired` reaches the
gate on every one.

**Downstream:**
- `--rrbs` byte-identity: untouched. No `--library` means `config.library.preset` is `None`, so the
  preset block is never written.
- MultiQC/JSON consumers: unaffected. No JSON key added, removed, or changed.
- Invocations removed: `--library` plus a specialty mode, which today exits 0 having ignored the
  preset. Same class the `--phred64` guard accepted.
- **No existing test breaks except the one predicted.** A searched suite-wide: no test anywhere
  combines `--library` with `--hardtrim5/3`, `--clock` or `--implicon`.

## Assumptions

1. #441 merges first and unmodified.
2. **The four specialty modes do not honour the clip flags** — `specialty.rs` has zero references to
   `clip_r1`/`clip_r2`/`three_prime_clip*`/`effective_clips`. They *do* clip on their own terms, which
   is why the wording is about the flags (§3).
3. **The text report and JSON are permitted to differ on single-end.** The decisive reason, found
   independently by *both* plan reviewers and absent from v1: **gating the JSON would break the F-5b
   agreement invariant that step 13 introduces.** The top-level `clip_r2` cannot change without a JSON
   schema break, so if the `library` block dropped R2 on single-end the two would disagree by
   construction. This makes the asymmetry close to a constraint rather than a preference. A machine
   consumer distinguishes the cases with top-level **`json["mode"]`** (`"single-end"` / `"paired-end"`,
   `report.rs:803-813`) — **there is no `paired` key**; v1 cited one that does not exist.
4. Override log lines stay ungated (B2 edge cases).
5. Four distinct override values (1/2/3/4) are acceptable test data even though no real preset has
   four distinct values — the test targets key/label wiring, not a kit. Named as such in its doc
   comment.
6. Step 13's agreement invariant assumes `PresetClips`'s fields stay plain `usize` (B). If one ever
   becomes an `Option`, the RRBS auto-clip at `main.rs:1136` could fire under a preset and the two
   JSON locations could diverge.
7. Docs changes are limited to what our change falsifies plus one factual bullet.
8. `--clump_only` short-circuits before the new block, so the refusal is unreachable there (B1).

*(v1's assumption 6 about `--implicon` is deleted — it is now settled fact in step 4.)*

## Validation

Every new test must be **proved able to fail**. **Each row names a test that demonstrably runs** —
v1 had two rows that could not fail (A's I-4, O-4).

| # | Verify | How | Expected |
|---|---|---|---|
| V1 | The refusal fires on all four modes | `cargo test library_plus_` | 4 pass |
| V2 | Those tests can fail | Delete the `if self.library.is_some()` block | all 4 red |
| V3 | Ordering is pinned | Move the new block above `cli.rs:1016` | `hardtrim5_range_error_precedes_the_library_refusal` red |
| V4 | No output on refusal | Binary in an empty dir, `--hardtrim5 30 --library emseq -o .` | exit≠0, `ls` empty |
| V5 | **The SE gate works** (replaces v1's vacuous `cargo test flag_summary`) | `cargo test test_text_report_omits_r2_clipping_on_single_end flag_summary_` | pass, and the names match real tests |
| V6 | The SE gate can fail | Revert `flag_summary` to ignore `paired` | step 8's test + integration step 12 red |
| **V6a** | **The stderr half can fail** | Revert **only** `main.rs:1121` to `flag_summary(true)`-equivalent, leaving `report.rs` gated | step 14's `assert!(!err.contains("--clip_R2"))` red |
| V7 | The paired half can fail | `if read_2` without `&& !paired` | step 16 + `report.rs` paired test red |
| V8 | T-1 closed | Swap the `clip_r2` / `three_prime_clip_r1` **key labels** | step 7 red (today's tests stay green — A proved this) |
| V9 | T-2 closed | Swap the `--clip_R2` / `--three_prime_clip_R1` **flag labels** in `pick` | step 7 red |
| V10 | T-3 closed | Revert `main.rs:391` to the explicit `cli.clip_*` list | step 18 red |
| V11 | F-5a closed | Delete the JSON `library` key, both arms | steps 10 and 15 red (A proved both stay green today) |
| V12 | F-5b closed | +100 to the JSON top-level `three_prime_clip_r2` | step 13 red |
| V13 | **Cardinality 1 still covered** (replaces v1's near-tautological `--rrbs` row) | Remove the `if i + 1 < len` guard so a trailing comma is always emitted | step 6 red (cardinality-1 JSON becomes invalid) |
| V14 | Suite and lints | `cargo test`, `cargo clippy --all-targets --release -- -D warnings`, `cargo fmt --all -- --check` | 751 + ~13 new, green; **fmt clean only if step 2's reformat was done** |
| V15 | Docs transcripts are real | Run the page's commands; assert **every quoted line appears verbatim** in the output | all present |
| V16 | `--rrbs` unaffected | Base-vs-new differential over `--rrbs --paired` (the shape whose `.txt` must stay byte-identical), `-o` a fresh dir, **counting files compared** | every FASTQ and `.txt` identical |

**V8-V12 must be run against the *new* tests** (A's §7): each mutation's target test is the one added
or amended in steps 6-18, not v1's originals.

**V15** was "byte-identical" in v1, which could not pass — the transcripts are excerpts and
`sample.fq.gz` is not a fixture (B).

**V16** note: pass `-o .` from a scratch directory and count the files compared. A comparison script
once silently overwrote tracked fixtures in `test_files/` and reported a green `0 differing` over
`0 files`. V16 is weaker evidence than it looks (A's O-4: no `--library` means the new line is never
written, so this is near-tautological) — it is retained as a cheap regression guard, not as the
argument that `--rrbs` is safe.

**V14** note: never gate on one test number — macOS and Linux differ by two.

## Questions or ambiguities

**Open (assumption taken, no blocker):**

- **O-1** — Should the JSON also omit R2 on single-end? **No** — §8.3's invariant argument now makes
  this close to settled rather than a preference.
- **O-2** — `--clump_only --library` reports `"--clump_only does not clip; …"` while specialty modes
  report `"… does not honour the clip flags; …"`. Different families, both accurate; no change.
- **O-3** — Whether the excluded editorial items become an issue or stay plan notes. Felix's filing
  bar sends review spin-offs to plan notes, so they stay in §11.
- **O-4** — `flags.md:14` refusal clause (step 23). Recommendation: add it.

**Critical:** none.

## §11 — Disposition of every finding NOT fixed here

A's I-6: seven code-review findings were silently omitted in v1. Recording a decision for each, so
absence is a choice rather than an oversight.

| Finding | Source | Disposition |
|---|---|---|
| F-5c — no `--library` differential test on the uBAM / interleaved paths | code review A | **Out.** All four paths verified APPLIED by both code reviewers via record-body md5; step 18 pins the one guard that could regress. Revisit if `setup_trimming` stops being the single funnel |
| C-1 — source comments justify the change | code review B | **Partly done.** Step 1 removes the two instances in `flag_summary`'s doc comment. The rest is editorial, excluded by scope |
| L-1 — text report records clip values only when a preset was used | code review B | **Out.** A real asymmetry (explicit `--clip_R1` gets no text line) but a pre-existing one this change does not worsen. Worth its own issue if a user asks |
| F-6 — preset names are case-sensitive | code review A | **Out.** `--library XGEN` is refused while `--clip_R1`/`--clip_r1` both work. One `ignore_case = true` would fix it; not in scope |
| F-9 — preset log line repeats per input file, after the adapter scan | code review A | **Out.** Consistent with the adjacent `--rrbs` messages |
| F-7 — `json_null()` exists but `"library": null` is hand-written | code review A | **Out.** Correct as written; consistency only |
| L-2 — `match … and_then` where the file's idiom is a let-chain | code review B | **Out.** Style |
| L-4 — a failing `Vec<u8>` byte comparison dumps ~39.7 MB | code review B | **Out**, but noted for whoever next edits those assertions — our new tests compare report text, not read bytes |

## Self-Review

**v1's own catches, both confirmed by the reviewers as real:** the two `flag_summary` call sites, and
`both_reports_record_the_preset_and_its_expanded_values` being single-end.

**What v2 fixed after review:**

- **A false statement in a user-facing error.** v1's "does not clip" contradicted `specialty.rs`'s own
  code and docs. The lesson generalises: the code review's md5 evidence proved *"the clip flags are
  inert here"*, and v1 wrote down the stronger, different claim *"this mode does not clip"*. Evidence
  for a narrow proposition was silently promoted to a broad one.
- **Half of v1's headline fix was untestable.** v1 named the two call sites its load-bearing detail,
  then covered one and left V7 as a manual eyeball. Naming a risk is not the same as testing it.
- **v1 destroyed coverage while adding it** (step 6) and **shipped two validations that could not
  fail** (V5, V13). All three are the same error: a validation asserted rather than executed.
- **v1 cited a JSON key that does not exist.**

**Remaining risks:**

- **R-1** — Removing working invocations (`--library` + specialty mode). Mitigated: the flag is
  unreleased; the `--phred64` guard set the precedent.
- **R-2** — The text/JSON asymmetry. Now argued from the F-5b invariant rather than from preference.
- **R-3** — Step 7's four-distinct-values data matches no real kit; named as wiring-only in its doc
  comment.
- **R-4** — **`CHANGELOG.md` is a three-way convergence.** #439 (`fix/439-clumpify-peak-memory` @
  `70e75e3`, one ahead of `dev`), #441, and this work all touch the same `### Unreleased` /
  `#### Changes` block. §9's decision to *edit in place* rather than *insert* reduces but does not
  remove the overlap. On any conflict, resolve by grepping for every issue number that should survive
  and confirming the diff has no `^-` lines — no test reads `CHANGELOG.md`, so a dropped entry is
  invisible to CI. **The `main.rs` anchors are safe**: #439's `main.rs` hunk is `@@ -200,10 +200,10 @@`,
  a net-zero replacement, so every line number below it is preserved.
- **R-5** — Two of four agents died mid-task this session (host sleeping). Any future delegated step
  should verify its own artifacts with a command rather than trusting a completion claim.

---

## §12 — Implementation notes (2026-08-17)

**Branch:** `fix/440-library-refusals`, worktree `/private/tmp/claude-501/impl440`, based on `dev` @
`872c4ae` (#441 merged). **Uncommitted**; patch preserved at
`scratchpad/plan-artifacts-snapshot/impl440_uncommitted.patch` and verified to apply cleanly.

**Result:** 9 files, +382 / −25. **764 tests pass** (752 base + **12 new**), 0 failures across 17
binaries. `cargo fmt --all -- --check` clean, `cargo clippy --all-targets --release -- -D warnings`
clean.

### Every plan step landed

Steps 1–23 all implemented. `#441` merged first as agreed: rebased `145432e` onto `dev`, resolved the
`CHANGELOG.md` conflict by keeping both entries (verified: all 51 pre-existing refs survive, zero `^-`
lines), force-pushed to the contributor's branch, CI 9/9 by name, squash-merged as `872c4ae` with
`Benjamin Demaille` preserved as author.

### Mutation matrix — every new test proved able to fail

| V | Mutation | Tests that went RED |
|---|---|---|
| V2 | refusal block disabled | all 4 `library_plus_*_rejected` |
| V3 | `--hardtrim5` range bail suppressed | `hardtrim5_range_error_precedes_the_library_refusal` (+ its #421 sibling) |
| V6 | `flag_summary` ignores `paired` | 4: both new unit tests + `both_reports_record_…` + `explicit_clip_flag_overrides_…` |
| **V6a** | **only `main.rs:1121` reverted** | **exactly 1: `explicit_clip_flag_overrides_the_preset_and_is_logged`** — the C-2 gap is closed |
| V7 | R2 dropped unconditionally | 3, incl. `paired_report_lists_all_four_values` |
| V8 | JSON `clip_r2`/`three_prime_clip_r1` **key** labels swapped | `test_json_report_records_all_four_overrides` — T-1 closed (this mutation survived all 751 tests before) |
| V9 | `--clip_R2`/`--three_prime_clip_R1` **flag** labels swapped in `pick` | 2 — T-2 closed |
| V10 | #408 guard reverted to the explicit clip list | `rename_into_ubam_refused_for_fastq_input_via_library_preset` — T-3 closed |
| V11 | JSON `library` key deleted | both presence tests — F-5a closed (both stayed green before) |
| V12 | top-level `three_prime_clip_r2` +100 | `both_reports_record_…` agreement loop — F-5b closed |
| V13 | overrides comma guard removed | 3, incl. `test_json_report_records_library_block` — **vindicates keeping the cardinality-1 test** |

All four mutated files restored and `cmp`-verified byte-identical after each run.

**V4** — all four modes: exit 1, correct message, **0 files written**.
**V15** — every line quoted in `library-presets.md` reproduced verbatim; the paired equivalence claim
holds (2/2 trimmed FASTQ byte-identical); 0 occurrences of `--clip_R2` in a single-end run's stderr or
text report.
**V16** — base-vs-new over 5 shapes (24 files compared, none produced zero output): every FASTQ and
every `.txt` report byte-identical; JSON identical once `command_line` is excluded, which records the
two binary paths.

### Deviations

- **D1** — rustfmt reflowed three places the plan did not anticipate: the four `bail!` strings, one
  test array in `cli.rs`, and one `assert_eq!` in `library.rs`. A's I-7 predicted `report.rs`; the
  class was wider. No behaviour change.
- **D2** — step 23 (`flags.md` refusal clause) was a "decide" item; implemented per the plan's own
  recommendation.
- **D3** — V3 was run as an equivalent mutation (suppress the range `bail!`) rather than physically
  relocating the refusal block. The property proved is identical: the range error must be reported
  first.
- **D4** — V16 required excluding `command_line` from the JSON comparison. Same harness artifact
  Reviewer A hit in its own differential.

### Iteration log

1. `cargo fmt` dirty in 3 places after the code edits (D1). Ran `cargo fmt --all`; 764 still green.
2. Mutation runs with several `cargo test` filters at once matched **zero** tests and reported
   "0 failed" — a false green that would have scored V6/V7 as harmless. Switched to whole-suite runs.
3. Whole-suite runs without `--no-fail-fast` stopped after the lib target, hiding every integration
   result. V6 read as 2 failures when it is 4.
4. Disk hit 100% (119 MiB free of 926 GiB); the release build failed with `os error 28`. `cargo clean`
   on three spent worktrees reclaimed **13.3 GB**.
5. The saved base binary was **silently truncated** (1.4 MB of 8 MB) because the disk filled during the
   copy, and the original had just been cleaned. Rebuilt from source in `rebase441`, whose source is
   `git diff`-identical to `dev`.
6. V4's `$MODE` unquoted stayed a single argument under zsh, so `--hardtrim5 30` was a clap error
   (exit 2), not the refusal. Redone with a function taking `"$@"` — all four then verified.
7. V15's `grep` patterns beginning with `--` were parsed as options and reported MISSING. Redone with
   `grep -qF --`.

Four of the seven iterations were **harness false negatives, not code defects** — items 2, 3, 5, 6 and
7 each produced a confident wrong answer. Every one was caught by asking whether the check could have
failed, not by reading its output.

---

## §13 — Post-implementation review round (2026-08-17)

Phase 5 ran dual `code-reviewer` on `9a738a3` plus `plan-manager` Mode B.

**`COVERAGE.md` — Verdict: COMPLETE.** 100 items, 95 DONE, 0 PARTIAL, 0 MISSING, 5 DEVIATED (V3 plus
D1–D4, each description checked against the code). It independently measured 764 passing / 17
binaries, +12 tests counted per file, `fmt` clean, and re-verified V1, V4, V5, V15. Two declared gaps:
clippy not run (disk at 3.6 GB), and the mutation rows unverifiable by a read-only audit.

**Reviewer A — READY TO PUSH.** 0 Critical, 0 High, 2 Medium, 5 Low. All four mutation rows it re-ran
reproduced. **Reviewer B — ready with changes,** 1 High, 4 Medium, 5 Low; it reproduced 6 of 11 rows,
including the *honest* V3 that D3 admits was run as an equivalent.

### Fixed in the amended commit (`2142577`)

| Finding | Fix |
|---|---|
| **B's F-1 (High)** — a permuted **top-level** JSON key label survived all 764 tests | Reproduced it (764 green, confirmed), then asserted the four top-level keys *and* their agreement with the block inside the distinct-value test. Re-running the mutation now kills that test. **This was the T-1 lesson applied in one place and not the other**: the agreement loop added by step 13 runs on `pbat` (8/8/8/8), so it could not see a permutation on either side, and the distinct 1/2/3/4 fixture was built and then only half used |
| **B's F-4 (Medium)** — `cli.rs:1071` still claimed the `--fastqc` block was "Last of the fallible checks" | Reworded; our block took that position |
| **B's F-3 (Medium)** — `--library`'s `--help` did not record the refusal, though `--fastqc` and `--rename` do | One clause naming `--clump_only` and the four specialty modes |
| **B's F-2 (Medium)** — a new test comment claimed a general invariant the tool does not hold | Narrowed to what the test actually pins. The assertion stays: it is the C-2 guard |
| **B's F-7 (Low)** — docs listed four of the five modes that refuse | `--clump_only` added |

Re-verified after the fixes: **764 passing, 0 failed**, `fmt` clean, `clippy -D warnings` clean.

### Open — for Felix, not blockers

- **A's F-1 (Medium)** — on a single-end run, `--library emseq --clip_R2 5` prints a gated
  `Clipping in force:` line and then an ungated override line naming `--clip_R2 5`, which reads as
  though it took effect. Both sentences are true; the juxtaposition misleads, and it is the same class
  this commit fixes. Strictly better than the parent, which printed a wrong-but-consistent line. A
  found the sharpest argument: `cli.rs:1056` already warns for `-a2` on a single-end run while
  `--clip_R2` gets nothing at all. **Recommended follow-up: give the clip-R2 family the `-a2`
  treatment** — one `validate()`-time warning — which fixes the preset and non-preset cases together.
  A explicitly said not to hold the push for it. B's F-6 is the same ground.
- **A's F-2 (Medium)** — only one of the four precedence relations the refusal depends on is pinned by
  a test; its *position* is not. Latent, cheap to close with three more ordering tests.
- **B's F-5 (Medium, scope decision)** — the refusal covers `--library` but not the four primitive clip
  flags it expands into, which remain silent no-ops under the specialty modes (measured). Consistent
  with the plan's scope, and it is the pre-existing behaviour the `--phred64` comment describes; worth
  a decision rather than silence.
