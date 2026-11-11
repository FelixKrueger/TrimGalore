# PLAN v2 — #440 follow-ups: give the clip flags the full `-a2` treatment

**Feature dir:** `plans/08172026_clip-r2-se-warning/`
**Upstream:** `plans/08162026_pr441-library-presets/PLAN.md` **§13** — Reviewer A's F-1 and F-2, plus
Reviewer B's F-5, whose user-visible half this now absorbs.
**Base:** `dev` @ `fdade7e`. Branch `fix/440-clip-r2-se-warning`, worktree `/private/tmp/claude-501/fu440`.
**Status:** **implemented** as **`ddc8230`** (amended from `3161b8a` after the review round). All 18
validation rows verified in §12; the review round and its fixes are in §13. Unpushed.

## Revision history

**v2 (2026-08-17)** — v1 was reviewed by two independent `plan-reviewer` agents; both returned READY
WITH CHANGES, A with one Critical and B with two. Felix chose the **full `-a2` treatment**.

| Change | Source |
|---|---|
| **The R2 override lines are now gated too**, on both surfaces (stderr banner and text report), via a new `read_2: bool` on `ClipOverride` | A's C-1 + B's C-2. v1 shipped exactly the state `-a2` calls a bug |
| **The reason is now mode-accurate**, via an `unusable_reason` ladder | B's C-1 / A's I-2. v1's bare `!self.paired` told a two-file `--clock` run "this is a single-end run" |
| **All four clip flags** are named under the specialty modes, not just the two R2 ones | falls out of the ladder; this is F-5's user-visible half |
| One warning **per flag**, not one combined line | A's I-1 / B's I-1 — v1's line read "`--clip_R2 5` and `--three_prime_clip_R2 7` **was** given" |
| Tests move to a **new file** with a mode table | A's Q-1 / B's O-3 — the precedent lives in `integration_adapter2.rs`, not among preset tests |
| Baseline corrected **764 → 768** | A's I-3 / B's I-3, measured independently by both. #446 landed after v1's figure |
| Docs step redirected; its stated precedent was false | A's I-5 / B's I-2 — no `--clip_R1/R2` line exists in `flags.md`, and `-a2`'s warning is documented nowhere |
| CHANGELOG **folds into #440's entry** instead of adding a second bullet | B's O-4 — #440's entry is itself still unreleased, so two bullets citing #440 would ship |
| Anchor drift fixed; assumption 4 narrowed | A's O / B's O-1, I-5 |
| New V-row pinning the stderr/text/JSON triple | B's I-4, which is a better row than v1's near-tautological V10 |

---

## Goal

Make the clip flags behave like `-a2` already does: **never claim a clip that a run cannot perform, and
never say something false about why.** Three deliverables:

1. **Warn** when a clip flag the user typed cannot be honoured — with a reason true for that mode.
2. **Gate the override lines** so no surface announces a Read 2 clip on a single-end run.
3. **Pin the refusal's precedence** with the three missing ordering tests (F-2).

## Context

### Why v1's design was wrong, in the words of the convention it cited

v1 justified itself as completing the `-a2` convention, citing `cli.rs:1048-1071`. It never read the
test that records *why* that implementation looks the way it does:

```rust
// tests/integration_adapter2.rs:414-422
/// `--small_rna` is the load-bearing case: … An earlier version of this
/// fix gated the `Adapter 2 (Read 2):` line but not the NOTE, so a single-end run
/// warned that `-a2` was unused and then announced that `-a2` had displaced the
/// preset default.
```

Substitute the clip vocabulary and that is v1 exactly: `2142577` gated the `Clipping in force:` line,
v1 added a warning, and the override line still announced that `--clip_R2 5` had displaced the preset
value. `-a2` treats that as the bug and fixes it by gating **both**, pinned with `!stderr.contains(…)`.

The same file pins reason accuracy, table-driven over five modes, with the failure message
*"warning must give a reason true for this mode"* (`:405-410`) — and its `clock`/`implicon` rows pass
**two files and no `--paired`**, which is the shape that made v1's reason false.

### The two inert-flag sets are different, and that shapes the code

| Mode | Which clip flags are inert | Why |
|---|---|---|
| `--hardtrim5/3`, `--clock`, `--implicon` | **all four** | `specialty.rs` reads no clip flag; these arms return before `setup_trimming` |
| single-end (trimming path) | the **two R2** flags only | R1 clipping applies normally |

So the reason ladder decides *both* the sentence and which flags are candidates. A single flat
"is it R2?" test cannot express this, which is why v1 collapsed the two cases and got the reason wrong.

### Verified facts this plan rests on (both reviewers, independently)

- `:1138` is after the last `bail!`/`return Err` in `cli.rs` (line **1132**).
- A preset always populates all four `Option`s (`library.rs:167-193`, `pick` returns `Some` on both
  arms), so the warning must read the **raw** `self.clip_r*` fields. Both reviewers ran the mutation:
  reading `effective_clips()` makes a bare `--library emseq` single-end run warn about flags the user
  never typed.
- `--rrbs`'s auto-`clip_r2 2` is paired-gated and writes a **local** binding, never `cli.clip_r2`
  (`main.rs:1141`) — it cannot interact.
- `--clump_only` bails on all three flags long before `:1138` (`cli.rs:839`, `:848`, `:852`).
- All three F-2 relations hold, verified by control flow **and** live binary. `--clock --paired
  --library pbat R1` does *not* slip `validate_paired_input`'s early `Ok`, because that requires
  `mode_label == "Paired-end"` and `:1029` passes `"--clock"`.
- `check_restartable_input`'s message is `"Input file not found"` (`cli.rs:519-525`).
- Baseline on `fdade7e`, macOS, debug: **768 passed, 0 failed, 17 binaries**.

## Behavior

### B1 — the warning

1. After every `bail!` in `validate()`, compute a reason via a ladder in the `-a2` order — specialty
   modes **before** `!paired`, so the more specific truth wins.
2. The reason determines the candidate flag set: all four under the specialty modes, the two R2 flags
   under single-end.
3. Emit **one warning per candidate flag the user actually set**, then continue. Never an error.
4. Paired trimming runs, and runs where the user set no candidate flag, emit nothing.

Reason wording reuses the vocabulary already shipped in the `--library` refusal — *does not honour the
clip flags* — because these modes **do** clip, at fixed positions of their own (`specialty.rs:107,219`
call `clip_5prime`; `:241`,`:352` document it). "Performs no clipping" would be false.

Edge cases:
- `--clip_R2 0` warns. An explicit zero is still an instruction the mode cannot honour.
- Multi-file single-end runs warn once per flag for the invocation, not per file.
- `--paired` with one interleaved uBAM: `self.paired` is true and R2 clipping genuinely applies after
  de-interleaving, so silent. Verified.
- A preset's values are invisible here by construction (raw fields).

### B2 — gating the override lines

1. `ClipOverride` gains `read_2: bool`.
2. Both loops that print override lines skip an override where `read_2 && !paired`:
   the stderr banner (`main.rs:1128`) and the text report (`report.rs:284`).
3. Paired output is unchanged. Single-end output loses only override lines for R2 flags.
4. **This changes report bytes** for single-end runs that override an R2 flag — the one behaviour v1
   claimed it would not touch.

Result, the invariant to hold: *no surface — banner, text report — announces a Read 2 clip on a
single-end run.* The `Clipping in force:` line already obeys it (`2142577`); this extends it to the
override lines.

### B3 — the ordering tests

Three `Cli::validate()` unit tests asserting the *other* error reports first and `--library` is not
mentioned.

## Implementation outline

### 1. `src/library.rs` — `ClipOverride` learns which read it belongs to

```rust
/// One clip flag given on the command line over a preset that also sets it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClipOverride {
    /// Flag name as the user would type it, e.g. `--clip_R1`.
    pub flag: &'static str,
    pub preset_value: usize,
    pub user_value: usize,
    /// A Read 2 flag, which a single-end run cannot honour.
    pub read_2: bool,
}
```

`pick` gains a `read_2: bool` parameter, passed at its four call sites as `false, true, false, true` —
the same order and the same values as `flag_summary`'s tuple array.

**The hazard this creates, and the guard for it.** There are now two lists encoding which flag is a
Read 2 flag: `flag_summary`'s array and `pick`'s four calls. They can drift. Deriving `read_2` from
`flag.ends_with("_R2")` in production is rejected — string-matching a flag name is what Reviewer B's
T-2 mutation caught last round. Instead, step 6 adds a test using the flag *name* as an independent
oracle, so a mis-set bool fails a test rather than shipping.

### 2. `src/main.rs:1128` and `src/report.rs:284` — gate both override loops

```rust
for o in clips.overrides.iter().filter(|o| cli.paired || !o.read_2) {
```

and in `report.rs`, the same with `config.paired`. Nothing else in either loop changes.

### 3. `src/cli.rs:1138` — the warning, with the `-a2` ladder

Insert before the `// Deprecation warnings for Perl-era flags` block:

```rust
// A mode that cannot honour a clip flag says so, rather than ignoring it silently.
let inert = if self.hardtrim5.is_some() || self.hardtrim3.is_some() {
    Some(("--hardtrim5/--hardtrim3 do not honour the clip flags", true))
} else if self.clock {
    Some(("--clock does not honour the clip flags", true))
} else if self.implicon.is_some() {
    Some(("--implicon does not honour the clip flags", true))
} else if !self.paired {
    Some((
        "it applies to Read 2 of a pair, and this is a single-end run",
        false,
    ))
} else {
    None
};
if let Some((reason, all_four)) = inert {
    for (flag, value) in [
        ("--clip_R1", if all_four { self.clip_r1 } else { None }),
        ("--three_prime_clip_R1", if all_four { self.three_prime_clip_r1 } else { None }),
        ("--clip_R2", self.clip_r2),
        ("--three_prime_clip_R2", self.three_prime_clip_r2),
    ] {
        if let Some(n) = value {
            eprintln!(
                "WARNING: {flag} {n} was given but is not used in this mode ({reason}). Ignoring."
            );
        }
    }
}
```

One line of comment, stating the fact. The sentence is `-a2`'s verbatim, per flag, so number agreement
is automatic.

### 4. `src/cli.rs` tests — the ordering trio (F-2)

Beside `hardtrim5_range_error_precedes_the_library_refusal` (`:2480`). Each asserts the expected
substring **and** `!err.contains("--library")`.

1. `hardtrim_conflict_precedes_the_library_refusal` — `--hardtrim5 20 --hardtrim3 20 --library emseq R1`
   → `"cannot be combined in one invocation"`.
2. `clock_arity_error_precedes_the_library_refusal` — `--clock --paired --library pbat R1` (one file)
   → `"requires an even number of input files"`.
3. `missing_input_precedes_the_library_refusal` — `--hardtrim5 20 --library emseq no_such_file.fq.gz`
   → `"Input file not found"`.

### 5. `src/cli.rs` tests — the warning is a warning

4. `inert_clip_flag_is_a_warning_not_an_error` — `--clip_R2 5 R1` and `--hardtrim5 20 --clip_R1 3 R1`
   both `validate()` to `Ok`. (`eprintln!` is invisible to the test harness, so the message itself is
   tested at the binary level in step 7.)

### 6. `src/library.rs` tests — the two-list drift guard

5. `every_override_knows_which_read_it_belongs_to` — resolve `accel` with all four flags overridden,
   then for each produced `ClipOverride` assert `o.read_2 == o.flag.ends_with("_R2")`. The flag name is
   the oracle; production code never inspects it. Kills a mis-set bool at either list.

### 7. `tests/integration_inert_clip_flags.rs` (new) — the mode table

Modelled on `tests/integration_adapter2.rs:386-412`, which is where this exact shape already lives.
Named for what it covers, since it is no longer R2-only.

6. `inert_clip_flag_warns_with_a_mode_accurate_reason` — table over `hardtrim5`, `hardtrim3`, `clock`,
   `implicon`, `single_end`. Each row asserts the run succeeds, stderr contains
   `"is not used in this mode"` **and** that row's reason substring. Specialty rows pass **two files and
   no `--paired`** (the shape that broke v1) and use `--clip_R1`, proving the all-four set. Use
   `clock_10K_R{1,2}.fastq.gz` for the clock/implicon rows, as `-a2`'s table does.
7. `paired_clip_flags_do_not_warn` — `--paired --clip_R2 5 --three_prime_clip_R2 7 R1 R2`: no warning.
8. `preset_alone_on_single_end_does_not_warn` — `--library emseq R1`: stderr has **no**
   `"is not used in this mode"`. The V4 regression guard.
9. `each_inert_flag_gets_its_own_line` — `--clip_R2 5 --three_prime_clip_R2 7 R1`: **two** warning
   lines, each naming its own flag. Covers the grammar fix and step 3's per-flag loop.

Every invocation passes `-o <tempdir>`. **Never omit `-o`** — trimming writes beside the input and
specialty modes write to the CWD; a v1 step that omitted it would have written into tracked
`test_files/`.

### 8. `tests/integration_library.rs` — the gating, end to end

10. `single_end_announces_no_read_2_clip` — `--library emseq --clip_R2 5 R1`: **neither** stderr **nor**
    the text report contains `--clip_R2`. This is the assertion that closes F-1; it fails today.
11. `paired_still_announces_r2_overrides` — `--paired --library emseq --clip_R2 5 R1 R2`: the override
    line **is** present in both. Guards an over-eager gate.

### 9. `src/report.rs` tests — the report half in isolation

12. `test_text_report_omits_r2_override_lines_on_single_end` — `paired = false`, `emseq` with
    `clip_r2: Some(5)`: the report has no `--clip_R2` line. The unit-level twin of step 10.

### 10. Docs

13. `docs/src/content/docs/guide/flags.md` — a **new bullet** in the cross-flag-interactions list,
    shaped like the `--fastqc` one. v1's claim that it would match how `-a2` is described was false:
    `grep -rn "not used in this mode" docs/` returns zero, so `-a2`'s warning is documented nowhere.
    Documenting the clip case and not `-a2`'s is a deliberate asymmetry, recorded in Q-1.

### 11. `CHANGELOG.md` — fold one clause into #440's existing entry

**No second bullet.** #440's entry is itself still under `### Unreleased`, so a new entry would ship two
bullets citing the same issue. Add one clause to the paragraph that already describes the refusal:

> … For the same reason `--library` is refused with `--hardtrim5`, `--hardtrim3`, `--clock` and
> `--implicon`, which do not honour the clip flags — and a clip flag those modes cannot honour, or a
> Read 2 clip flag on a single-end run, now warns instead of passing silently.

## Efficiency

Not on any hot path. The ladder is at most four `Option::is_some()` checks and one four-element array
walk, once per invocation. The two gated loops gain one `bool` test per override, at most four per run.

## Integration

**Reads:** `Cli::paired`, `clip_r1`, `clip_r2`, `three_prime_clip_r1`, `three_prime_clip_r2`,
`hardtrim5`, `hardtrim3`, `clock`, `implicon`; `ClipOverride::read_2`.
**Writes:** stderr warnings; the stderr banner's override lines; **the text report's override lines**.

**The text report changes** for single-end runs that override an R2 flag — stated plainly, because v1
listed "touches no report" as a virtue when it was the boundary of the defect.

**JSON is unchanged**, and that is now a recorded decision rather than an omission: the JSON records
what was *requested* (`clip_r2: 5` on a single-end run), the text report what was *in force*. A
consumer distinguishes them with the top-level `mode` key. V8 pins the resulting stderr/text/JSON
triple so the asymmetry is asserted, not assumed.

**Downstream:** `--rrbs` unaffected (paired-gated, local binding). `validation` CI unaffected — it
passes no clip flag on any shape, and warnings are stderr-only. No existing test breaks: the only
single-end-plus-explicit-R2 invocation in the suite is an oracle run asserting bytes, which a stderr
line cannot affect.

## Assumptions

1. The warning reads the user's **raw** `self.clip_r*` fields, never `effective_clips()`.
2. `--clip_R2 0` warrants a warning.
3. One line per flag, matching `-a2`'s per-flag construction.
4. Placing the warning after every `bail!` keeps a **`validate()`-level** refusal quiet. Later refusals
   — `sanity_check_any`, uBAM format validation, the output-collision pre-flight — all run after
   `cli.validate()?` at `main.rs:242` and will still be preceded by the warning. Same property `-a2`
   already has.
5. The specialty modes cannot honour **any** of the four clip flags (`specialty.rs` reads none).
6. `read_2` is set at `pick`'s call sites, not derived from the flag string; step 6's test is what keeps
   the two lists in agreement.

## Validation

Baseline **768**. Expect roughly **780**; measure it, do not assume. Use `--no-fail-fast` — without it
cargo stops after the lib target and hides all 283 integration results. Never gate on one number
(macOS and Linux differ by two). One `cargo test` filter per invocation: multiple filters match nothing
and report a false "0 failed".

| # | Verify | How | Expected |
|---|---|---|---|
| V1 | The reason is true in every mode | step 6's table | 5 rows pass |
| V2 | It can fail | delete the ladder's three specialty arms, leaving bare `!self.paired` | the clock, implicon, hardtrim5 and hardtrim3 rows go red |
| V3 | All four flags are named under the specialty modes | step 6 rows use `--clip_R1` | pass |
| V4 | It can fail | hard-code `all_four = false` | the specialty rows go red |
| V5 | No noise on a preset-only single-end run | step 8 | pass |
| V6 | It can fail | read `self.effective_clips().clip_r2` instead of the raw field | step 8 red — both reviewers reproduced this |
| V7 | No surface announces an R2 clip on single-end | steps 10 and 12 | pass |
| V8 | It can fail | revert either gate — `main.rs` only, then `report.rs` only | step 10 red on each, separately |
| V9 | The paired half survives | step 11 | pass |
| V10 | It can fail | change the filter to `|o| !o.read_2` (drop R2 always) | step 11 red |
| V11 | `read_2` matches the flag name | step 5 | pass |
| V12 | It can fail | swap two `read_2` arguments at `pick`'s call sites | step 5 red |
| V13 | One line per flag | step 9 | pass, exactly 2 lines |
| V14 | The three precedence relations hold | steps 1-3 | pass |
| V15 | Each can fail | move the `--library` refusal above `cli.rs:1011` | all three red, plus the pre-existing `hardtrim5_range_error_…` |
| V16 | The stderr/text/JSON triple is what we intend | single-end `--clip_R2 5`: stderr warns, text names no R2 clip, JSON still records `clip_r2: 5`, `mode` is `single-end` | all four hold |
| V17 | Suite and lints | `cargo test --no-fail-fast`; `cargo clippy --all-targets --release -- -D warnings`; `cargo fmt --all -- --check` | green. **Run clippy in release** — both reviewers ran debug only, for disk reasons, and flagged release as still owed |
| V18 | Paired output is byte-identical | base-vs-new differential on `--paired --library emseq --clip_R2 5`, `-o` a fresh dir, counting files | every FASTQ, `.txt` and `.json` identical |

## Questions or ambiguities

- **Q-1 (Open, decided).** Documenting the clip warning while `-a2`'s stays undocumented is an
  asymmetry. Assumption: document the clip case anyway; a missing doc for `-a2` is not a reason to
  leave two features undocumented. Filing an issue for `-a2`'s docs is left to Felix's filing bar.
- **Q-2 (Open).** The JSON still records an inapplicable `clip_r2` on single-end runs (§Integration).
  Assumption: leave it; changing it is a schema change for every MultiQC consumer. V16 pins the
  current triple so a future change is a deliberate act.
- **Critical:** none.

## Self-Review

**What v2 fixed, and the pattern behind it:**

- **v1 cited a precedent without reading how it was resolved.** The `-a2` implementation was cited from
  `cli.rs`; the reasoning lived in a test doc comment in `integration_adapter2.rs`, and it explicitly
  rejects v1's design. In this repo, hard-won reasoning is in the tests — a citation that stops at the
  implementation is half a citation.
- **v1 collapsed two different inert-flag sets into one predicate.** Once the reason ladder exists, the
  "which flags" question answers itself per branch. v1's single `!self.paired` test is what made both
  the reason false and the flag set wrong.
- **v1 described its scope limit as a virtue.** "Touches no report, no JSON, no output file" read as a
  safety property; it was the boundary of the defect. §Integration now states the report change plainly.
- **v1's docs step rested on something that does not exist**, and its baseline was stale. Both were
  caught by two reviewers independently, which is the argument for the dual convention.
- **The new `read_2` field creates a two-list drift risk**, which v2 answers with a test using the flag
  name as an oracle rather than by string-matching in production.

**Remaining risks:**

- **R-1** — New output on previously silent invocations, now on more shapes than v1 (all four flags
  under the specialty modes). That is the point, and it matches `-a2`; nothing parses stderr in CI.
- **R-2** — The text report changes for single-end R2-override runs. No golden fixture covers it (the
  `validation` job passes no clip flags), but V18 pins the paired path as unchanged.
- **R-3** — F-5's remaining half: under the specialty modes the flags now *warn*, but they are still
  accepted rather than refused, unlike `--library`. Deliberate — refusing them would remove working
  invocations for flags that have shipped for years, which is a different decision from refusing a
  brand-new flag.

---

## §12 — Implementation notes (2026-08-17)

**Commit `3161b8a`** on `fix/440-clip-r2-se-warning`, base `dev` @ `fdade7e`. 8 files, **+437 / −18**.
**780 tests pass** (768 base + **12 new**), 0 failures, 18 binaries. `cargo fmt --all -- --check` clean.
**`cargo clippy --all-targets --release -D warnings` clean** — run in release, which both plan reviewers
had left owed for disk reasons.

All **13 numbered outline items** across §1–§11 landed, and all **18 validation rows** are verified.

### Mutation matrix — every new test proved able to fail

| V | Mutation | RED |
|---|---|---|
| V2 | ladder loses its three specialty arms (bare `!paired`) | `inert_clip_flag_warns_with_a_mode_accurate_reason` |
| V4 | `all_four` hard-coded `false` | same table |
| V6 | warning reads `effective_clips()` instead of the raw field | **2**: `preset_alone_on_single_end_does_not_warn` + `explicit_clip_flag_overrides_the_preset_and_is_logged` |
| **V8a** | **only** the `main.rs` gate reverted | `single_end_announces_no_read_2_clip` |
| **V8b** | **only** the `report.rs` gate reverted | **2**: the report unit test + the integration test |
| V10 | gate drops R2 unconditionally | `paired_still_announces_r2_overrides` |
| V12 | two `read_2` arguments swapped at `pick`'s call sites | **6**, incl. the oracle test |
| V15 | the three earlier checks suppressed | **9**: all three new ordering tests + 6 pre-existing that share those checks |

V8a and V8b matter most: each surface is covered **independently**, so neither gate can be reverted
silently. All four mutated files `cmp`-verified byte-identical after every run.

### V16 — the stderr/text/JSON triple, measured

Single-end `--clip_R2 5`: stderr warns and names the flag; the text report contains **zero** mentions of
`--clip_R2`; the JSON records `clip_r2: 5` with `mode: single-end`. The asymmetry is now asserted rather
than assumed.

### V18 — blast radius, measured against a base-commit build

Five shapes, 24 files compared, none producing zero output. Every FASTQ and every `.txt` identical on all
paired shapes and on single-end without an R2 override. All eight JSON files identical once
`command_line` — which records the two binary paths — is excluded. **The only behavioural output change
in the entire differential is one line:**

```
14d13
< --clip_R2 5 was given on the command line and overrides the emseq preset value 10
```

Exactly the intended change, on exactly the intended shape.

### Deviations

- **D1** — `cargo build` first failed with `error: cached cgu … should have an object file, but doesn't`,
  a corrupted incremental cache left by the earlier disk exhaustion, not a code fault. `cargo clean`
  fixed it. Worth recognising: it reads like a compiler error and is a disk artefact.
- **D2** — rustfmt reflowed one test's argument array (`clock_arity_error_precedes_the_library_refusal`).
  Applied `cargo fmt`; no behaviour change.
- **D3** — step 10's assertion had to be narrowed while writing it. The plan said the single-end test
  asserts stderr contains no `--clip_R2`, but **the warning names the flag on purpose** — it says the
  value was *not* used. Asserting the flag's absence would have contradicted the feature. The test
  targets the *override sentence* (`"--clip_R2 5 was given on the command line"`) and separately asserts
  the warning **is** present. Caught before running, by reading the assertion against the feature.
- **D4** — V15 was executed as an equivalent (suppress the three earlier checks) rather than by
  relocating the refusal block. Same property proved: each earlier check must report first.

---

## §13 — Post-implementation review round (2026-08-17)

Phase 5 ran dual `code-reviewer` on `3161b8a` plus `plan-manager` Mode B. **Commit amended to
`ddc8230`** — 8 files, +495/−22, **13 new tests, 781 passing**, `fmt` clean, `clippy --release
-D warnings` clean.

**`COVERAGE.md` — Verdict: COMPLETE.** 65 items, 60 DONE, 0 PARTIAL, 0 MISSING, 5 DEVIATED. It
independently measured 780/0/18 at `3161b8a`, confirmed D1–D4 accurate line by line, confirmed
Assumption 6's two-list guard is real (no flag-string inspection in production), and confirmed Q-2 and
R-3 were honoured as deferred rather than silently implemented.

**Both reviewers: READY WITH CHANGES, no defect in shipped behaviour.** A enumerated **41** live
invocations, B **32**; neither found a false reason, an unwarned inert flag, or the two surfaces
disagreeing. Both reproduced V6/V8a/V8b/V12, V16 and V18 independently — A's differential ran **43
files, 39 identical, 4 differing, all four the same single line**, and covered the two uBAM report sites
§12's own differential had missed.

### Fixed in `ddc8230`

| Finding | Fix |
|---|---|
| **A's H-1 (High)** — dropping the `all_four` narrowing shipped **green**; a plain `--clip_R1 3 <se>` run would print a sentence false in both halves | New `read_1_clip_flags_do_not_warn_on_single_end`. Reverting the narrowing now turns it red |
| **M-1 (both)** — `--three_prime_clip_R1` had **no test**; deleting its whole entry from the warning array left all 780 green. B noted this is the class A's plan-review I-6 raised, closed for the R2 3' flag and left open for the R1 one | The mode table now carries **both** flags of each inert set, so all four array entries are exercised. Deleting the entry now turns the table red |
| **M-2 (both)** — `flags.md` said "The trimming report lists only the clipping that applied", false of the **JSON** report | Says "the text trimming report", and names the JSON's `mode` key |
| A's docs items — the CHANGELOG's "either surface" is false if JSON counts | "neither the log nor the text report" |
| Coverage's observation 1 — §12 claimed "23 outline steps"; the outline has **13** numbered items | Corrected |
| Coverage's observation 2 — one CHANGELOG line ran 111 chars against ~95 neighbours | Reflowed; max is now 97. Took three attempts, because patching one line of a wrapped block just moves the overflow — reflow is a paragraph operation |
| Clippy regression introduced by the M-1 fix | The 5-tuple table tripped `clippy::type_complexity`; replaced with a named `Case` struct rather than `#[allow]`-ing it, per CLAUDE.md |

### Open — for Felix

- **B's M-3 (Medium), reproduced.** A single-end `--rename --output-format ubam --clip_R2 5` run prints
  *"`--clip_R2 5` … is not used in this mode … Ignoring."* and is then **refused because `--clip_R2` is
  set** — the #408 guard at `main.rs:398` reads `effective_clips().any_set()`, which counts an inert
  Read 2 flag. So the tool says it is ignoring the flag and then refuses on account of it.

  This is **the third place** encoding "which read does a clip flag belong to", and the commit did not
  teach it — B found it by answering the review question directly. The over-refusal is pre-existing;
  this commit's invariant is what makes it read as a contradiction.

  Not fixed here because relaxing a shipped refusal needs its own decision: the #408 guard exists to
  stop `--rename` silently dropping annotations, and narrowing it to ignore inert flags is a behaviour
  change to a data-safety guard, not a wording fix.
