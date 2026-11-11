# PLAN REVIEW — Reviewer B (generation g1)

**Plan:** `plans/08172026_clip-r2-se-warning/PLAN.md`
**Reviewer:** B (independent; Reviewer A's report not read)
**Worktree:** `/private/tmp/claude-501/prB440` @ `fdade7efb8b579fd9b25899d04b1609dba860b20` (current `dev` tip)

```
$ git rev-parse HEAD
fdade7efb8b579fd9b25899d04b1609dba860b20
$ git worktree list
/Users/fkrueger/Github/TrimGalore  fdade7e [dev]
/private/tmp/claude-501/fu440      fdade7e [fix/440-clip-r2-se-warning]
/private/tmp/claude-501/prA440     fdade7e (detached HEAD)
/private/tmp/claude-501/prB440     fdade7e (detached HEAD)
```

---

## Verdict — **READY WITH CHANGES** (two Critical)

The plan is unusually well-built: I applied its step-1 snippet and its three F-2 tests verbatim and
everything compiled, passed, `fmt` clean, `clippy -D warnings` clean, and produced exactly the stderr
the plan predicts for all four cases it names. The central design decision (raw fields, never
`effective_clips()`) is **correct**, and V4's mutation is **not** vacuous — I reproduced the regression.

Two things must change before implementation:

1. **C-1 — the warning states a false reason on three legal invocations.** `--clock`, `--implicon` and
   `--hardtrim5/3` do not require `--paired`, so a genuinely paired two-file run reaches the
   `!self.paired` branch and is told *"this is a single-end run"*. This repo already asserts the
   opposite invariant in a test whose failure message is literally **"warning must give a reason true
   for this mode"** (`tests/integration_adapter2.rs:408`). Q-2's exclusion is coherent as a *coverage*
   decision but not as a *correctness* one — the plan's snippet needs either the `-a2` reason ladder
   (which is Q-2) or a two-condition suppression guard.
2. **C-2 — the plan does not remove the juxtaposition F-1 is about.** The warning lands 12 lines above
   the override line it corrects, separated by the version banner. The `-a2` family already faced this
   exact situation and its fix documents warning-alone as the insufficient half:
   `tests/integration_adapter2.rs:414-422`.

Everything else is Important-or-below. The baseline figure is stale (768, not 764) and step 12's docs
anchor does not exist.

---

## Claim-verification table

### Anchors (task item 1)

| Claim | Result | Settled by |
|---|---|---|
| `:1013` — #386 hardtrim5+hardtrim3 bail | VERIFIED (message line; `if` is at 1011, `bail!` at 1012) | `src/cli.rs:1011-1017` |
| `:1021` — `--hardtrim5` range bail | VERIFIED exact | `src/cli.rs:1021` |
| `:1026` — `--hardtrim3` range bail | VERIFIED exact | `src/cli.rs:1026` |
| `:1029` — `validate_paired_input("--clock")` | VERIFIED exact | `src/cli.rs:1029` |
| `:1032` — `validate_paired_input("--implicon")` | VERIFIED exact | `src/cli.rs:1032` |
| `:1045` — `check_restartable_input(path, "Input file not found")` | VERIFIED exact | `src/cli.rs:1045` |
| `:1048-1070` — the `-a2` warning block | VERIFIED (block is 1048-**1071**; off by one at the end) | `src/cli.rs:1048-1071` |
| `:1114` — the `--library` refusal | VERIFIED (comment at 1114, `if self.library.is_some()` at 1115) | `src/cli.rs:1114-1136` |
| `:1138` — deprecation-warning block, the insertion point | VERIFIED exact | `src/cli.rs:1138` |
| `:2480` — `hardtrim5_range_error_precedes_the_library_refusal` | VERIFIED exact (`fn` line) | `src/cli.rs:2480` |
| `:2565` — the #421 ordering-test section header | VERIFIED exact | `src/cli.rs:2565` |
| `:837` — the `--clump_only` `--clip_r2` bail | **REFUTED as an exact anchor** — 837 is the closing brace of the `clip_r1` block; the guard is 838, the `bail!` 839. Claim substantively holds | `src/cli.rs:838-839` |

```
$ grep -n "clip_R2\|clip_r2" src/cli.rs | head
src/cli.rs:838:            if self.clip_r2.is_some() {
src/cli.rs:839:                anyhow::bail!("--clump_only does not clip; --clip_r2 is not compatible");
src/cli.rs:846:            if self.three_prime_clip_r2.is_some() {
src/cli.rs:848:                    "--clump_only does not clip; --three_prime_clip_r2 is not compatible"
```

### The substantive claims (task items 2-9)

| # | Claim | Result | Settled by |
|---|---|---|---|
| 2 | "There is no warning today", zero matches in `cli.rs`/`main.rs`; nothing equivalent elsewhere | **VERIFIED** | `grep -rn "is not used in this mode" src/ tests/` → exactly one hit, `src/cli.rs:1065` (the `-a2` warning). The only other `clip_r2` diagnostics anywhere are the `--clump_only` bails (838/846) and the RRBS auto-set note (`main.rs:1143`). Nothing in `setup_trimming` |
| 3 | A preset always populates all four `Option`s, so `effective_clips()` would warn on every SE `--library` run | **VERIFIED**, twice | `src/library.rs:167-193` — `pick` returns `Some(v)` or `Some(preset_value)`, never `None`; and empirically by mutation (see V4 below) |
| 3b | Stress test: does reading the raw fields ever give the wrong answer? | **YES — three cases.** See **C-1**. Not because the fields are raw, but because `!self.paired` is the wrong sole predicate | CASES 5/7/8 below |
| 4 | `:1138` is after every `bail!` in `validate()` | **VERIFIED** — the last `bail!`/`return Err` in the whole of `cli.rs` is line 1132 | `grep -n "bail!\|return Err" src/cli.rs \| tail -1` → `1132` |
| 4b | "a refused run stays quiet" | **VERIFIED but overstated** — true only for `validate()`-level refusals. `cli.validate()?` is `main.rs:242`; `sanity_check_any`, the uBAM format validation and the output-collision pre-flight all bail after it | `src/main.rs:242, 263` |
| 5 | `--rrbs` auto-`clip_r2 2` is in `setup_trimming` and paired-only, so it cannot reach `!self.paired` | **VERIFIED** | `src/main.rs:1141` — `cli.rrbs && !cli.non_directional && cli.paired && clips.clip_r2.is_none()`; the result is a *local* `let clip_r2`, never written back to `cli.clip_r2` |
| 6 | `--clump_only` is unreachable via the bail | **VERIFIED** | `src/cli.rs:839` (`--clip_r2`), `:848` (`--three_prime_clip_r2`), `:852` (`--library`) — all long before 1138 |
| 7 | The three F-2 tests will pass, incl. `clock_arity_error_...` not slipping through the early `Ok` | **VERIFIED empirically.** The early `Ok` requires `mode_label == "Paired-end"` (`cli.rs:557`); `--clock` passes `"--clock"`, so the arity bail fires | ran all four; see output below |
| 7b | Assumption 5 — `check_restartable_input`'s message is `"Input file not found"` | **CONFIRMED** | `src/cli.rs:519-525` formats `"{not_found}: {}"`; `missing_input_precedes_the_library_refusal` passes |
| 8 | Baseline "764 + 8 new" | **REFUTED — stale.** Actual on `fdade7e`, macOS, debug: **768 passed, 0 failed** | see below |
| 9 | Snippet compiles, clippy `-D warnings` clean, correct stderr for the four named cases | **VERIFIED** (all four, plus 7 extra shapes) | see below |
| 9b | Revert verified | **VERIFIED** with `git status --short`, not asserted | see below |

---

## Evidence

### Baseline (task item 8)

```
$ cargo test --no-fail-fast          # 17 targets
485 (lib) + 0 (main) + 15 + 13 + 16 + 11 + 8 + 1 + 8 + 57 + 14 + 21 + 3 + 57 + 19 + 40 + 0 (doc)
SUM passed = 768
SUM failed = 0
```

`768`, not `764` — `#446` landed after the plan's figure was measured. **V9 should read `768 + 8 = 776`
on macOS.** The plan's own caution ("never gate on one test number; macOS and Linux differ by two")
still applies. The plan's `--no-fail-fast` advice is correct: without it cargo stops after the lib
target and you never see the 283 integration results.

### The four new tests (task item 7)

Applied plan steps 4, 9, 10, 11 verbatim:

```
test cli::tests::clip_r2_on_single_end_is_a_warning_not_an_error ... ok
test cli::tests::hardtrim_conflict_precedes_the_library_refusal ... ok
test cli::tests::clock_arity_error_precedes_the_library_refusal ... ok
test cli::tests::missing_input_precedes_the_library_refusal ... ok
test result: ok. 489 passed; 0 failed  (485 → 489)
```

`clock_arity_error_precedes_the_library_refusal` is the one the task flagged as doubtful. It **passes**,
and the reason is the `mode_label` discriminator:

```rust
// src/cli.rs:557
if self.input.len() == 1 && self.paired && mode_label == "Paired-end" {
    return Ok(());
}
```

`--clock --paired --library pbat R1` runs `validate_paired_input("Paired-end")` first (`:661`), which
*does* take the early `Ok`; then `--clock` calls it again with `mode_label == "--clock"` (`:1029`),
which does not, and bails with `--clock mode requires an even number of input files (R1/R2 pairs),
got 1`. The plan's expected substring `"requires an even number of input files"` matches. There is a
near-identical existing precedent at `src/cli.rs:2600` (`odd_input_count_precedes_the_fastqc_refusal`)
whose doc comment records the same subtlety from the other direction.

### Lints and behaviour (task item 9)

```
$ cargo fmt --all -- --check      → exit 0
$ cargo clippy --all-targets -- -D warnings   → Finished; exit 0
```
*(Deviation: I ran clippy in **debug**, not `--release`. The volume dropped to 10 GiB free during the
run and a release build of `resvg`/`fastqc-rust` risked filling it. Same lint set for this code; a
release clippy still owes verification at implementation time.)*

The four cases the task named:

| Case | Result |
|---|---|
| `--clip_R2 5 <se>` | warns, exit 0 |
| `--library emseq --clip_R2 5 <se>` | warns, exit 0 |
| `--library emseq <se>` | **silent** (0 matches for `is not used in this mode`) |
| `--paired --clip_R2 5 R1 R2` | **silent** (0 matches) |

```
$ trim_galore --dont_gzip --clip_R2 5 -o . BS-seq_10K_R1.fastq.gz
WARNING: --clip_R2 5 was given but is not used in this mode (it applies to Read 2 of a pair, and this is a single-end run). Ignoring.

Trim Galore v2.3.0
fdade7e — macos/aarch64 — built 2026-08-17T15:41:23Z
==================================================
...
EXIT=0
```

### V4 — the important mutation (confirmed **not** vacuous)

Mutated `self.clip_r2` → `self.effective_clips().clip_r2`, then ran the case step 7 guards:

```
$ trim_galore --dont_gzip --library emseq -o . BS-seq_10K_R1.fastq.gz
WARNING: --clip_R2 10 was given but is not used in this mode (...). Ignoring.
exit=0 warn=1
```

A flag the user never typed, reported as "was given". Step 7
(`preset_alone_on_single_end_does_not_warn`, asserting stderr does **not** contain
`is not used in this mode`) would go red. **V4 is the plan's best row** and it does what it claims.

### Revert (task item 9b)

```
$ git checkout -- src/cli.rs
$ git status --short
?? BASELINE_TEST.log
?? scratch_b/
$ git rev-parse HEAD
fdade7efb8b579fd9b25899d04b1609dba860b20
$ grep -c "inert" src/cli.rs           → 0
$ grep -c "REVIEWER-B SCRATCH" src/cli.rs → 0
```

No tracked file modified. `cargo clean` run (`Removed 10509 files, 1.8GiB total`). **Left behind:**
`scratch_b/` (8.8M), `scratch_c/` (1.8M), `BASELINE_TEST.log` (132K), all untracked inside this
throwaway worktree — the `rm -rf` was refused by the permission gate and I did not retry it. They go
away with `git worktree remove`.

---

## Critical

### C-1 — the warning states a reason that is false in three modes, breaking an invariant this repo already asserts

`--clock`, `--implicon` and `--hardtrim5/3` are **inherently paired but do not require `--paired`**.
Nothing in `validate()` sets or demands `self.paired` for them — `--clock` only calls
`validate_paired_input("--clock")` (`cli.rs:1028-1033`), which requires an *even file count*, not the
flag. So a two-file clock run reaches the `!self.paired` branch:

```
$ trim_galore --clock --clip_R2 5 R1 R2          # two files, no --paired
exit=0
WARNING: --clip_R2 5 was given but is not used in this mode (it applies to Read 2 of a
pair, and this is a single-end run). Ignoring.
```

It is not a single-end run. Same for `--implicon --clip_R2 5 R1 R2` (warn=1) and
`--hardtrim5 20 --clip_R2 5 R1 R2` (warn=1). The *conclusion* ("not used … Ignoring") is true — I
confirmed `src/specialty.rs` reads no clip flag and the specialty arms return before any
`setup_trimming` call (`main.rs:880/988/1038/2059/2117/2168`) — but the parenthetical reason is false,
and it misdirects: a user told "this is a single-end run" will add `--paired`, which **silences the
warning without making the flag do anything**:

```
$ trim_galore --clock --paired --clip_R2 5 R1 R2
exit=0 warnings=0
```

So under the plan as written, the warning's presence is *anti-correlated* with whether the user has a
problem.

This is not a new standard I am inventing. The repo already tests this invariant for the warning family
the plan is deliberately copying:

```rust
// tests/integration_adapter2.rs:405-410 — table covers clock, implicon, hardtrim5, hardtrim3, single_end
let (ok, stderr) = run(&args);
assert!(ok, "[{tag}] should warn, not fail: {stderr}");
assert!(
    stderr.contains("not used in this mode") && stderr.contains(reason),
    "[{tag}] warning must give a reason true for this mode; got: {stderr}"
);
```

Note that the table's `clock`/`implicon` rows pass **two files and no `--paired`** (lines 392-394) —
exactly the shape above. The `-a2` block's `unusable_reason` ladder (`cli.rs:1052-1062`) exists *for
this reason*, testing the specialty modes **before** `!self.paired`. And §13 of the upstream plan
records Reviewer A's recommendation as "give the clip-R2 family **the `-a2` treatment**" — the ladder
*is* the treatment; the plan implements a simplified version that drops the discrimination.

**On Q-2's scope.** Excluding Q-2 as a *coverage* decision is coherent — leaving `--clip_R1`/
`--three_prime_clip_R1` silent under the specialty modes keeps F-5 exactly as recorded, and R-2 says so
honestly. Excluding it as a *message-correctness* decision is not, because the plan's own snippet emits
a false sentence in those modes. Two ways out:

- **Option A (strictly in scope, ~2 conditions).** Suppress the warning under the specialty modes, so
  the message is never false and F-5 stays open untouched:
  ```rust
  if !self.paired && !self.clock && self.implicon.is_none()
      && self.hardtrim5.is_none() && self.hardtrim3.is_none() {
  ```
  Adds no coverage. Costs: a clock/hardtrim user with `--clip_R2` still gets nothing (status quo).
- **Option B (= Q-2, ~4 more lines).** Reuse the `-a2` ladder shape, so every mode gets a true reason
  and F-5's user-visible half closes.

I'd take **B** — it is the smaller lie-free change *and* it makes the two warnings genuinely one
convention, which is the plan's own stated goal ("the R2 family's missing half of an established
convention"). But A is enough to clear this finding, and A is the one that needs no scope decision from
Felix.

### C-2 — the warning does not remove the juxtaposition F-1 is about, and `-a2` already rejected warning-alone as the fix

The plan's §Context defines the defect as the two lines misleading *together*. Here is what the user
actually gets, F-1's exact invocation:

```
$ trim_galore --dont_gzip --library emseq --clip_R2 5 -o . BS-seq_10K_R1.fastq.gz
WARNING: --clip_R2 5 was given but is not used in this mode (...). Ignoring.   ← line 1

Trim Galore v2.3.0
fdade7e — macos/aarch64 — built 2026-08-17T15:41:23Z
==================================================

Auto-detecting adapter type...
No adapters detected in the first 10000 reads. Defaulting to Illumina adapter.
Adapter: Illumina (AGATCGGAAGAGC)

Library preset 'emseq' selected: --clip_R1 10 --three_prime_clip_R1 10        ← line 12
--clip_R2 5 was given on the command line and overrides the emseq preset value 10   ← line 13
```

`cli.validate()?` is `main.rs:242`, above the banner; the preset lines are printed in `setup_trimming`
(`main.rs:1122-1136`). The corrective sentence is **twelve lines and one banner** away from the two
lines it corrects, and it is the first thing on stderr — before the tool has even said its own name.
The misleading pair is unchanged.

The `-a2` family hit precisely this and its fix says so in a doc comment:

```rust
// tests/integration_adapter2.rs:414-422
/// Single-end must not announce a Read 2 adapter it is about to ignore, and must
/// not contradict its own warning. …
/// `--small_rna` is the load-bearing case: … An earlier version of this
/// fix gated the `Adapter 2 (Read 2):` line but not the NOTE, so a single-end run
/// warned that `-a2` was unused and then announced that `-a2` had displaced the
/// preset default.
#[test]
fn single_end_does_not_announce_a_read2_adapter() { … }
```

"Warned that it was unused and then announced that it had displaced the preset default" **is** the
state this plan ships. `-a2` fixed it by gating *both* the announcement and the NOTE, and the assertion
`!stderr.contains("NOTE: Read 2 adapter taken from -a2")` pins it. The plan's §Context justifies
keeping the override line because it "is a true statement about user input" — but `-a2`'s NOTE was
equally true about user input and was gated anyway.

**The completing fix is ~3 lines in `main.rs:1128-1136`**: skip an override whose flag names Read 2 when
`!cli.paired`, exactly as `ResolvedClips::flag_summary(paired)` already skips them for the banner
(`library.rs:122-130`). That gives you one consistent rule — *nothing on stderr claims a Read 2 clip on
a single-end run* — instead of two lines that disagree.

I am flagging this Critical because the plan's stated goal is not met, not because the approved item is
undelivered: the `validate()`-time warning Felix approved works. If Felix wants only the approved item,
say so in the plan and record the residue as a risk beside R-1; do not leave §Context claiming the
juxtaposition is fixed.

---

## Important

### I-1 — the message does not agree in number, and the plan's assumption 3 is what causes it

```
$ trim_galore --dont_gzip --clip_R2 5 --three_prime_clip_R2 7 -o . <se>
WARNING: --clip_R2 5 and --three_prime_clip_R2 7 was given but is not used in this mode
(it applies to Read 2 of a pair, and this is a single-end run). Ignoring.
```

"**was** given but **is** not used … **it** applies" with a plural subject. The `-a2` sentence the plan
reuses was written for one flag; assumption 3 (one combined line) makes it plural without adapting the
verb. Either pluralise conditionally, or emit one line per flag — which is also closer to `-a2`, whose
warning is per-flag by construction. Not a blocker, but it is a brand-new user-facing string and this
is exactly the class of wording a maintainer bounces.

### I-2 — step 12's docs anchor does not exist, and its stated precedent is false

The plan says: *"`docs/src/content/docs/guide/flags.md` — one clause on the `--clip_R1/R2` line …
matching how `-a2` is already described."* Both halves are wrong.

```
$ grep -n -- "clip_R1\|clip_R2" docs/src/content/docs/guide/flags.md
13:- **`--rrbs`** in paired-end directional mode auto-sets `--clip_R2 2` …
14:- **`--library <preset>`** expands into the four clip flags …
```

There is **no `--clip_R1/R2` line** — only the `--rrbs` bullet and the `--library` bullet. And:

```
$ grep -rn "not used in this mode\|Ignoring" docs/src/content/docs/
(no matches)
```

`-a2`'s warn-and-ignore behaviour is **documented nowhere** in `docs/`, so there is no description to
match. Pick one deliberately: add a new bullet to `flags.md` (and accept that the clip warning is then
documented while `-a2`'s is not), or drop step 12 and stay consistent with `-a2`. Either is fine; the
plan's rationale for step 12 currently rests on something that does not exist.

### I-3 — V9's baseline is stale

768, not 764 (measured above). Expect **776** on macOS after the 8 new tests.

### I-4 — the JSON report contradicts the new warning in the same run

Pre-existing and out of scope, but the warning turns a silent inconsistency into an explicit one:

```
$ trim_galore --dont_gzip --clip_R2 5 --three_prime_clip_R2 7 -o . <se>
stderr:  WARNING: … is not used in this mode … Ignoring.
JSON  :  clip_r2 = 5      three_prime_clip_r2 = 7      library = None
text  :  (no "Clipping in force" line at all)
```

`2142577` gated the *text* report; the JSON still records the raw request. A MultiQC-side consumer
reading `parameters.clip_r2` on a single-end run still sees `5` on the very run whose log says
"Ignoring". Defensible ("JSON records what was requested, text records what was in force") but it
should be a recorded decision, and it makes the plan's §Integration line "Touches no report, no JSON"
technically true and rhetorically misleading. **This would also make a much better V10 than the one the
plan has** (see O-2).

### I-5 — assumption 4 is right but claims more than it delivers

"Placing it after all bails is right, so a refused run stays quiet" holds for `validate()` refusals
only. `cli.validate()?` is `main.rs:242`; `sanity_check_any` (`:263`), the uBAM format validation and
the output-collision pre-flight all bail afterwards, so e.g. a colorspace or truncated input still gets
the warning before its refusal. Same property as `-a2`, so nothing to fix — just narrow the wording so
the divergence from `-a2` is stated accurately.

---

## Optional

- **O-1 — anchor drift.** `:837` is the closing brace, not the `--clip_r2` guard (838/839); `:1013` is
  the message line, not the `if` (1011). Trivial, but the plan is otherwise anchor-exact, so worth
  correcting for the implementer.
- **O-2 — V10 is near-tautological.** An `eprintln!`-only change cannot alter output bytes; the
  differential can fail only if the implementer accidentally edits a report path. Keep it if it is
  cheap, but the row that would actually earn its place is I-4's: assert what the JSON says about
  `clip_r2` on a single-end run, so the text/JSON/stderr triple is pinned rather than assumed.
- **O-3 — Q-1: put the four binary tests in a new file.** Recommendation: **new file**,
  `tests/integration_clip_r2_warning.rs`, modelled on `tests/integration_adapter2.rs` (which is where
  the analogous warning is tested and which already has the mode table C-1 needs). Two reasons beyond
  the plan's own: `integration_library.rs`'s header declares itself *"Binary-driven integration tests
  for `--library` presets"* and two of the four cases involve no preset at all; and if C-1 is fixed via
  Option B, the tests become a mode table, which belongs next to `-a2`'s, not among preset tests. Keep
  `explicit_clip_flag_overrides_the_preset_and_is_logged`'s existing stderr assertion in
  `integration_library.rs` as the preset-side guard.
- **O-4 — CHANGELOG.** One sentence, matches the file's `- **Bold** ([#NNN](url)) — …` shape: compliant.
  A new entry is defensible. But note #440's existing entry is *itself* still under `### Unreleased`
  (`CHANGELOG.md:4-8`), so the release notes will carry two bullets citing #440 — the plan's rationale
  ("the previous entries describe released-in-this-cycle behaviour") does not distinguish the two cases,
  since neither has shipped. Folding one clause into the existing bullet is the cleaner read; Felix's
  call, and nothing breaks either way.
- **O-5 — `--clip_R2 0`.** Warns as designed (verified: `WARNING: --clip_R2 0 was given but is not used
  …`). Assumption 2 is decided explicitly, which is right. Minor: `0` is a no-op on a paired run too, so
  "is not used in this mode" implies a paired run would use it, which is vacuous for `0`. Not worth code.
- **O-6 — efficiency.** Agreed, nothing on a hot path. `validate()` runs exactly once (`main.rs:242`),
  and the cost is at most one `Vec` of two `String`s. Correct as written; no change.
- **O-7 — the comment.** `// Read 2 clipping cannot apply to a single-end run; the flags are inert here.`
  One line, states the fact, no history, no justification. Compliant. If C-1 is fixed via Option B the
  comment will need to stop saying "single-end run".

---

## Verified as correct (so a missing finding means a check, not an omission)

- **No existing test breaks.** `grep -rn -- "clip_R2\|clip_r2" tests/` hits only
  `tests/integration_library.rs`. The single single-end-plus-explicit-R2 invocation in the whole suite
  is the oracle run at `integration_library.rs:156-169`, which asserts only `assert!(ok, …)` and
  `assert_eq!(read(…), read(…))` — a new stderr line cannot affect either. The stderr assertion at
  `:149-152` (`!err.contains("--clip_R2")`) applies to the *first* run, which sets no raw R2 flag, so it
  stays silent. Both remaining clip_R2 tests (`:258`, `:287`) are `--paired`. *Verified by grep + read
  plus 11 manual binary runs; I did not execute the integration suite with the patch applied.*
- **CI is unaffected.** `grep -n -- "clip" .github/workflows/*.yml` matches only "clippy" — the
  `validation` job passes no clip flag at all, on any shape.
- **`--rrbs` cannot interact** (`main.rs:1141`, paired-gated, local binding).
- **`--clump_only` cannot reach the branch** (`cli.rs:839/848/852`).
- **`--paired` + one interleaved uBAM stays silent** (`self.paired` true; R2 clipping genuinely
  applies after de-interleaving): `trim_galore --dont_gzip --paired --clip_R2 5 -o . ubam_paired_test.bam`
  → `exit=0 warn=0`.
- **Multi-file single-end** (`--clip_R2 5 --three_prime_clip_R2 7 a.fq b.fq`, no `--paired`) → one
  warning for the invocation, correct and not per-file.
- **`--passthrough` cannot interact** — it bails unless `--paired` (`cli.rs:932-934`).
- **The "Ignoring" half of the message is true under the specialty modes** — `src/specialty.rs` reads no
  clip flag, and specialty dispatch returns before every `setup_trimming` call site.
- **V8's mutation is non-vacuous** (by reading, not executed): moving the `--library` block above
  `:1013` puts it above the hardtrim conflict bail (1011), the range bails (1018-1027) and the input-file
  loop (1044-1046), so all three new tests plus the pre-existing `hardtrim5_range_error_…` go red on
  `!err.contains("--library")`.
- **V1/V3/V5/V6 name real, runnable checks** — `run()` in `integration_library.rs:28` returns stderr, so
  steps 5-8 are observable; step 4 is a genuine `Ok`-assertion and passed. The plan's warning about
  multiple `cargo test` filters matching nothing is correct and worth keeping.
- **The plan's claim that a `validate()` `eprintln!` cannot be unit-tested** is right, and reducing the
  unit test to the `Ok`-assertion is the honest response.
