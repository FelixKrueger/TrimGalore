# PLAN REVIEW — Reviewer A (generation g1)

**Plan:** `plans/08172026_clip-r2-se-warning/PLAN.md`
**Worktree:** `/private/tmp/claude-501/prA440` @ `fdade7efb8b579fd9b25899d04b1609dba860b20` (current `dev` tip — verified)
**Method:** every anchor read at `file:line`; the plan's step-1 snippet applied, built, linted, run against
8 invocations, the V4 mutation applied and reproduced, then reverted and the revert verified.

---

## Verdict — **READY WITH CHANGES**

The plan's core is right and I could not break it. The central design decision (read the raw
`self.clip_r2`, never `effective_clips()`) is **verified correct and load-bearing** — I reproduced the
exact regression it guards. All three F-2 precedence relations hold at `fdade7e` and the three proposed
tests will pass as written, including the `--clock` case the brief flagged as a possible slip. The
snippet compiles, `clippy -D warnings` clean, `fmt` clean, and produces exactly the four documented
outcomes.

What needs changing before implementation:

1. **C-1** — F-1 is only **half** closed. The misleading pair the plan exists to fix survives
   **verbatim in `_trimming_report.txt`**, which is the artifact that persists; the plan lists "touches
   no report" as a virtue. Either fold in the two-line gate or say plainly that the report half stays
   open. Reproduced below.
2. **I-1** — the snippet emits ungrammatical output in the two-flag case: *"--clip_R2 5 and
   --three_prime_clip_R2 3 **was** given"*. Reproduced.
3. **I-2** — the `-a2` convention the plan invokes has a **mode-accurate-reason** rule with a test
   pinning it (`tests/integration_adapter2.rs:369`). The new warning fires on single-end **specialty**
   runs and states the single-end reason there. This is *not* Q-2 (which is about extending coverage);
   it is the reason text being wrong on a shape the warning already covers.
4. **I-3** — the test baseline is **768**, not 764. V9 should read 776.
5. **I-4** — steps 5–8 omit `-o`, so they would write into tracked `test_files/`. Every existing test
   in that file passes `-o <tempdir>`.
6. **I-5** — step 12's docs target **does not exist**, and its stated precedent is false: `-a2`'s
   warning is documented nowhere in `docs/`.
7. **I-6** — three behaviours the plan decides in prose (`--three_prime_clip_R2` alone, the combined
   two-flag line, `--clip_R2 0`) have no test and no V-row. I-1 is exactly what that gap hides.

---

## Claim-verification table

| # | Claim | Result | Settled by |
|---|---|---|---|
| 1a | `cli.rs:1013` — #386 hardtrim5+hardtrim3 bail | **VERIFIED** (drift 1) | `cli.rs:1011-1017`; `bail!` at `:1012`, quoted text at `:1013` |
| 1b | `:1021` / `:1026` — hardtrim5 / hardtrim3 range bails | **VERIFIED** exact | `cli.rs:1021`, `cli.rs:1026` |
| 1c | `:1029` / `:1032` — `validate_paired_input` for `--clock` / `--implicon` | **VERIFIED** exact | `cli.rs:1028-1033` |
| 1d | `:1045` — `check_restartable_input(path, "Input file not found")` | **VERIFIED** exact | `cli.rs:1043-1046` |
| 1e | `:1048-1070` — the `-a2` warning block, wording as quoted | **VERIFIED** | `cli.rs:1048-1071`; message at `:1064-1066`, reason at `:1059` |
| 1f | `:1114` — the `--library` refusal | **VERIFIED** | comment `:1114`, `if self.library.is_some()` `:1115`, bails `:1117/1122/1127/1132` |
| 1g | `:1138` — the deprecation block, proposed insertion point | **VERIFIED** exact | `cli.rs:1138` `// Deprecation warnings for Perl-era flags` |
| 1h | `:2480` — `hardtrim5_range_error_precedes_the_library_refusal` | **VERIFIED** exact | `cli.rs:2480` |
| 1i | `:2565` — the #421 ordering-test section | **VERIFIED** exact | `cli.rs:2565` `// ── #421: the guards sit below the value and input-list checks ──` |
| 1j | `:837` — the `--clump_only` `--clip_r2` bail | **REFUTED** (drift 2) | `cli.rs:839` is the bail; `if self.clip_r2.is_some()` is `:838`. `:837` is a closing brace |
| 1k | a `── #440` test section exists to place F-2's tests in | **VERIFIED** | `cli.rs:2417` `// ── #440: --library refused where no clip site is reachable ──` |
| 2 | "There is no warning today", zero matches in `cli.rs`/`main.rs`, nothing equivalent elsewhere | **VERIFIED** | `grep -rn "WARNING" src/*.rs \| grep -i clip` → **zero matches** (exit 1). `grep -n "clip_r2\|clip_R2" src/main.rs` → only the RRBS auto-set note and config plumbing |
| 3 | A preset always populates all four `Option`s, so `effective_clips()` would warn on every SE preset run | **VERIFIED** | `library.rs:167-178` — `pick` returns `Some(v)` in the `Some` arm and `Some(preset_value)` in the `None` arm, unconditionally. Reproduced live: see V4 below |
| 4 | `:1138` is after every `bail!` in `validate()` | **VERIFIED** | `grep -n "bail!" src/cli.rs \| tail` — last is `:1132`; `:1139-1164` are pure `eprintln!`, `Ok(())` at `:1166`. No `?` either |
| 4b | "so a refused run stays quiet" | **PARTLY REFUTED** | true for `validate()`; `main.rs` bails at `:290, :319, :330, :340, :347, :403` and in every `preflight_output_collisions` call, all after `cli.validate()?` at `main.rs:242` |
| 5 | `--rrbs` auto-`clip_r2 2` is in `main.rs::setup_trimming` and paired-only | **VERIFIED** | `main.rs:1141` — `if cli.rrbs && !cli.non_directional && cli.paired && clips.clip_r2.is_none()`. Writes a local, never `self.clip_r2`, so it is invisible to the raw-field check |
| 6 | `--clump_only --clip_R2` is unreachable | **VERIFIED** | `cli.rs:838-840` bails `--clump_only does not clip; --clip_r2 is not compatible` (line number differs — see 1j) |
| 7a | test 9 gets `"cannot be combined in one invocation"` | **VERIFIED** live | `Error: --hardtrim5 and --hardtrim3 cannot be combined in one invocation; …` |
| 7b | test 10 hits the arity error, not the early `Ok` | **VERIFIED** live | `mode_label` is `"--clock"` at `cli.rs:1029`, and the early `Ok` at `:556` requires `mode_label == "Paired-end"` → falls to `:559`. Binary: `Error: --clock mode requires an even number of input files (R1/R2 pairs), got 1` |
| 7c | test 11 gets `"Input file not found"` (assumption 5) | **VERIFIED** live | `Error: Input file not found: no_such_file.fq.gz` |
| 7d | no clap-level `conflicts_with` aborts `parse_from` for the three invocations | **VERIFIED** | `grep -n "conflicts_with" src/cli.rs` → only the five adapter presets and `no_poly_g`. None on hardtrim/clock/paired/library |
| 8 | baseline "764" | **REFUTED** | `cargo test --no-fail-fast` on `fdade7e`: **768 passed, 0 failed, 0 ignored, 17 binaries**, `EXIT=0`. #446 (`fdade7e`) is the +4 |
| 9a | the snippet compiles, clippy-clean, fmt-clean | **VERIFIED** | `cargo build` ok; `cargo fmt --all -- --check` → FMT CLEAN; `cargo clippy --all-targets -- -D warnings` → `CLIPPY_EXIT=0`, zero warning/error lines |
| 9b | the four documented stderr outcomes | **VERIFIED** all four | table below |
| 9c | revert verified | **VERIFIED** | `git status --short` → empty; `git diff --stat HEAD` → empty; HEAD still `fdade7e`. Empty status also proves `test_files/` was not polluted |
| — | `run()` in `integration_library.rs` returns stderr | **VERIFIED** | `tests/integration_library.rs:28-37` returns `(status.success(), stderr)` |
| — | `validate()` runs once per invocation | **VERIFIED** | single call site, `main.rs:242` |

**Drift 1** (`:1013` → `bail!` at `:1012`) and **Drift 2** (`:837` → `:839`) are cosmetic; both point
into the right block. Worth fixing so V8's "move the refusal block above `cli.rs:1013`" is unambiguous.

### 9b — the behavioural matrix (patched debug binary, `-o <scratch>`, `--dont_gzip`)

```
se_clip_r2                 exit=0  warnlines=1
      > WARNING: --clip_R2 5 was given but is not used in this mode (it applies to Read 2 of a pair, and this is a single-end run). Ignoring.
se_preset_clip_r2          exit=0  warnlines=1
      > WARNING: --clip_R2 5 was given but is not used in this mode (…single-end run). Ignoring.
se_preset_only             exit=0  warnlines=0
pe_clip_r2                 exit=0  warnlines=0
--- extra probes ---
se_both_r2_flags           exit=0  warnlines=1
      > WARNING: --clip_R2 5 and --three_prime_clip_R2 3 was given but is not used in this mode (…). Ignoring.
se_clip_r2_zero            exit=0  warnlines=1
      > WARNING: --clip_R2 0 was given but is not used in this mode (…). Ignoring.
se_hardtrim_clip_r2        exit=0  warnlines=1
      > WARNING: --clip_R2 5 was given but is not used in this mode (it applies to Read 2 of a pair, and this is a single-end run). Ignoring.
pe_hardtrim_clip_r2        exit=0  warnlines=0
```

All four plan-specified cases behave exactly as documented. The four extra probes are where the
findings come from.

---

## Findings

### CRITICAL

#### C-1 — F-1's misleading pair survives verbatim in the trimming report, and the plan presents that as a feature

The plan's Goal is F-1: `--library emseq --clip_R2 5` on a single-end run reads as though the 5 took
effect. Reviewer A's finding quoted **two** artifacts — stderr *and* the text report. The plan fixes
neither of them; it prepends a warning elsewhere in the stream. Patched binary,
`--library emseq --clip_R2 5 --dont_gzip -o <dir> BS-seq_10K_R1.fastq.gz`:

```
=== stderr ===
1:WARNING: --clip_R2 5 was given but is not used in this mode (…single-end run). Ignoring.
11:Library preset 'emseq' selected: --clip_R1 10 --three_prime_clip_R1 10
12:--clip_R2 5 was given on the command line and overrides the emseq preset value 10

=== BS-seq_10K_R1.fastq.gz_trimming_report.txt ===
12:Library preset: emseq
13:Clipping in force: --clip_R1 10 --three_prime_clip_R1 10
14:--clip_R2 5 was given on the command line and overrides the emseq preset value 10
```

Two things follow, neither acknowledged in the plan:

- **The report file contains the misleading pair and no warning at all.** `grep "not used in this mode"`
  over the report returns nothing. The report is the artifact a collaborator opens next month; stderr is
  gone by then. `report.rs:284-293` writes the override loop with **no `paired` gate**, while
  `flag_summary(config.paired)` on `:282` gates the R2 values out of `Clipping in force`. That
  asymmetry inside one block is the whole of F-1, and it is untouched.
- **On stderr the pair is still adjacent** (lines 11–12) and the warning is ten lines *above* it, at
  line 1 — so the caveat is read before the claim it qualifies, not with it.

Reviewer A's F-1 offered two options: (a) suppress R2 override lines when `!paired`, or (b) the
`validate()` warning, with the claim that (b) "needs no gating in the report at all". **The artifact
above refutes that claim**, and the plan inherits it — §Integration lists "**Touches no report, no
JSON, no output file**" as a virtue when it is in fact the boundary of what the fix achieves.

For completeness, the JSON keeps recording the value as in force on a single-end run:

```json
"library": { "preset": "emseq", "clip_r1": 10, "clip_r2": 5, "three_prime_clip_r1": 10,
             "three_prime_clip_r2": 10,
             "overrides": [{"flag": "--clip_R2", "preset_value": 10, "user_value": 5}] }
```

(`report.rs:863-865` describes these as "the ones in force". Whether the JSON should change is a
separate machine-readability argument I am **not** making — MultiQC consumers may want the raw record.)

**What to do — pick one and say so in the plan:**

1. **Preferred:** fold in the two-line gate. Filter the override loops at `main.rs:1128` and
   `report.rs:284` so R2 overrides are skipped when `!paired`, keyed off the same read-2 flag the
   `flag_summary` table already uses at `library.rs:122-127` — ideally by adding `read_2: bool` to
   `ClipOverride` rather than string-matching the flag name. This genuinely closes F-1. Note it **changes
   report bytes** for SE preset+override runs, so it needs its own test and it invalidates V10 as
   written — that is a real scope increase, not a free two lines.
2. **Acceptable:** keep the plan stderr-only, but replace the "touches no report" line with an explicit
   statement that the report half of F-1 remains open, and carry it back to
   `08162026_pr441-library-presets/PLAN.md` §13 as a third open item so it is not silently dropped.

Ranked Critical not because any code here is wrong, but because the plan would be signed off as
"closes F-1" when the persistent half of F-1 is untouched and the plan reads as though it is not.

### IMPORTANT

#### I-1 — the snippet's two-flag line is ungrammatical

Reproduced verbatim from the plan's own snippet:

```
WARNING: --clip_R2 5 and --three_prime_clip_R2 3 was given but is not used in this mode …
```

`inert.join(" and ")` builds a plural subject under a singular verb. Assumption 3 (one combined line)
is fine; the sentence has to agree with it. Either pick the verb by `inert.len()`, or restructure so
the number does not matter, e.g. `"WARNING: {} cannot be used in this mode (…). Ignoring."` — which
also drops the `-a2` sentence's "was given", the part that does not survive pluralisation.

This is the concrete cost of I-6: the plan decides the combined-line behaviour in prose and then never
tests it, so the defect ships.

#### I-2 — the warning's reason is not mode-accurate, and the convention it cites has a test for exactly that

The plan's argument is that this is "the R2 family's missing half of an established convention". That
convention's rule is stronger than the plan's wording captures. `cli.rs:1052-1062` picks a reason by
mode — hardtrim, then clock, then implicon, then single-end — and
`tests/integration_adapter2.rs:369`, *`unusable_a2_warns_with_a_mode_accurate_reason`*, is a
table-driven test over all five whose assertion is `stderr.contains(reason)` with the comment "a reason
true for that mode".

The plan's snippet has one reason. Measured:

```
se_hardtrim_clip_r2   → WARNING: --clip_R2 5 … (it applies to Read 2 of a pair, and this is a
                        single-end run). Ignoring.
pe_hardtrim_clip_r2   → (silent)
```

On the first line the operative reason is that `--hardtrim5` performs no clipping at all; the stated
reason implies that adding `--paired` would make the flag apply — and the second line shows it would
not, it would just go quiet. A user following the warning's own logic ends up worse off.

**This is not Q-2.** Q-2 is about *extending* the warning to shapes it does not currently cover, and I
accept that exclusion. This is about the reason text being wrong on a shape the new warning **already
fires on**. It is fixable inside the plan's scope by giving `inert` the same `unusable_reason` chain
`-a2` uses (about four lines), and it is what makes deferring Q-2 coherent: with a mode-accurate reason
the specialty-paired gap is a *gap*; with a single reason it is an *incorrect statement*.

#### I-3 — the baseline is 768, not 764

```
$ cargo test --no-fail-fast          # on fdade7e
768 passed; 0 failed; 0 ignored; 17 binaries; EXIT=0
```

The lib target alone is 485. The +4 is `fdade7e` itself (`#446`). V9 should read **776 (768 + 8)** on
macOS, with the plan's existing caveat that Linux differs. The "+8" arithmetic is correct: 1 unit test
(step 2) + 4 binary (step 3) + 3 ordering (step 4).

#### I-4 — steps 5–8 must pass `-o <tempdir>`; as written they write into `test_files/`

The plan's step 3 gives argument lists with no `-o`. Trimming writes beside the input, so
`run(&["--clip_R2", "5", R1])` drops `BS-seq_10K_R1_trimmed.fq.gz` and two report files into the
tracked fixture directory. All eight existing tests in `tests/integration_library.rs` pass
`-o <tempdir>` (lines 68, 86, 113, 138, 166, 186, 229, 265, 302), including the refusal test at
`:250` that produces no output at all. The plan's §Validation *does* say "Run binary cases with `-o .`
in a scratch dir" — but that instruction is in the V-block, not in the step that will be copied into
the test file. Spell `-o` into steps 5–8 and reuse the file's `tempdir(tag)` helper.

#### I-5 — step 12's docs target does not exist, and its precedent is false

There is no `--clip_R1/R2` line in `docs/src/content/docs/guide/flags.md`. The file is not a flag
reference: its sections are `## Cross-flag interactions`, `## RRBS-specific guidance`,
`## Adapter specification recap`, `## Perl-era flag forms`, `## Where to find what`. The only clip
mentions are the `--rrbs` bullet (`:13`) and the `--library` bullet (`:14`).

And the stated precedent — "matching how `-a2` is already described" — is not there:

```
$ grep -rn "not used in this mode" docs/src/content/docs/    → (no matches)
$ grep -rn "adapter2" docs/… | grep -i "single\|ignor\|not used"  → (no matches)
```

`-a2`'s single-end warning is documented **nowhere** in the docs, so there is nothing to match.

**Recommendation:** add a new bullet to `## Cross-flag interactions` in `flags.md`, in the shape of the
`--fastqc` bullet already there (which documents a refusal the same way). Also note `--help` already
carries the fact — `cli.rs:125` reads "Remove N bases from the 5' end of Read 2 (paired-end only)" — so
the docs clause is reinforcement, not the primary record. If a second home is wanted,
`docs/src/content/docs/guide/paired-end.md` already discusses `clip_R2`.

#### I-6 — three decided behaviours have no test and no V-row

| Decided in | Behaviour | Test | V-row |
|---|---|---|---|
| B1 edge case 1 / assumption 3 | both flags → **one** line listing both | none | none |
| B1 edge case 2 / assumption 2 | `--clip_R2 0` still warns | none | none |
| step 1 code | `--three_prime_clip_R2` alone warns | none — steps 5–8 all use `--clip_R2` | none |

Half the feature (`--three_prime_clip_R2`) is never exercised by any proposed test, so deleting its
`if let` block would keep all 776 tests green. Add one case covering both flags together (which pins
the combined line *and* would have caught I-1) and one for `--three_prime_clip_R2` alone. `--clip_R2 0`
is cheap to fold into the same table.

#### I-7 — assumption 4 over-claims

"Placing it after all bails is right, so **a refused run stays quiet**" is true only of `validate()`.
`main.rs` refuses after `cli.validate()?` (`:242`) at `:290, :319, :330, :340, :347, :403` and in every
`preflight_output_collisions` call — so an output-collision refusal, a §3.4b `--preserve-tags`
refusal, or an aligned-BAM rejection all print this warning first. The **placement is still right** and
the divergence from `-a2` is still an improvement; only the sentence needs narrowing to "no
`validate()`-level refusal also emits this warning".

### OPTIONAL

- **O-1 — V1's command.** The "How" cell is `cargo test explicit_r2_clip_on_single_end_warns preset_plus_explicit`
  followed by the instruction to run separately. Just make it two commands; the current cell is the
  broken form a reader will paste. (The filters themselves are fine — cargo does substring matching, so
  `preset_plus_explicit` matches test 6.)
- **O-2 — V10 is weak but not vacuous.** It guards the plausible wrong implementation of zeroing
  `self.clip_r2`, so keep it. But run it on the **preset** SE shape too, which is the report-richest
  case — and note that if C-1 option 1 is taken, V10's "both reports byte-identical" inverts and must
  become an assertion about *which* lines changed.
- **O-3 — V8's mutation is sound.** I traced it: moving the refusal block above `:1013` makes all three
  new tests fail on `!err.contains("--library")` — test 9 and 11 via the `hardtrim5` arm (`:1117`),
  test 10 via the `clock` arm (`:1127`). A real killing mutation for all three.
- **O-4 — the comment is convention-compliant.** One line, present tense, no history. Mildly redundant
  — "cannot apply … the flags are inert here" says it twice; `// Read 2 clipping cannot apply to a
  single-end run.` is enough.
- **O-5 — CHANGELOG.** A new entry is the right call, and the proposed text is one sentence in the
  file's `- **Bold lead** ([#N](url)) — …` shape. Specify that it goes under the existing
  `#### Changes` subheading of `### Unreleased` (`CHANGELOG.md:4-6`), not directly under `### Unreleased`.
- **O-6 — efficiency: agreed, nothing on a hot path.** Once per invocation, at most one two-element
  `Vec`. Worth noting the raw-field design is also the cheaper one: `effective_clips()` allocates a
  `Vec<ClipOverride>` per call and the plan's V4 mutation calls it twice.
- **O-7 — Q-1 recommendation: a new file.** The plan offers `integration_library.rs` vs a new file. The
  decisive fact is that neither option is where the precedent lives: `tests/integration_adapter2.rs`
  holds `unusable_a2_warns_with_a_mode_accurate_reason` (`:369`) and
  `single_end_does_not_announce_a_read2_adapter` (`:424`) — this exact test shape, table-driven over
  modes. Three of the four proposed cases involve no preset, and with I-2 and I-6 folded in the case
  list grows to six or seven, which is a table, not four ad-hoc tests. **Recommend
  `tests/integration_clip_r2_warning.rs`**, copying `integration_adapter2.rs`'s table-driven layout and
  its `binary()`/`tempdir()`/`run()` helpers, and leaving in `integration_library.rs` only
  `preset_plus_explicit_r2_clip_on_single_end_warns` and `preset_alone_on_single_end_does_not_warn`,
  which are genuinely about preset adjacency. Cost: ~15 duplicated helper lines, which is already the
  convention across all 14 integration files.
- **O-8 — Q-2's exclusion is coherent, with one caveat.** Deferring the *paired* specialty coverage is
  defensible: it is Reviewer B's F-5, a scope decision Felix owns, and nothing in this plan depends on
  it. The caveat is I-2 — the reason-accuracy half cannot be deferred with it, because the warning
  already fires on single-end specialty runs today.
- **O-9 — anchor drift.** `:837` → `:839`; `:1013` → `bail!` at `:1012`. Fix so V8's instruction is
  unambiguous.

---

## Stress-testing the central design decision (the brief's item 3)

**The decision is correct.** `library.rs:167-178`: `pick` returns `Some` in *both* arms, so a preset
makes all four fields `Some` unconditionally. I applied the V4 mutation (`self.effective_clips().clip_r2`
in place of `self.clip_r2`) and ran a bare `--library emseq` single-end run:

```
WARNING: --clip_R2 10 and --three_prime_clip_R2 10 was given but is not used in this mode …
```

Flags the user never typed, on every single-end preset run. Test 7's assertion
(`!stderr.contains("is not used in this mode")`) would go red. **V4 is a real killing mutation, not a
comment** — the plan's most important claim about its own validation holds.

I could not find an invocation where the raw fields give the wrong answer:

| Invocation | Raw-field behaviour | Correct? |
|---|---|---|
| `--rrbs --clip_R2 2 <se>` | warns | ✅ the RRBS auto-clip is paired-only (`main.rs:1141`) and writes a local, so nothing is hidden |
| `--rrbs --paired … ` (auto-set fires) | silent | ✅ `self.clip_r2` stays `None`; `self.paired` is true anyway |
| `--paired` + one interleaved uBAM | silent | ✅ `self.paired` is true and R2 exists after de-interleaving (`validate_paired_input`'s early `Ok`, `cli.rs:556`) |
| `--clip_R2 5 a.fq b.fq c.fq` | warns **once** | ✅ `validate()` runs once (`main.rs:242`); one line for the run, not per file |
| `--clip_R2 5 R1.fq R2.fq` (forgot `--paired`) | warns | ✅ **the best case for the feature** — two SE runs, and the warning is the one signal that the pair was not paired |
| `--clump_only --clip_R2 5` | unreachable | ✅ bails at `:839` |
| `--library emseq <se>` | silent | ✅ |

The only shape where a user might expect a warning and not get one is `--paired` + a specialty mode —
which is Q-2, and out of scope by construction.

**Is `validate()` the right layer?** Yes. It is where `-a2` warns, it runs exactly once before any I/O,
and the check needs no format knowledge (it keys on `self.paired`, not on a detected input shape, so
the uBAM interleaved case is handled for free). `main.rs` would be the wrong layer for the warning —
though note C-1's report gate, if taken, does belong in `main.rs`/`report.rs`.

---

## Verified correct — so absence of a finding is not absence of a check

- Every one of the eleven anchors in the brief, read at `file:line`. Nine exact, two off by 2 (O-9).
- `:1138` is after the last `bail!` (`:1132`) and after the last `?` in `validate()`; the divergence
  from `-a2`'s placement is deliberate and better (narrow the claim per I-7).
- No pre-existing warning anywhere in `src/` mentions a clip flag — the plan's "no warning today" is
  exactly right, including in `main.rs::setup_trimming`.
- `--rrbs`' auto-`clip_r2 2` cannot reach the `!self.paired` branch, and does not write `self.clip_r2`.
- `--clump_only` + any of the four clip flags bails before the insertion point.
- All three F-2 relations verified **twice** — by reading control flow, and by running the binary at
  `fdade7e`. The `--clock --paired --library pbat R1` case the brief flagged as a possible slip does
  **not** slip: the early `Ok` at `cli.rs:556` requires `mode_label == "Paired-end"`, and `:1029`
  passes `"--clock"`.
- Assumption 5 confirmed: the message is exactly `Input file not found` (`cli.rs:1045`,
  formatted `"{not_found}: {path}"` at `:524`).
- No clap `conflicts_with` will abort `parse_from` on the three F-2 invocations.
- The snippet builds, `fmt` clean, `clippy --all-targets -- -D warnings` exit 0 with zero diagnostics.
- All four plan-specified stderr outcomes reproduced exactly; exit code 0 in every warning case (V6).
- V4 and V8 are both real killing mutations. V2, V3, V5, V6, V7 name tests that will exist and will
  run. No vacuous row found — V10 is weak (O-2), V1's command form is wrong but its filters do match
  (O-1).
- Revert verified by command, not assertion: `git status --short` empty, `git diff --stat HEAD` empty,
  HEAD still `fdade7e`. The empty status is also the proof that no run polluted `test_files/`.

**Not checked:** release-profile clippy (disk at 11 GiB free after the debug build; debug clippy with
`-D warnings` was clean, and the snippet is profile-independent). The `--library` + interleaved-uBAM
paired shape was reasoned from `cli.rs:556` rather than run against a uBAM fixture. `cargo clean` was
run in the worktree on completion.
