# CODE REVIEW — Reviewer A — `9a738a3` (`fix/440-library-refusals`)

**Target:** commit `9a738a3` — `fix(cli): refuse --library where clipping cannot apply, and gate R2 reporting (#440)`
**Parent:** `872c4ae` (PR #441, merged into `dev`)
**Diff reviewed:** `git diff 872c4ae 9a738a3` — 9 files, +382 / −25
**Plan:** `plans/08162026_pr441-library-presets/PLAN.md` v2 (incl. §12 implementation notes)
**Worktree:** `/private/tmp/claude-501/pr441-A`

```
$ git rev-parse HEAD
9a738a333f1c2c75219b4f15ce2ebff418447495
```

---

# Verdict: **READY TO PUSH**

Every behavioural claim in the commit message and §12 that I checked is true, and all four mutation
rows I re-ran reproduced exactly as recorded. The implementation is sited correctly, the wording is
factually accurate, both `flag_summary` call sites are genuinely covered, and the JSON is untouched.

I found **no defect that should block the push**. I found **one Medium follow-up** (a user-visible
confusion of the same class the commit fixes, surviving in one shape and untested) and **one Medium
test-coverage gap** (the block's *position* is not pinned, only one of the four precedence
relations it depends on). Both are latent, not regressions. Everything else is Low style.

| Priority | Count | Area |
|---|---|---|
| Critical | 0 | — |
| High | 0 | — |
| Medium | 2 | 1 behavioural follow-up (F-1), 1 test gap (F-2) |
| Low | 5 | wording/consistency (S-1 … S-5) |

---

# Defects

## F-1 — Medium — the R2-override log line still asserts an override that cannot apply, one line below the line that says so

**Where:** `src/main.rs:1125-1132` (stderr banner loop), `src/report.rs:283-291` (text report loop).
Deliberate per PLAN §B2 edge case 1 and §Assumptions 4 — so this is a *documented decision*, and I am
reporting it because the decision leaves the commit's own headline problem alive in one shape, not
because it was overlooked.

`--library emseq --clip_R2 5` on a **single-end** run, real binary at HEAD:

```
=== stderr (preset lines) ===
Library preset 'emseq' selected: --clip_R1 10 --three_prime_clip_R1 10
--clip_R2 5 was given on the command line and overrides the emseq preset value 10
=== text report (preset lines) ===
Library preset: emseq
Clipping in force: --clip_R1 10 --three_prime_clip_R1 10
--clip_R2 5 was given on the command line and overrides the emseq preset value 10
```

Read as a pair, this tells a user their `--clip_R2 5` *took effect* (it "overrides the preset value"),
immediately after a line that correctly says no R2 clipping is in force. Neither sentence is false —
the plan's argument for leaving the override line ungated is sound in isolation — but the juxtaposition
invites exactly the wrong conclusion, and this is the same failure mode the commit exists to fix:
report text asserting clipping that cannot apply.

Note it is a strict improvement on the parent, which printed
`Clipping in force: --clip_R1 10 --clip_R2 5 …` — wrong *and* internally consistent. So: not a
regression, an incomplete fix.

**What sharpens this into a real finding rather than a taste call** is that the repo already has a
precedent for the opposite behaviour on an R2-family flag, three sections up the same `validate()`:

```
$ trim_galore -a2 AGATCGGAAGAGC --dont_gzip -o . BS-seq_10K_R1.fastq.gz
WARNING: -a2/--adapter2 was given but is not used in this mode (it applies to Read 2 of a pair,
and this is a single-end run). Ignoring.

$ trim_galore --clip_R2 5 --dont_gzip -o . BS-seq_10K_R1.fastq.gz
(no warning at all)
```

`cli.rs:1056-1057` warns for `-a2` on a single-end run; `--clip_R2` gets nothing, preset or not. That
half is pre-existing and correctly out of scope by the plan's rules, but this commit creates the
adjacency that makes it visible.

**Recommendation:** don't hold the push. Either (a) suppress R2 override lines when `!paired`, or
better (b) file a follow-up to give the clip-R2 family the `-a2` treatment — one warning at
`validate()` time when `--clip_R2` / `--three_prime_clip_R2` is given on a single-end run, which fixes
the bare-flag case and this one together and needs no gating in the report at all. Option (b) is the
one consistent with the file's own idiom.

**Also untested:** no test covers `--library` + an R2 override on a single-end run in either
direction, so whichever behaviour is intended is currently unpinned.

## F-2 — Medium — the refusal block's *position* is not pinned; only one of the four precedence relations it relies on is tested

**Where:** `src/cli.rs:1112-1134` (the block) and `src/cli.rs:2462-2474`
(`hardtrim5_range_error_precedes_the_library_refusal`).

PLAN §B1.1 rests the whole design on the block being **the last fallible check**, and lists four
precedence guarantees that "follow automatically" (§B1 edge cases): `--hardtrim5 0` range error,
`--hardtrim5`+`--hardtrim3` (#386), `--clock`/`--implicon` arity, missing input files.

I confirmed all four hold at HEAD:

```
$ trim_galore --clock --paired --library pbat BS-seq_10K_R1.fastq.gz
Error: --clock mode requires an even number of input files (R1/R2 pairs), got 1
$ trim_galore --hardtrim5 20 --library emseq /nonexistent/nope.fq
Error: Input file not found: /nonexistent/nope.fq
$ trim_galore --hardtrim5 20 --hardtrim3 20 --library emseq BS-seq_10K_R1.fastq.gz
Error: --hardtrim5 and --hardtrim3 cannot be combined in one invocation; run the two trims as
separate invocations, feeding the first trim's output to the second.
```

But only the *range* relation is guarded by a test. §12 **D3** records that V3 was run as an
equivalent mutation (suppressing the range `bail!`) rather than relocating the block — and the
substitution is not equivalent for this property. It proves "the range error is reported"; it does not
prove "the block sits after everything else".

I ran the mutation D3 skipped: moved the block from `cli.rs:1112` up to just after the `--hardtrim3`
range check (before `validate_paired_input` and before the `check_restartable_input` loop).

```
$ cargo test --no-fail-fast      # with the block relocated
EXIT=0
(no FAILED lines)
TOTAL PASSED: 764
```

**The whole 764-test suite stays green**, while the diagnostics degrade for real:

```
=== MUTATED: --clock --library pbat with ONE input ===
Error: --clock does not honour the clip flags; --library is not compatible
=== MUTATED: --hardtrim5 20 --library emseq /nonexistent/nope.fq ===
Error: --hardtrim5 does not honour the clip flags; --library is not compatible
```

A user with a typo'd filename is told about `--library`. That is a user-visible property the plan
claims and nothing defends.

**Recommendation:** one more test, mirroring the existing one, is enough to close the gap that matters
most — a missing-input case, since it is the one an ordinary user hits:

```rust
/// The refusal is the last fallible check, so an unreadable input list is
/// diagnosed before the flag combination.
#[test]
fn missing_input_precedes_the_library_refusal() {
    let cli = Cli::parse_from([
        "trim_galore", "--hardtrim5", "20", "--library", "emseq", "no_such_file.fq",
    ]);
    let err = cli.validate().unwrap_err().to_string();
    assert!(err.contains("Input file not found"), "got: {err}");
    assert!(!err.contains("--library"), "got: {err}");
}
```

Cheap, and it is the assertion that makes the relocation mutation fail. Optionally a second for
`--clock` arity. Not a blocker — the code is correct today.

---

# Style observations (not defects)

- **S-1 — Low — the new block is the only one in its family without an issue marker.**
  `src/cli.rs:1112` reads `// The specialty arms return before setup_trimming, so no preset reaches a
  clip site.` Every sibling carries the number: `cli.rs:1007` `// #386 —`, `cli.rs:1071` `// #421 —`,
  `main.rs:1114` `// #440 —`, `main.rs:384` `// #408 —`, `report.rs:863` `// #440 —`. The tests section
  header does say `#440`, so it is traceable; prefixing `// #440 —` would match the file. Content is
  otherwise fully CLAUDE.md-compliant: one line, states the fact, current state, no justification.

- **S-2 — Low — `CHANGELOG.md:17-19`, "For the same reason" points at the wrong antecedent.**
  The nearest preceding reason is the Tecan/NuGEN clause ("diversity trimming is a different
  algorithm, and it must run *without* `--rrbs`"), which is not the reason the specialty modes refuse
  a preset. The intended antecedent is "A preset is only a clip bundle", two clauses back. Suggest
  "Because a preset is only a clip bundle, `--library` is also refused with …".

- **S-3 — Low — `docs/.../library-presets.md:42` bullet header contradicts its own body.**
  "**No effect in the specialty modes.**" followed by "`--library` is refused … rather than accepted
  and silently ignored." "No effect" is what the commit *stopped* being true. The surrounding bullets
  are "No adapter sequence." / "No `--length`." / "No mode changes." — all about what a preset does not
  *set*, which this one is not. Suggest "**Not available in the specialty modes.**"

- **S-4 — Low — `docs/.../flags.md:14` omits `--clump_only`, which also refuses `--library`.**
  The added clause names the four specialty modes. `--clump_only` refuses it too, from `cli.rs:850`:
  ```
  $ trim_galore --clump_only --library emseq BS-seq_10K_R1.fastq.gz
  Error: --clump_only does not clip; --library is not compatible
  ```
  A flag-interaction page listing four of the five refusals is incomplete. One extra name.

- **S-5 — Low — commit message uses `--` as a dash where the repo is em-dash-conscious.**
  "Both of its call sites are updated -- the stderr banner and the text report -- or the same run…".
  CLAUDE.md flags em-dashes as a content-targeted CI concern in user-facing strings; this is a commit
  message so nothing greps it, but it is the only `--`-as-dash in the message.

## On the commit message's own argument (brief item 4)

Asked to check whether §8.3's reasoning holds. **The conclusion is right; the load-bearing clause is
the second one, not the first.** The message says gating the JSON "would break the block/top-level
agreement invariant added here" — but that invariant is introduced by a test in this same commit, so
as an argument it is circular: a self-imposed invariant cannot constrain the choice that imposed it.
What does constrain it is the second clause, "the top-level `clip_r2` cannot change without a schema
break" — that is a real, external constraint and it settles the matter on its own.

There is also a cleaner argument available than either the plan's or the message's, which I'd suggest
if anyone revisits this: the `library` block records what the **preset resolved to**, while
`Clipping in force:` claims what was **applied**. Those are different propositions, so the asymmetry
is principled rather than merely tolerated — and it is exactly what the changelog's own new wording
("the clipping it resolved to — all four values in JSON … and in the text report the values in force
for the run") already says. No change needed; noting it because it is a stronger footing than R-2.

---

# What I verified as correct

Absence of a finding below means I ran the check, not that I skipped it.

## 1. Correctness of the refusal

- **Sited where the plan says.** `src/cli.rs:1112-1134`, immediately after the #421 `fastqc_requested`
  block closes at `:1110`. Confirmed it is the **last fallible check**: `grep -n "bail!"` shows the
  final four `bail!`s in the file are the new ones (`:1115 :1120 :1125 :1130`), and everything from
  `:1136` to `Ok(())` at `:1164` is `eprintln!` — no `?`, no `bail!`, no `return Err`.
- **All four precedence guarantees hold at HEAD** — command output quoted under F-2 above: #386
  combination first, `--clock` arity first, missing input first; range error first is additionally
  test-covered.
- **No input can reach a specialty arm with a preset still set.** `validate()` runs before
  `ensure_output_dir` and before every dispatch arm; V4 below confirms zero files written.
- **`--clump_only` short-circuits before the block**, as PLAN §B1 last bullet and §O-2 claim —
  verified for both `--clump_only --library` and `--clump_only --hardtrim5 --library` (S-4 quote).

**V4 re-run, all four modes, real binary, fresh empty dirs, `-o .`, function taking `"$@"`:**

```
rc=1 files=0 msg=Error: --hardtrim5 does not honour the clip flags; --library is not compatible
rc=1 files=0 msg=Error: --hardtrim3 does not honour the clip flags; --library is not compatible
rc=1 files=0 msg=Error: --clock does not honour the clip flags; --library is not compatible
rc=1 files=0 msg=Error: --implicon does not honour the clip flags; --library is not compatible
```

Positive control, proving `files=0` is capable of detecting a write:

```
control (no --library): rc=0 files=1   <- the check can fail
BS-seq_10K_R1.30bp_5prime.fq.gz
```

## 2. The wording — verified independently, and it is accurate

The claim to check was that all four modes **do** clip while none reads `--clip_*`.

```
$ grep -n "clip_5prime|clip_r1|clip_r2|three_prime_clip|effective_clips|then clips" src/specialty.rs
107:            record.clip_5prime(start);          # inside hardtrim3
219:            record.clip_5prime(start);          # inside hardtrim3_to_bam
241:/// as `:R1:<umi1>:R2:<umi2>:F1:<fix1>:F2:<fix2>`, then clips:      # --clock
352:/// UMI barcode (`:BARCODE`), then clips R2 by `umi_length` bases.  # --implicon
```

Zero hits for `clip_r1` / `clip_r2` / `three_prime_clip*` / `effective_clips` anywhere in
`specialty.rs`. `hardtrim5` truncates (`record.truncate(keep)` at `specialty.rs:64`), `hardtrim3`
clips 5', `--clock` and `--implicon` clip at fixed positions per their doc comments.

**So "does not clip" would have been false and "does not honour the clip flags" is exactly right.**
Consistent with both neighbouring families: it keeps the `--clump_only` shape
`<mode> <incapability>; <flag> is not compatible` (`cli.rs:834-850`) and it is the clip-flag analogue
of the `-a2` block's per-mode phrasing (`cli.rs:1050-1055`, "`--clock` performs no adapter trimming").

## 3. The gate — both call sites, and `paired` is the right source of truth at each

- **Exactly two call sites**, no others anywhere including docs:
  ```
  src/main.rs:1121:            clips.flag_summary(cli.paired)
  src/report.rs:282:            library.flag_summary(config.paired)
  ```
  (plus the definition at `library.rs:120` and the two new unit tests).
- **They cannot disagree.** `main.rs` dispatch is a clean partition on `cli.paired`:
  `cli.paired && input.len()==1` → `run_paired_ubam_single_file`; `cli.paired` → `run_paired`; both
  `return`/end before the single-end loop at `:1028-1035`. Same shape in `run_ubam_output`
  (`:2041/:2067/:2159`). Every paired arm funnels to `write_paired_reports`
  (`main.rs:1801, 2008, 2429, 2548`) whose config is `paired: true` (`:2628`); the two single-end
  drivers hardcode `paired: false` (`:1459`, `:2249`). There is no path where `run_single_file`
  executes with `cli.paired == true`.
- **Empirically confirmed across every run shape**, banner vs. text report, real binary:

  | shape | banner | report | |
  |---|---|---|---|
  | `--paired --library accel` (2 FASTQ) | 4 flags | 4 flags | AGREE |
  | `--paired --library accel --output-format ubam` | 4 flags | 4 flags | AGREE |
  | `--paired --library accel ubam_paired_test.bam` (interleaved in) | 4 flags | 4 flags | AGREE |
  | `--library accel --output-format ubam` (SE) | 2 flags | 2 flags | AGREE |
  | `--library accel` (SE FASTQ) | 2 flags | 2 flags | AGREE |

  The two-vs-four contrast proves the extraction is content-sensitive, and V6a below proves the banner
  is independently mutable, so `AGREE` is not vacuous.
- **Multi-pair:** `--paired --library accel` over 2 pairs → 2 banners, 4 text reports, all four
  identical and all listing four flags.
- **`--clump_only`:** unreachable, refused earlier at `cli.rs:850` (quoted in S-4).
- **RRBS interaction is safe.** The auto-clip at `main.rs:1136` sits *after* the banner, but it can
  only fire when `clips.clip_r2.is_none()`, and `PresetClips`'s fields are plain `usize`
  (`library.rs:22-27`) so a selected preset always yields `Some` for all four. The banner therefore
  cannot be stale under a preset — this is PLAN assumption 6, and it currently holds.

## 4. What the gate did NOT touch

- **JSON untouched.** The only non-test change in `report.rs` is the `flag_summary(config.paired)`
  call inside `write_report_header`. The JSON writer (`report.rs:866-905`) is byte-unchanged; no key
  added, removed, renamed, or reordered.
- **A run with no `--library` is unreachable by this change, by construction.** `resolve(None, …)`
  returns `preset: None` (`library.rs:153-163`), and *both* consumers are guarded by a preset check —
  `main.rs:1118` `if let Some(preset) = clips.preset` and `report.rs:274-276`
  `&& let Some(preset) = library.preset`. So `flag_summary` is never called and neither line is
  written. The one difference the changelog admits (`"library": null` in JSON) predates this commit.
  Confirmed empirically: a no-preset run has no `Library preset` line and `parameters.library` is
  present-and-null.
- **Override log lines correctly left ungated** — the two loops are outside the gate at
  `main.rs:1125-1132` and `report.rs:283-291`. Design intent holds; see F-1 for the consequence.
- **Docs' asymmetry explanation is factually right.** `mode` really is the discriminator, and there is
  no `paired` key:
  ```
  mode = "single-end"
  library = {"preset":"emseq","clip_r1":12,"clip_r2":10,"three_prime_clip_r1":10,
             "three_prime_clip_r2":10,"overrides":[{"flag":"--clip_R1","preset_value":10,"user_value":12}]}
  top-level clip keys = {'clip_r1':12,'clip_r2':10,'three_prime_clip_r1':10,'three_prime_clip_r2':10}
  ```

## 5. Test quality — four matrix rows re-run from scratch, all four reproduce

I did not take the green run on trust. Each mutation was applied to my worktree, run with
`cargo test --no-fail-fast` (whole suite, no filters), then reverted and `cmp`-verified byte-identical.

| Row | Mutation applied | Result | §12 claim |
|---|---|---|---|
| **V2** | `if self.library.is_some()` → `if false` | `FAILED. 477 passed; 4 failed` — exactly `library_plus_{hardtrim5,hardtrim3,clock,implicon}_rejected` | matches |
| **V6a** | **only** `main.rs:1121` → `flag_summary(true)` | `FAILED. 7 passed; 1 failed` — **exactly** `explicit_clip_flag_overrides_the_preset_and_is_logged` | matches, exactly 1 |
| **V8** | JSON `clip_r2` / `three_prime_clip_r1` **key labels** swapped | `FAILED. 480 passed; 1 failed` — `report::tests::test_json_report_records_all_four_overrides` | matches |
| **V11** | `None => writeln!(… "library": null …)` → `None => {}` | 2 binaries red — `test_json_report_library_is_null_without_a_preset` **and** `a_run_without_a_preset_reports_no_library_block` | matches |

**V6a is the important one and it is exactly as advertised.** The C-2 gap really is closed: reverting
the stderr half alone, leaving `report.rs` gated, kills precisely one test. Half the headline fix is no
longer silently revertible.

Two corroborations of §12's "survived the whole suite before" claims, which I could not test directly
without a parent binary (disk):
- Under **V8** the *only* failure is the new test. `both_reports_record_the_preset_and_its_expanded_values`
  stayed green — because `pbat` is 8/8/8/8, so a key swap is invisible there. **This vindicates the
  four-distinct-values fixture design (T-1 / R-3) empirically**, not just by argument.
- Under **V11** the two failures are attributable to the *new* `contains_key` assertions specifically:
  `serde_json` returns `Value::Null` for a missing key, so each test's pre-existing `is_null()` line
  still passes under the mutation. F-5a is genuinely closed.

**Suite, lints, formatting at HEAD:**

```
$ cargo test --no-fail-fast
EXIT=0 ... 17 test binaries, 0 FAILED
TOTAL PASSED: 764
$ cargo fmt --all -- --check
FMT_EXIT=0
```

**Clippy — and I verified the check was live rather than a cached false green.** The run reported clean
in ~2 s with no `Checking trim_galore` line, which is the pattern that has burned this work before, so
I injected a canary (`pub fn _clippy_canary(v: &Vec<u8>) -> usize`):

```
error: writing `&Vec` instead of `&[_]` involves a new object where a slice will do
error: could not compile `trim-galore` (lib) due to 1 previous error
CLIPPY_EXIT=101
```

Canary removed, `git status` clean, re-run → `CLIPPY_EXIT=0`. **The green is real.**

> **Deviation, disclosed:** CI runs `cargo clippy --all-targets --release`; I ran `--all-targets`
> (debug) because the volume was at 4.4 GB free and the brief forbade a speculative release build. To
> establish the proxy is sound I confirmed there is **no** `cfg(debug_assertions)` or feature-gated
> code in any of the six changed files, so both profiles lint the same items. Residual risk: a
> release-only lint, which I judge negligible here. Flagging it so this is a disclosed gap, not an
> unchecked claim.

### Seams the matrix does not cover

Beyond F-2 (the significant one), two minor ones, both Low and neither worth a test on its own:

- `preset_with_multiple_overrides_round_trips_through_json` asserts the `overrides` array but **not**
  the four resolved `clip_r*` keys, so V8 leaves it green (confirmed — V8's sole failure was the unit
  test). T-1 is closed at unit level only. Adding the four key assertions to the integration test
  would give it binary-level coverage for ~4 lines.
- The four `library_plus_*_rejected` tests call `cli.validate()` directly, so nothing in this commit
  pins that the refusal reaches a user through `main()`. Covered transitively by
  `an_unknown_preset_name_is_refused`, which runs the real binary and gets a `validate()` error, plus
  my V4 above. Noted for completeness only.
- `test_clump_only_rejects_library` (`cli.rs:1747`) asserts only `err.contains("--library")`, so it
  would stay green if the new block ever preceded the `--clump_only` block and changed the message.
  Only reachable via a doubly-invalid invocation; not worth tightening.

### Test-quality positives worth recording

- All 12 new tests plus the 4 amended assertions are attributable to a mutation that kills them
  (4 verified by me directly, the rest by §12's matrix which reproduced 4/4 where I sampled it).
- Keeping `test_json_report_records_library_block` (cardinality 1) was the right call, and V13's
  recorded result confirms it is the only detector of the comma branch.
- Doc comments on the new tests are genuinely useful, notably the `--implicon` one at `cli.rs:2445`
  explaining `require_equals`, which otherwise reads like a missing argument.

## 6. Repo conventions

- **`cargo fmt --check`** clean; **clippy `-D warnings`** clean (canary-verified live).
- **CLAUDE.md comment rule — compliant throughout.** I read every added comment in the diff. All are
  1-2 lines, state the fact, describe current state, and none narrates the state being fixed or
  pre-empts an objection. Best of them:
  `// \`Value::Null\` is also what a missing key yields, so pin the key's presence.` — fact, one line,
  explains a non-obvious assertion. The rationale that CLAUDE.md sends to the commit message (the
  "does not clip" argument, the JSON asymmetry) is indeed in the commit message and **not** in the
  source. Only nit is S-1, a missing issue marker.
- **`library.rs:117-119`** — the two justification sentences the plan asked to be removed from
  `flag_summary`'s doc comment are gone; the replacement is two lines of fact. PLAN step 1 satisfied,
  and it closes the part of code-review B's C-1 that §11 promised.
- **Commit message fact-checked.** Long by design per CLAUDE.md; I checked it for *errors* rather than
  length and found none:
  - "the only place a preset reaches a clip site" — accurate. `effective_clips()` is also called at
    `main.rs:391` (#408 guard) and the three report configs, but the values that drive trimming come
    only from the `trimmer::TrimConfig` built inside `setup_trimming` (`main.rs:1211`).
  - "Only `--clump_only` rejected it" — accurate (`cli.rs:850`).
  - the `specialty.rs` clipping evidence — accurate, reproduced above.
  - the `--phred64` / #358 precedent, "a mode-dependent rule is worse to document and test" — accurate,
    it is nearly verbatim from `main.rs:274-277`, and that comment likewise admits it "does remove
    working invocations".
  - "12 new tests, 764 passing" — accurate, 764 measured.
  - The one soft spot is the JSON-asymmetry argument; see the note under Style above. It is a weak
    first clause supporting a correct conclusion, not an error.

## 7. The changelog and docs

- **#441's entry corrected in place, no second entry** — confirmed; the hunk is 5 `-` lines / 8 `+`
  lines inside the existing paragraph, and every removed sentence reappears intact in the added text.
  The falsified clause ("record the preset name and **all four values** in force") is now
  "record the preset name and the clipping it resolved to — all four values in JSON …, and in the text
  report the values in force for the run, so a single-end report lists Read 1 only" — **which is true**,
  and matches the binary's real output (verified in V15 below).
- **Nothing lost.** The unique issue-reference set is *identical* to the parent's:
  ```
  $ comm -3 <parent refs> <new refs>
  (empty — identical sets)
  ```
  **52 unique linked GitHub issue/pull refs**, matching the brief's expected number exactly
  (49 of them TrimGalore's own, 49 in the parent too). For the record, a bare `grep -oE '#[0-9]+'`
  returns 70, but 18 of those are unlinked or artefacts — `#3529` is a real MultiQC reference in prose,
  `#8722` is the `&#8722;` minus-sign entity. The URL-based count is the meaningful one.
- **Docs match the binary — V15 re-run independently.** Every line quoted in `library-presets.md`
  reproduced verbatim from a real run of `--library emseq --clip_R1 12`:
  ```
  PRESENT  Library preset 'emseq' selected: --clip_R1 12 --three_prime_clip_R1 10
  PRESENT  --clip_R1 12 was given on the command line and overrides the emseq preset value 10
  PRESENT  Library preset: emseq
  PRESENT  Clipping in force: --clip_R1 12 --three_prime_clip_R1 10
  PRESENT  --clip_R1 12 was given on the command line and overrides the emseq preset value 10
  stderr --clip_R2 hits: 0
  report --clip_R2 hits: 0
  ```
  Harness control, proving the check can report MISSING:
  ```
  === a line that must NOT be present ===
  MISSING (good - the check can fail)
  ```
  `grep -qF --` used throughout, per the known `--`-as-option trap.
- **The docs' JSON block (`library-presets.md:76-83`) matches the real JSON exactly** — all four keys,
  the single override entry, and the values 12/10/10/10 (output quoted in §4 above).
- **`flags.md:14`** added clause is accurate (PLAN step 23 / O-4 implemented per its own
  recommendation, §12 D2). See S-4 for the one omission.
- One observation, not a defect: on a single-end run the **top-level** `clip_r2` / `three_prime_clip_r2`
  keys are still present and non-null (10 in the transcript above), so a naive consumer reading
  `parameters.clip_r2` sees a value that was not applied. Pre-existing, unchanged by this commit, and
  the new docs paragraph plus `mode` is the documented answer. Mentioning it only so it is on record
  as considered.

---

# Worktree state

Every mutation reverted and verified byte-identical against a `_snap/` copy taken before any edit
(`cli.rs`, `main.rs`, `report.rs` — `cmp` clean each time). Scratch artefacts removed. Final state,
command run rather than asserted:

```
$ git rev-parse HEAD
9a738a333f1c2c75219b4f15ce2ebff418447495
$ git status --short
(clean)
```

`cargo clean` run in this worktree on completion, per the disk warning.

**No GitHub writes of any kind were made.** No push, no comment, no PR. Nothing outside
`/private/tmp/claude-501/pr441-A` was modified except this report file.
