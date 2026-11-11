# PLAN REVIEW — Reviewer B (generation g1)

**Plan:** `plans/08162026_pr441-library-presets/PLAN.md`
**Base verified:** `/private/tmp/claude-501/pr441-B` at `145432e179ae614814dd9754846b76c77137127b`
**Tree state at start and end of review:** tracked tree clean; only the pre-existing untracked
leftovers (`mutM.log`, `mutcheck.sh`, `testsummary.log`, `verify1.log`).
**Reviewer:** B. Independent — the other reviewer's report was not read.

```
$ cd /private/tmp/claude-501/pr441-B && git rev-parse HEAD && git status --short
145432e179ae614814dd9754846b76c77137127b
?? mutM.log
?? mutcheck.sh
?? testsummary.log
?? verify1.log
```

The release binary already present in the worktree self-identifies as the right base, so every
runtime observation below is from `145432e`:

```
$ ./target/release/trim_galore --implicon=8 --paired --library pbat -o … R1 R2 2>&1 | head -5

Trim Galore v2.3.0
145432e — macos/aarch64 — built 2026-08-16T11:48:44Z
```

---

## Verdict

**READY WITH CHANGES.**

The plan is unusually well-grounded: every line-number anchor I checked resolves to what it
says it does, the `flag_summary` two-call-site catch is correct and load-bearing, the step-10
"this existing test must be updated" catch is correct, and the four-distinct-values test lever
does what the plan claims. Nothing needs rework.

Four things must change before implementation, none of them structural:

1. **The refusal message is factually wrong for `--clock` and `--implicon`.** The repo's own doc
   comments describe both modes as clipping. `"--clock does not clip"` contradicts
   `specialty.rs:241-243`. There is an in-file precedent (`cli.rs:1050-1055`) that fixes this
   inside Felix's chosen wording family.
2. **Step 4's `--implicon 8` test does not parse the way the plan assumes.** Proven, not
   guessed. Assumption 6 can be closed now: `--implicon=8` (or bare `--implicon`) is required.
3. **The stderr log line has no automated test at all** — the exact desync the plan calls "the
   load-bearing detail" is checked only by V7, a manual eyeball. One assertion closes it.
4. **#441's own CHANGELOG entry becomes false** ("Both trimming reports record the preset name
   and all four values in force"). Both entries ship in the same `Unreleased` section, so the
   release ships a self-contradiction. This is factual correction, not the editorial trimming
   Felix excluded.

Plus one coverage regression the plan introduces (the single-element `overrides` JSON array
loses its only test) and two validation checks that cannot fail as written (V5's first command,
V15's "byte-identical").

---

## Claim-verification table

Every row was settled by reading or running the code, not by trusting the plan.

### Line-number anchors

| Plan claim | Verdict | Settled by |
|---|---|---|
| `cli.rs:1071-1097` = #421 `--fastqc` refusal on the four modes | **VERIFIED** | `cli.rs:1071` comment `// #421 — these arms reach no `fastqc::run` call site.`; four bails at 1075/1081/1087/1093. The enclosing `if self.fastqc_requested()` block actually runs 1073-1110 (a fifth check, paired-single-file, at 1098-1109) — the plan's step 3 already says "the whole `fastqc_requested` block", so the citation is a deliberate narrowing, not an error |
| `cli.rs:2309-2350` = that block's four tests | **VERIFIED** | 2309 is the section comment, tests at 2314/2324/2334/2344; the fourth test's closing brace is 2351. Off by one line, immaterial |
| `cli.rs:2466` = `hardtrim5_range_error_precedes_the_fastqc_refusal` | **VERIFIED** | exact fn name at 2466, asserts `err.contains("between 1 and 999")` |
| `cli.rs:833-850` = `--clump_only` family whose wording was selected | **VERIFIED** | clip family 833-851; `--library` bail is `cli.rs:850`: `"--clump_only does not clip; --library is not compatible"` |
| `cli.rs:1009` = #386 `--hardtrim5`+`--hardtrim3` bail | **VERIFIED** | `1009: if self.hardtrim5.is_some() && self.hardtrim3.is_some()`, bail 1010-1014 |
| `cli.rs:1016-1024` = the two range checks | **VERIFIED** | 1016-1020 hardtrim5, 1021-1025 hardtrim3, both `must be between 1 and 999` |
| `cli.rs:1026-1031` = `validate_paired_input` for clock/implicon | **VERIFIED** | exact |
| `cli.rs:891-901` = `--clump_only` specialty mutual exclusions | **VERIFIED** | comment at 891, bails 893/896/899/902 |
| `main.rs:273-277` = the `--phred64` uniform-refusal precedent | **VERIFIED** | 273-277 is verbatim *"Rejected uniformly rather than per-mode… this does remove working invocations — but a mode-dependent rule … is worse to document, worse to test"* |
| `main.rs:393` = the widened `--rename`/uBAM guard | **VERIFIED** | `393: && (cli.effective_clips().any_set() \|\| cli.hardtrim5.is_some() \|\| cli.hardtrim3.is_some())` |
| `main.rs:1121` = stderr log call site, `cli.paired` in scope | **VERIFIED** | `1121: clips.flag_summary()`, inside `fn setup_trimming(cli: &Cli, …)` which opens at `main.rs:1046` |
| `report.rs:279` = text-report call site, `config.paired` in scope | **VERIFIED** | `279: writeln!(w, "Clipping in force: {}", library.flag_summary())?;`; `report::TrimConfig.paired` is `report.rs:172` |
| `library.rs:121-134` = `flag_summary` body | **VERIFIED** | exact. Note: the doc comment is 118-120 and contains **one** justification sentence, not the "two" the plan says — cosmetic |

### Substantive claims

| # | Plan claim | Verdict | Settled by |
|---|---|---|---|
| 2 | `flag_summary` has exactly two call sites | **VERIFIED** | `grep -rn "flag_summary" . --include="*.rs"` → 3 hits total: `library.rs:121` (def), `report.rs:279`, `main.rs:1121`. No test calls it — see finding V-1 |
| 3 | `both_reports_record_the_preset_and_its_expanded_values` is single-end, asserts all four, must be updated | **VERIFIED** | `tests/integration_library.rs:174-204`: args are `--dont_gzip --library pbat -o <dir> R1` (no `--paired`), and line 189 asserts `text.contains("--clip_R1 8 --clip_R2 8 --three_prime_clip_R1 8 --three_prime_clip_R2 8")`. The gate falsifies it; the suite goes red without step 10 |
| 4a | Insertion after the whole `fastqc_requested` block keeps ours last of the fallible checks | **VERIFIED** | the last `bail!` inside `validate()` is `cli.rs:1103`; from 1112 to the function's close at 1141 there are only `eprintln!` deprecation warnings. Inserting at 1111 makes ours the new last fallible check |
| 4b | `--hardtrim5 0 --library emseq` reports the range error | **VERIFIED** | range check 1016-1020 precedes 1111 |
| 4c | `--clock`/`--implicon` arity errors still win | **VERIFIED** | `validate_paired_input` 1026-1031 precedes 1111 |
| 4d | `--hardtrim5 + --hardtrim3` (#386) still wins | **VERIFIED** | bail at 1010 precedes 1111 |
| 4e | (extra, not claimed) missing input files still report first | **VERIFIED** | `check_restartable_input` loop at `cli.rs:1042-1044` precedes 1111 — consistent with the #421 comment's promise about "invalid input lists" |
| 4f | `--clump_only --hardtrim5 --library emseq` message | **VERIFIED, see F-6** | inside `if self.clump_only`, the `--library` bail (850) precedes the `--hardtrim5` mutual-exclusion (892), so the user gets `"--clump_only does not clip; --library is not compatible"`. Pre-existing family ordering; our block at 1111 is unreachable for clump-only runs |
| 5 | "Does not clip" is accurate for all four specialty modes | **REFUTED for `--clock` and `--implicon`** | `specialty.rs:241-243`: *"then clips: R1 at position 13 …, R2 at position 15"*; `specialty.rs:352`: *"then clips R2 by `umi_length` bases"*; `cli.rs:416`: *"then clip R2 by N bases"*. See C-2 |
| 6a | `accel` is 10/15/10/10 | **VERIFIED** | `library.rs:71` |
| 6b | `overrides` carry the preset's displaced values while resolved values are the user's | **VERIFIED** | `library.rs:165-176`: `pick` pushes `ClipOverride { flag, preset_value, user_value: v }` and returns `Some(v)`. With 1/2/3/4 overridden: resolved = 1/2/3/4, overrides = `[(--clip_R1,10,1), (--clip_R2,15,2), (--three_prime_clip_R1,10,3), (--three_prime_clip_R2,10,4)]`, in flag order. Independently assertable |
| 6c | Four distinct values kill *every* key permutation | **VERIFIED** | any pairwise swap of the four JSON key labels changes at least one of the four asserted values |
| 6d | T-1 is real (today's test survives the key swap) | **VERIFIED by inspection** | today's data is `emseq` + `clip_r1: 12` → resolved 12/10/10/10. Swapping the `clip_r2` / `three_prime_clip_r1` key labels leaves both asserted values at 10 → `report.rs:1453-1454` still pass |
| 6e | T-2 is real (today's test survives the flag-label swap) | **VERIFIED by inspection** | today's only override is `--clip_R1`; `report.rs:1456-1458` assert `overrides[0]` only, so swapping the `--clip_R2` / `--three_prime_clip_R1` literals at `library.rs:181/183` is invisible |
| 6f | Four overrides exercise the ≥2 comma path for the first time | **VERIFIED** | hand-written separator at `report.rs:894-896`; the only override-bearing test today has one entry, and the integration test asserts `overrides` is *empty* (`integration_library.rs:198-203`) |
| 7 | Two docs transcripts are single-end and get falsified | **VERIFIED** | `library-presets.md:48` is `trim_galore --library emseq --clip_R1 12 sample.fq.gz` (no `--paired`); its transcript at 51-54 shows all four on line 52 (that is the `main.rs:1121` log line). The text-report block at 62-66 shows all four on line 64. Both falsified. See F-3 for what the plan misses about the *adjacent* JSON block |
| 8 | Refusing breaks no existing test | **VERIFIED** | `grep -rn '"--library"' --include="*.rs" src/ tests/` → 5 hits in `tests/integration_library.rs`, none with a specialty mode; `grep -n library src/cli.rs` → the only exclusion test is `test_clump_only_rejects_library` (`cli.rs:1723`). Nothing combines `--library` with `--hardtrim5/3`, `--clock` or `--implicon` |
| 9a | Specialty modes write no trimming report | **VERIFIED** | `grep -n "report" src/specialty.rs` → **zero hits**. `main.rs:497` is `return Ok(())` at the end of the hardtrim5 arm, so the specialty arms never reach `setup_trimming` |
| 9b | Exactly three `report::TrimConfig` construction sites | **VERIFIED** | `main.rs:1449` (`run_single_file`, `paired: false`), `main.rs:2239` (`run_ubam_output_single`, `paired: false`), `main.rs:2618` (`write_paired_reports`, `paired: true`); plus one test helper at `report.rs:1345`. All three call `write_report_header` (1476 / 2266 / 2645), so the gate reaches line 279 on all of them, with the correct `paired` value hardcoded per arm |
| 9c | `validate()` precedes `ensure_output_dir`, so a refused run creates nothing | **VERIFIED** | `main.rs:237 cli.validate()?;` vs `main.rs:424 naming::ensure_output_dir(output_dir)?;` |
| 10 | Step 8: `test_text_report_records_preset_name_and_expanded_values` already sets `paired = true` | **VERIFIED** | `report.rs:1378 config.paired = true;` — stays green unchanged |
| 11 | (extra) T-3's premise: `main.rs:393` survives reversion green today | **VERIFIED** | `rename_into_ubam_refused_for_fastq_input` (`tests/integration_ubam_out.rs:1006-1039`) passes an explicit `--clip_R1 3`, so reverting 393 to the explicit `cli.clip_*` list keeps it green. `rename_into_ubam_refused_for_hardtrim_fastq_input` (1139) covers the hardtrim disjunct. Nothing covers the `--library`-only disjunct |
| 12 | (extra) assumption 6: does `--implicon 8` parse? | **REFUTED — `--implicon=8` is required** | `cli.rs:419-421` is `#[clap(long = "implicon", …, num_args = 0..=1, require_equals = true)]`. Run at `145432e`: `--implicon 8 --paired R1 R2` → `Error: Paired-end mode requires an even number of input files (R1/R2 pairs), got 3` — the `8` binds as a positional. `--implicon=8 --paired --library pbat …` starts running normally. See C-1 |
| 13 | (extra) is `--library` + specialty currently accepted and ignored? | **VERIFIED** | the `--implicon=8 --paired --library pbat` run above exits the validation stage and begins processing, confirming the bug the plan fixes |
| 14 | The suite is 751 tests at `145432e` | **UNCHECKED** | taken from the brief. I did not re-run `cargo test`; the pre-built binary proves the tree compiles at this commit |
| 15 | #441 preserves `--rrbs` byte-identity | **UNCHECKED** | out of this review's scope (that was the code-review round). The plan's V13 re-checks it |

---

## Findings

### Critical

Nothing blocks implementation. The four items below are Important; I have deliberately kept
"Critical" empty rather than inflating it.

### Important

**C-1 — Step 4, test 4 will not do what it says. `--implicon 8` must be `--implicon=8`.**

`cli.rs:419-421` sets `require_equals = true` with `num_args = 0..=1`, so a space-separated
value is not consumed by the flag. Measured at `145432e`:

```
$ ./target/release/trim_galore --implicon 8 --paired -o … R1 R2 2>&1 | head -5
Error: Paired-end mode requires an even number of input files (R1/R2 pairs), got 3
```

The plan's test asserts `err.contains("--library")`, so it fails **loudly** (not silently) —
but it fails, and the implementer then debugs the wrong thing. Assumption 6 can be closed now
rather than deferred: use `--implicon=8`, or copy the existing
`fastqc_plus_implicon_rejected` (`cli.rs:2345`), which passes bare `--implicon` and relies on
`default_missing_value = "8"`. The bare form is the better template — it matches the
neighbouring test it is mirroring, and it exercises the `is_some()` path the guard reads.

**C-2 — `"--clock does not clip"` and `"--implicon does not clip"` are false statements, and
the repo says so in its own doc comments.**

- `specialty.rs:241-243` — *"then clips: R1 at position 13 (removes UMI + CAGT + A), R2 at
  position 15 (removes UMI + CAGT + A + 2bp end-repair bias)"*
- `specialty.rs:352` — *"then clips R2 by `umi_length` bases"*
- `cli.rs:416` — *"then clip R2 by N bases"*

`--hardtrim5/3` are also not clip-free — they truncate to N bases. What is actually true is
that none of the four *honour the clip flags*, which is a different sentence. The
`--clump_only` wording is accurate for `--clump_only` precisely because that mode does nothing
to records at all (`cli.rs:878`: *"preserves record contents byte-identically"*); the property
does not transfer.

This does not require reopening Felix's decision 3 (the `--clump_only` wording *family*). The
same file already contains the correct phrasing for exactly these four modes, in the #369 `-a2`
block at `cli.rs:1050-1055`:

```rust
Some("--hardtrim5/--hardtrim3 perform no adapter trimming")
Some("--clock performs no adapter trimming")
Some("--implicon performs no adapter trimming")
```

The minimal-diff fix keeps the family shape and the terse register:

```
"--hardtrim5 performs no end clipping; --library is not compatible"
"--clock performs no end clipping; --library is not compatible"
```

If Felix prefers the literal `does not clip` cadence, `"--clock does not apply the clip flags;
--library is not compatible"` is the other option. Either way the assertion in the four tests
(`contains("--library") && contains("--clock")`) is unaffected, so this is a string change only.

**C-3 — The stderr log line, which the plan itself calls "the load-bearing detail", has no
automated coverage. V7 is a manual comparison.**

The plan's §Context is right that gating only `report.rs:279` would desync the log from the
report. But nothing in V1-V15 would *catch* that mistake: V7 says "Real SE and PE `--library
pbat` runs; compare the stderr line against the report line", which is an eyeball step. The one
existing test that reads stderr,
`explicit_clip_flag_overrides_the_preset_and_is_logged`
(`tests/integration_library.rs:130-169`), asserts only `err.contains("--clip_R1 12")` and
`err.contains("preset value 10")` — both satisfied by the **ungated override line**, so it stays
green whether or not `main.rs:1121` was changed. I checked this specifically: it is a near-miss,
not a catch.

Fix: add one assertion to that existing single-end test (or to a new one):

```rust
// the log line and the text report must agree about R2 on a single-end run
assert!(!err.contains("selected: --clip_R1 12 --clip_R2"), "log line still lists R2:\n{err}");
```

or, more robustly, assert the `Library preset '…' selected:` line's tail equals the report's
`Clipping in force:` tail. Without this, half of the plan's own headline catch is revertible in
silence.

**C-4 — #441's CHANGELOG entry becomes factually false, in the same `Unreleased` section as
ours.**

`CHANGELOG.md` (the #441 entry, currently at the top of `### Unreleased`) says:

> Both trimming reports record the preset name and all four values in force
> (`parameters.library` in JSON, `null` when no preset was used)

After the SE gate the text report records **two** values on a single-end run. Both entries ship
in the same release, so the changelog states a thing and then removes it two bullets later.

The plan excludes "trimming #441's CHANGELOG entry" from scope — correctly, that is editorial.
This is not trimming; it is the same category the plan explicitly admits for docs in §7
("Docs changes *required for factual correctness by our behaviour change* are in scope"). The
fix is a phrase: *"record the preset name and the values in force"*, dropping "all four". One
edit, no length change, no register change.

Two softer instances, both Optional rather than Important because they remain arguably true:
`library.rs:16-18` ("both trimming reports record the expanded values, so an old report stays
self-describing") and the `flag_summary` doc comment the plan is replacing anyway.

**C-5 — Step 6 silently drops the only test of the single-element `overrides` JSON array.**

`test_json_report_records_library_block` is the sole test with a non-empty `overrides` array,
and today it has exactly **one** entry (`report.rs:1456-1458`). Step 6 replaces its data with
four overrides. Afterwards:

- 0 overrides — covered (`integration_library.rs:198-203`)
- 4 overrides — covered (new unit test + new step-14 integration test)
- **1 override — covered nowhere**

The one-element case is the branch where `report.rs:894` (`if i + 1 < len`) must *not* emit a
comma. Trailing-comma bugs in hand-written JSON are exactly the class this test existed to
catch. Keep the one-override assertions: either leave the emseq/`clip_r1: 12` test in place and
add the accel/four-override test beside it, or add the JSON `overrides` length-1 assertion to
`test_text_report_names_both_values_for_an_override`'s neighbourhood. The plan's stated reason
for folding T-1 and T-2 into one test ("both need multiple overrides") does not require
*deleting* the single-override case.

### Optional

**F-1 — V5's first command matches zero tests, so it cannot fail.**

V5 reads `cargo test flag_summary`. There is no test whose name or body contains
`flag_summary` — the grep proving the two-call-site claim also proves this (3 hits, all
production code). `cargo test flag_summary` reports `0 passed; 0 failed` and exits 0. Per the
repo's standing rule this is precisely a check that cannot fail. Replace it with the real test
names.

Related and worth doing: **`flag_summary` has no direct unit test at all**, and it is the
function being changed. Two asserts in `library.rs`'s existing `mod tests` —
`resolve(Some(Pbat), default()).flag_summary(false) == "--clip_R1 8 --three_prime_clip_R1 8"`
and `…flag_summary(true)` with all four — are the cheapest possible guard for B2 and would sit
next to `preset_table_matches_the_agreed_numbers` (`library.rs:202`), which is the same kind of
"the numbers are the feature" test.

**F-2 — The paired half of the gate has no mutation listed, so an over-eager gate is not
provably caught.**

Step 13 (`paired_report_lists_all_four_values`) exists exactly to stop a future gate from
dropping R2 from *paired* runs. V5 lists it under "pass"; no V-row mutates anything that would
make it fail. Add one: change the guard from `if read_2 && !paired` to `if read_2`, and require
that `paired_report_lists_all_four_values` **and**
`test_text_report_records_preset_name_and_expanded_values` (`report.rs:1376`, `paired = true`)
both go red. Cheap, and it turns step 13 from a hopeful test into a proven one.

**F-3 — The docs page will show the asymmetry adjacently, and step 16-17's one sentence does
not cover it.**

`library-presets.md` lines 60-79 present, for the *same* single-end invocation, the text report
(soon without `--clip_R2`) and then, seven lines later, the JSON block **with** `"clip_r2": 10`.
The plan's added sentence ("R2 values are listed for paired runs only") sits with the text
transcript and reads as a global statement, which the JSON immediately below contradicts. A
reader's next move is a bug report.

The sentence needs to name both: something like *"The text report lists R2 values for paired
runs only; the JSON records the resolved configuration in full and consumers read `paired`
themselves."* Same length, closes the loop. Note the JSON transcript itself (70-79) needs **no**
value change — I checked, it stays correct under the plan's design.

**F-4 — V15's "byte-identical" is unachievable as written.**

The two transcripts are excerpts: a real run emits many more stderr lines and a much longer
report header, and `sample.fq.gz` is not a fixture in `test_files/`. `diff` against the page
cannot pass. The achievable check is "every quoted line appears verbatim in the real output"
(the quoted lines are filename-free, so substitution of a real fixture is safe). Say that,
otherwise the implementer either burns time or quietly downgrades the check and records a green.

**F-5 — Step 11's "enforced invariant" is conditional, and the condition is worth stating.**

The agreement `json.parameters[k] == json.parameters.library[k]` holds because the top-level
keys come from `config.clip_*` (`main.rs:1485-1488`, `2273-2276`, `2660-2663`) while the library
block comes from `cli.effective_clips()` (`main.rs:1471`, `2261`, `2640`) — two different
derivations. They agree only because the one place they can diverge, the `--rrbs` directional
auto-clip at `main.rs:1136`, fires only when `clips.clip_r2.is_none()`, and a preset always
supplies `clip_r2` because `PresetClips` fields are plain `usize` (`library.rs:23-26`), not
`Option`. So the invariant is real and type-enforced today — but it rests on "every preset sets
all four", which the plan never says. One line in §8 keeps a future reader from making
`PresetClips` fields optional and quietly breaking `--rrbs --paired --library <kit>`.

**F-6 — One sentence missing from the ordering analysis (`--clump_only` + specialty +
`--library`).**

Traced: for `--clump_only --hardtrim5 --library emseq` the user gets
`"--clump_only does not clip; --library is not compatible"` (`cli.rs:850`), because within the
`if self.clump_only` block the clip-flag family precedes the specialty mutual-exclusions
(`cli.rs:892`). Our new block at 1111 is unreachable for clump-only runs. That is sensible and
unchanged by this plan; the plan's O-2 covers the wording question but not the reachability, and
a reader checking "does the new block conflict with `:850`?" would benefit from the one line.

**F-7 — The new refusal shadows the #408 `--rename` guard for one combination.**

`--rename --output-format ubam --library emseq --hardtrim5 20` on FASTQ input reports the
#408 rename refusal today (`main.rs:391-397`, reached because `effective_clips().any_set()` is
true); after the change it reports the `--library`/`--hardtrim5` refusal, because
`validate()` runs at `main.rs:237` and the guard at 393. No test asserts that combination
(checked), and the hardtrim-only case is still covered by
`rename_into_ubam_refused_for_hardtrim_fastq_input` (`integration_ubam_out.rs:1139`), so
nothing breaks. Worth a line in §Integration's "currently-working invocations removed" so a
future bisect does not hunt for it.

**F-8 — `flag_summary` can now return an empty string; today unreachable, worth noting.**

With `paired = false` and only R2 values set, `parts` is empty and the report would print
`Clipping in force: ` with a dangling label. Unreachable today: `report.rs:275-276` prints the
line only when `library.preset` is `Some`, and every preset sets `clip_r1`. One line in §8
records why it is safe, so a future change that prints the line without a preset does not
produce the dangling label.

**F-9 — Test count.** The plan says "751 + ~12 new". Counting the plan's own steps: 5 in
`cli.rs`, 1 new in `report.rs` (2 modified), 2 new in `integration_library.rs` (2 modified), 1
in `integration_ubam_out.rs` → **9** new. Immaterial, and V14 correctly says not to gate on a
number.

---

## Assumptions

### Stated assumptions — assessment

| # | Plan's assumption | Assessment |
|---|---|---|
| 1 | #441 merges first, unmodified | Sound. I verified every anchor against `145432e`, which is #441's head; if #441 is revised, re-check steps 3, 10, 16, 17 |
| 2 | "Does not clip" is accurate for all four | **Wrong** — see C-2. The evidence the plan cites (`--clip_R1` is inert by md5) proves "does not honour the clip flags", which is a different proposition from "does not clip" |
| 3 | Text and JSON may differ on single-end | Defensible, and stronger than the plan argues — see the O-1 discussion below |
| 4 | Override log lines stay ungated | Sound, and the distinction is the right one. It also happens to be what keeps two existing tests green (`report.rs:1398`, `integration_library.rs:130`) — worth noting as a bonus, not a reason |
| 5 | 1/2/3/4 acceptable as test data | Sound |
| 6 | `--implicon` may need `--implicon=8` | **Now settled: it does.** See C-1 |
| 7 | Docs limited to falsified transcripts + one bullet | Sound in principle; incomplete in practice — F-3 (the adjacent JSON block) and C-4 (the CHANGELOG clause) are the same category and are missing |

### Assumptions the plan leaves implicit

- **A8 — Every `report::TrimConfig` site already knows its own paired-ness.** True (verified
  9b: three sites, `paired` hardcoded per arm), but the plan names only the two *`flag_summary`*
  call sites. The report call site is one function serving three configs, and the gate's
  correctness depends on all three being right. Now checked; worth recording so nobody re-checks.
- **A9 — Presets always supply all four values.** The load-bearing reason F-5's invariant holds.
- **A10 — `paired` is never carried on `ResolvedClips`, so the two consumers must independently
  hold the correct notion of paired-ness.** True today (`cli.paired` in `setup_trimming`,
  `config.paired` in the report), and the reason C-3's missing assertion matters: nothing
  mechanically ties them together.
- **A11 — The new block reads `self.library`, not `effective_clips()`.** Correct choice (the
  refusal is about the *flag*, not the resolved values), but unstated. If it read
  `effective_clips().any_set()` it would also refuse a bare `--clip_R1` under `--hardtrim5`,
  which today is an accepted no-op and is not in scope.

### The text/JSON asymmetry (O-1) — argued

**It is defensible, and the plan under-argues it.** The strongest reason is one the plan does
not give: **gating the JSON would contradict step 11.** Step 11 adds an enforced invariant that
`parameters.clip_r2 == parameters.library.clip_r2`. The top-level `clip_r2` key has echoed
inapplicable single-end values since before #441 and cannot be changed without a schema break.
So gating only `library.clip_r2` makes the JSON internally inconsistent *and* breaks the
invariant this same PR is adding; gating both is a breaking schema change for every MultiQC
consumer. The asymmetry is therefore close to forced, not merely tolerated. That belongs in
§8.3 — it converts a judgement call a future reviewer can reopen into a constraint they cannot.

**The cost is real but is a docs cost, not a code cost.** The one place the asymmetry is
visible side by side is `library-presets.md:60-79`, and F-3 fixes that with one sentence. The
two artefacts also have genuinely different jobs: the text report is prose about what happened
(`Clipping in force:` is an assertion of application), the JSON is a resolved-configuration
record whose consumer can read `paired` itself. That is the same split the plan already draws
correctly for the override lines.

**Is it worse than the inconsistency being fixed?** No, and the asymmetry is between two
*different* artefacts read by different audiences, while the bug being fixed is a false
statement inside one artefact read by a human. Fixing a false claim and accepting a documented
difference in scope between two formats is the right trade.

---

## Validation sufficiency

Per-finding regression coverage — would the fix be *silently* revertible?

| Finding | Regression lever | Verdict |
|---|---|---|
| Specialty refusal | V2 (delete the block → 4 red), V3 (reorder → ordering test red) | **Covered** |
| SE text line | V6 (revert `flag_summary` → unit + integration red) | **Covered** |
| SE **log** line (`main.rs:1121`) | none — V7 is manual | **NOT covered — C-3** |
| Paired half of the gate | none listed | **NOT covered — F-2** |
| T-1 (key permutation) | V8 | **Covered**; I independently confirmed today's test survives the swap |
| T-2 (flag labels) | V9 | **Covered**; likewise confirmed |
| T-3 (`main.rs:393`) | V10 | **Covered**; premise verified against `integration_ubam_out.rs:1006` |
| F-5a (present-vs-null) | V11 | **Covered**; `serde_json` really does return `Value::Null` for a missing key, so `is_null()` alone is blind — the `contains_key` assertion is the right lever |
| F-5b (top-level/library agreement) | V12 | **Covered** |
| Single-element `overrides` array | lost by step 6 | **NOT covered — C-5** |

Other validation notes:

- **V5's first command cannot fail** (F-1).
- **V15 cannot pass** (F-4).
- **V4 is good and specific** — "run in an empty dir, exit≠0, `ls` empty". The existing
  `rename_into_ubam_refused_for_fastq_input` shows a better mechanised form:
  `assert_eq!(dir_listing(&dir), vec!["sp.fastq"])` — the whole listing, not one filename. Step
  15 should copy that shape (`integration_ubam_out.rs:1032-1038`), and V4 could then be a test
  rather than a manual run.
- **V13's hazard note is exactly right** and worth keeping verbatim: `-o .` from a scratch
  directory, and count the files compared.
- **V13's scope is thin for what it claims.** "Every FASTQ and `.txt` identical" — the `.txt`
  trimming report is where the change lives, so a `--rrbs` **paired** case is the one that must
  be in the matrix (paired keeps all four, so the report text must be byte-identical there). A
  `--rrbs` single-end case with no `--library` is also unchanged (no preset → no
  `Clipping in force:` line at all, per `report.rs:275-276`). Naming those two shapes explicitly
  would make V13 falsifiable rather than a gesture.

---

## Alternatives

**Store `paired` inside `ResolvedClips` at construction instead of threading it through
`flag_summary(paired)`.** I would not. `ResolvedClips` is produced by `Cli::effective_clips()`
(`cli.rs:1162`), which has five production call sites (`main.rs:393`, `1116`, `1471`, `2261`,
`2640`) plus two in tests. Baking `paired` in would give the codebase *two* records of
paired-ness (`ResolvedClips.paired` and `report::TrimConfig.paired`) that can disagree — the
same class of desync the plan's "load-bearing detail" is guarding against, moved rather than
removed. It would also make the struct mean two things: `ResolvedClips` is currently a pure
value ("the four numbers plus their provenance"), and `main.rs:393` consumes it with no interest
in paired-ness at all. The threaded `bool` keeps one source of truth per artefact and costs two
one-word edits. **The plan's choice is correct.**

**Have `flag_summary` return the parts and let callers filter.** Worse for two callers: it moves
identical filtering logic into both, which is how the log and the report drift apart. Reject.

**A cheaper variant of the same idea, worth naming:** keep `flag_summary(paired)` *and* add the
two-line `library.rs` unit test from F-1. That gets the "the function itself is tested" benefit
without moving `paired` anywhere.

**Message wording:** reuse the #369 `-a2` phrasing precedent (C-2) rather than inventing new
wording. It is in the same file, forty lines above the insertion point, and already covers
exactly these four modes.

---

## Efficiency

Agreed with the plan: nothing here is on a hot path, and it says so honestly.
`flag_summary` runs once per report; the guard is four `Option::is_some()` calls once per
invocation. Two pedantic notes, neither worth acting on: the plan says "the array becomes a
3-tuple" (it becomes an array *of* 3-tuples), and `parts` could take `Vec::with_capacity(4)`.
Moving on.

---

## Conventions (CLAUDE.md)

- **`// The specialty arms return before setup_trimming, so no preset reaches a clip site.`** —
  one line, current state, no history, no justification. **Complies.** It states a mechanism
  rather than narrating the next line, which is the right side of that rule.
- **The `flag_summary` doc comment** ("R2 values are listed for paired runs only") — two lines,
  current state, and dropping the "so a run stays reproducible…" sentence is the correct call:
  that is rationale, and rationale belongs in the commit message. **Complies.** Minor: the plan
  says "the two justification sentences"; there is one (`library.rs:119-120`).
- **The section header** `// ── #440: --library refused where no clip site is reachable ──`
  matches the existing `// ── #421: …` shape at `cli.rs:2309`. **Complies.**
- **The CHANGELOG entry** — one sentence, bold lead, issue link. Matches the file's actual
  format (`CHANGELOG.md` links `[#440](…/issues/440)`, i.e. the issue not the PR, so the plan's
  link form is locally correct). **Complies with the one-sentence rule.** "no longer lists" is
  historical phrasing, which the rule bans in *source comments* — a changelog is inherently
  about change, so it is right here.
- One convention item the plan does not raise: **the plan's own framing line** ("a deliberate
  contrast with #441's 19-line entry") is plan prose, not source, so it is fine where it is —
  but it must not migrate into the commit message as a jab at a first-time contributor's entry.

---

## Action items

### Critical

None.

### Important

1. **C-2** — Change the four refusal messages. `--clock` and `--implicon` demonstrably clip
   (`specialty.rs:241-243`, `:352`). Use the in-file precedent at `cli.rs:1050-1055`:
   `"--clock performs no end clipping; --library is not compatible"`. String-only change; the
   tests' assertions are unaffected. Update assumption 2 to say what the md5 evidence actually
   proves ("does not honour the clip flags").
2. **C-1** — Fix step 4's fourth test: `--implicon=8`, or bare `--implicon` to mirror
   `cli.rs:2345`. Close assumption 6 as settled rather than "to be confirmed".
3. **C-3** — Add an automated assertion that the **stderr** log line is gated. Extend
   `explicit_clip_flag_overrides_the_preset_and_is_logged`
   (`tests/integration_library.rs:130`), which reads stderr already and today would stay green
   whether or not `main.rs:1121` was changed. Without this, half the plan's headline catch is
   silently revertible.
4. **C-4** — Add one factual edit to #441's CHANGELOG clause: *"record the preset name and the
   values in force"*, dropping "all four". Same scope rule the plan applies to docs in §7.
5. **C-5** — Do not delete the single-override `overrides` JSON case. Keep the emseq test and
   add the accel/four-override test beside it (or assert the length-1 array elsewhere). It is
   the only guard on `report.rs:894`'s no-trailing-comma branch.

### Optional

6. **F-1** — Replace V5's `cargo test flag_summary` (matches nothing, cannot fail) with the real
   test names, and add two direct `flag_summary(false)` / `flag_summary(true)` asserts in
   `library.rs`'s `mod tests` — the function being changed currently has no unit test.
7. **F-2** — Add a mutation for the paired half: `if read_2` without `&& !paired` must turn
   `paired_report_lists_all_four_values` and `report.rs:1376` red.
8. **F-3** — Make the docs sentence cover both formats, since `library-presets.md:60-79` shows
   the text and JSON transcripts for the same single-end run seven lines apart.
9. **F-4** — Downgrade V15 from "byte-identical" to "every quoted line appears verbatim".
10. **F-5** — Record in §8 that step 11's invariant holds because `PresetClips` fields are plain
    `usize`, which is what stops the `--rrbs` auto-clip at `main.rs:1136` from diverging them.
11. **F-6/F-7/F-8** — Three one-line additions: the `--clump_only --hardtrim5 --library`
    resolution (`cli.rs:850` wins); the #408 guard shadowing under
    `--rename --output-format ubam --library --hardtrim5`; and why an empty `flag_summary` tail
    is unreachable.
12. **V13** — Name the `--rrbs` shapes: a **paired** `--rrbs` case is the one that must prove the
    `.txt` report is byte-identical, since paired keeps all four values.
13. **Step 15** — Copy `integration_ubam_out.rs:1032-1038`'s whole-listing assertion
    (`assert_eq!(dir_listing(&dir), vec!["sp.fastq"])`) rather than a per-file check, and
    consider promoting V4 from a manual run to a test in the same shape.
14. **§8.3** — Add the "gating the JSON would contradict step 11" argument. It is the strongest
    defence of O-1 and the plan omits it.

---

## What I verified as correct (so absence of a finding means absence of a defect, not of a check)

Positively checked and found correct: all thirteen line-number anchors; the two-call-site claim
(and that no test calls `flag_summary`); the step-10 single-end-test claim; the placement's
ordering consequences including the four precedence cases the plan names plus the input-file
check it does not; the `--clump_only` interaction; the four-distinct-values lever including the
`overrides`/resolved independence and the ≥2 comma path; T-1's and T-2's "currently green under
mutation" premises by inspection; T-3's premise against the real test; both docs transcripts;
that no existing test combines `--library` with a specialty mode; that `specialty.rs` writes no
report; that there are exactly three `report::TrimConfig` production sites with correct
`paired` values reaching line 279; and that `validate()` precedes `ensure_output_dir`.

Not checked: the 751-test baseline (taken from the brief; the tree does compile at `145432e`),
and #441's `--rrbs` byte-identity (settled in the prior code-review round, re-checked by V13).

**Tree state on completion:**

```
$ cd /private/tmp/claude-501/pr441-B && git status --short
?? mutM.log
?? mutcheck.sh
?? testsummary.log
?? verify1.log
```

No source was mutated during this review. The two binary invocations used the pre-existing
`target/release/trim_galore` and an out-of-tree scratch `-o` directory; `git status --short
test_files/` is empty, confirming no fixture was touched.
