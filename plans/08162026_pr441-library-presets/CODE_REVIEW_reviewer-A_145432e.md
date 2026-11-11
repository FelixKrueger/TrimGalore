# Code review — PR #441 `feat(cli): --library presets expand into end clipping (#440)`

**Reviewer A** · commit `145432e179ae614814dd9754846b76c77137127b` · base `2ab0918`
Worktree: `/private/tmp/claude-501/pr441-A`

```
$ cd /private/tmp/claude-501/pr441-A && git rev-parse HEAD
145432e179ae614814dd9754846b76c77137127b
```

Nothing in this report was fixed in the PR. Every source mutation I made to test a hypothesis was
reverted; the worktree ends at a clean `git status` on `145432e`, and `test_files/` is byte-unchanged
(checksum `806b111292184f6a50f8eb5fb6ef2a59` before and after).

---

## Summary and verdict

**Mergeable with changes.** The design is implemented as agreed. No spec violation touches the five
preset numbers, the aliases, the membership criterion, the override semantics, or the
four-values-in-both-reports requirement — I verified each against the built binary, not against the
PR prose. `--rrbs` is untouched: a run with no `--library` produces byte-identical FASTQ output *and*
byte-identical `.txt` reports against a binary built from the base commit, with the JSON differing
only by the new `"library": null` key.

Three things want a maintainer decision before merge:

1. **F-1 (High)** — `--library` is accepted and silently ignored on all four specialty modes
   (`--hardtrim5`, `--hardtrim3`, `--clock`, `--implicon`). Not a regression (the four clip flags it
   stands for are inert there today, which I proved), but #421 in this repo's own CHANGELOG
   established that silently-ignored flags on the specialty modes are a defect worth a breaking
   change, and a preset turns one dropped flag into four.
2. **F-2 (Medium)** — on a **single-end** run the new text-report line reads
   `Clipping in force: --clip_R1 8 --clip_R2 8 --three_prime_clip_R1 8 --three_prime_clip_R2 8`.
   Two of those four are not in force; I proved the R2 values are inert by byte comparison. "In
   force" is a claim the pre-PR text report never made, and every SE preset run makes it. Since the
   report is what licenses changing preset values between releases, this weakens the property the
   whole versioning decision rests on.
3. **F-5 (Medium)** — two of the new report tests are insensitive to the defects they look like they
   guard. Deleting the JSON `library` key entirely leaves both "no preset → null" tests green, and an
   outright numeric disagreement between `parameters.library.three_prime_clip_r2` and
   `parameters.three_prime_clip_r2` goes undetected by the whole 751-test suite. Demonstrated with
   real mutated output below.

The rest is code quality: CHANGELOG entry length, three aliases that exist in code but in no
user-facing text, case-sensitive preset names, a hand-written `"library": null` where `json_null()`
exists, and comment verbosity against CLAUDE.md's stated rule.

The most important positive result: the PR's **one** differential integration test
(`preset_output_is_identical_to_the_flags_it_stands_for`) is the only test in the whole suite that
catches the silent-no-op class, and it does catch it. See "Mutation testing" — that test is
load-bearing and should not be weakened.

---

## Verification environment and commands

```
$ cargo build --release                                    # exit 0
$ cargo test                                               # 751 passed, 0 failed
$ cargo fmt --all -- --check                               # exit 0, no output
$ cargo clippy --all-targets --release -- -D warnings      # exit 0
```

Full-suite breakdown (17 test binaries; the brief's expected figure was 729, the actual on this
commit/platform is **751**, all green):

```
471 (lib) + 15 + 13 + 16 + 11 + 6 + 1 + 8 + 57 + 14 + 21 + 3 + 57 + 19 + 39 = 751 passed; 0 failed
```

**Clippy was proven able to fail before I trusted the clean pass.** A clean clippy run here prints no
`Checking trim-galore` line, which is the fingerprint-hit pattern the project notes warn about, so I
inserted a deliberate `needless_range_loop` + `useless_vec` into `src/library.rs`:

```
error: the loop variable `i` is only used to index `v`
error: useless use of `vec!`
error: could not compile `trim-galore` (lib) due to 2 previous errors
```

The probe was then reverted and the clean run reproduced with `Compiling trim-galore v2.3.0`
in the output. The clean result is real.

---

## Spec clauses verified SATISFIED

Absence of a finding for a row below means the check was run and passed.

| Clause | How verified | Result |
|---|---|---|
| Shape: one `--library` enum, not one boolean per kit | `src/library.rs:35-53`, `src/cli.rs:143-144` | SATISFIED — single `Option<LibraryPreset>` |
| Exactly five values | `src/library.rs:35-53` — `EmSeq`, `Accel`, `Zymo`, `ScBs`, `Pbat` | SATISFIED |
| No `rrbs` / `clock` / `implicon` / NuGEN value | same | SATISFIED |
| Exact numbers: emseq 10/10/10/10, accel 10/15/10/10, zymo 10/10/10/10, scbs 6/6/6/6, pbat 8/8/8/8 | `src/library.rs:69-75` read directly, **and** all five confirmed from live binary log lines | SATISFIED — all five rows match |
| Aliases `swift`, `xgen` → accel; `single_cell` → scbs | `src/library.rs:37-52`; binary run for each | SATISFIED (plus three extras, see F-4) |
| Scope: 5'/3' clipping only — no adapter set, no default `--length` | `PresetClips` has exactly four `usize` fields; `resolve()` touches nothing else; no `length`/`adapter` reference anywhere in `src/library.rs` | SATISFIED |
| Override WARNS and WINS, never errors | `src/library.rs:151-194`; `cli.validate()` accepts the combination; binary run | SATISFIED |
| Warning names BOTH the preset value and the user value | binary stderr: `--clip_R1 12 was given on the command line and overrides the emseq preset value 10` | SATISFIED |
| Preset **name** in the `.txt` report | `Library preset: pbat` observed in a real report | SATISFIED |
| All **four** expanded values in the `.txt` report | `Clipping in force: --clip_R1 8 --clip_R2 8 --three_prime_clip_R1 8 --three_prime_clip_R2 8` observed | SATISFIED (SE wording is F-2) |
| Preset **name** + all **four** values in the `.json` report | `parameters.library` observed with `preset` + four keys + `overrides` | SATISFIED |
| JSON stays valid: empty / single / multiple overrides | `serde_json` parse of real reports with 0, 1 and 4 overrides | SATISFIED |
| No user-supplied preset file / config discovery / env-var override | `grep -rn "env::var\|toml\|presets/\|read_to_string\|config_dir\|dirs::"` over changed files → 0 hits; `Cargo.toml` not in the diff | SATISFIED |
| `--rrbs` behaviour untouched for a run with no `--library` | base-vs-PR binary byte comparison over 6 invocation shapes | SATISFIED |
| Versioning documented with a CHANGELOG entry | `CHANGELOG.md` Unreleased §Changes | SATISFIED (length is F-3) |
| Two presets at once impossible | `--library emseq --library pbat` → `error: the argument '--library <PRESET>' cannot be used multiple times` (no silent last-wins) | SATISFIED |
| Unknown name refused at parse time | `--library nugen` → non-zero, error names `nugen` | SATISFIED |
| `cargo fmt` / `cargo clippy -D warnings` clean | commands above, with a failure probe | SATISFIED |

### `--rrbs` / no-preset byte-identity, against a base-commit build

I built the base commit into a separate tree (`git archive 2ab0918 | tar -x` into the scratchpad, then
`cargo build --release`) and compared every produced file:

```
=== se_plain       : compared 3 files, identical 2, differing 1   (json only)
=== se_clips       : compared 3 files, identical 2, differing 1   (json only)
=== pe_rrbs        : compared 6 files, identical 4, differing 2   (json only)
=== pe_rrbs_clipr2 : compared 6 files, identical 4, differing 2   (json only)
=== pe_rrbs_nondir : compared 6 files, identical 4, differing 2   (json only)
=== pe_plain       : compared 6 files, identical 4, differing 2   (json only)
```

Every FASTQ output and every `.txt` trimming report is byte-identical. `diff` on the JSON shows the
whole delta:

```
24a25
>     "library": null,
```

(the other hunk is `command_line`, which records the two different binary paths — my harness, not the
PR). **The CHANGELOG's claim that "the only difference anywhere is the new `"library": null` key in
the JSON report" is accurate.**

This also settles the review question about `src/main.rs:1136`, which changes the `--rrbs`
auto-`--clip_R2 2` guard from `cli.clip_r2.is_none()` to `clips.clip_r2.is_none()`. **The change is
correct and cannot alter a run with no `--library`:** `library::resolve(None, user)`
(`src/library.rs:152-161`) early-returns the user's four `Option`s verbatim, so
`clips.clip_r2 == cli.clip_r2` identically when no preset is given. The `--rrbs` (auto-clip fires),
`--rrbs --clip_R2 7` (user suppresses it) and `--rrbs --non_directional` (skipped) shapes above
exercise both sides of the guard and all came out byte-identical.

With a preset, `--paired --rrbs --library emseq` gives `clip_r2: 10` and no
`Setting the option '--clip_r2 2'` line — the suppression the spec and `docs/guide/flags.md`
describe. Worth noting as reassurance: every preset's `clip_r2` is ≥ 6, so suppressing the RRBS
auto-clip can never *reduce* R2 5' clipping below the 2 bp the end-repair mask wanted.

---

## Silent-no-op trace — every dispatch path

Each path run once with `--library pbat` (or `accel`), always with `-o <fresh dir>`; for the trimming
paths byte-compared against the four explicit flags **and** against no clipping (so the equivalence
result is non-vacuous).

| Dispatch path | Preset | Evidence |
|---|---|---|
| Single-end trim | **APPLIED** | log line + output ≠ no-clip run |
| Paired trim (two files) | **APPLIED** | `preset == four flags` byte-identical, `≠ no-clip` |
| `--paired` single interleaved uBAM → FASTQ | **APPLIED** | `_val_1.fq` md5: preset `de7ac194…` = flags `de7ac194…` ≠ none `bb256c8b…` |
| `--output-format ubam`, SE | **APPLIED** | `samtools view` record-body md5: preset `cc2a4fc4…` = flags `cc2a4fc4…` ≠ none `24fb8d5f…` |
| `--output-format ubam`, PE two files | **APPLIED** | body md5: preset `efb73f5a…` = flags `efb73f5a…` ≠ none `2be55b6c…` |
| `--output-format ubam`, single interleaved uBAM | **APPLIED** | preset log line present, 5 outputs |
| `--demux` | **APPLIED** | trimming runs first inside `run_single_file` (`main.rs:1511`) and demux splits the trimmed output; preset log present, 11 outputs |
| `--clumpify` | **APPLIED** | `_trimmed.fq.gz` md5 `c518a18a…` identical between `--library pbat` and the four flags |
| `--clump_only` | **REJECTED at CLI** | `Error: --clump_only does not clip; --library is not compatible` (`src/cli.rs:848-850`) |
| `--hardtrim5` | **SILENTLY IGNORED** | F-1 |
| `--hardtrim3` | **SILENTLY IGNORED** | F-1 |
| `--clock` | **SILENTLY IGNORED** | F-1 |
| `--implicon` | **SILENTLY IGNORED** | F-1 |

Note the uBAM whole-file md5s do differ between preset and flags — that is only the `@PG CL:` line
recording the two different command lines, which is expected per CLAUDE.md (`assert_ubam_eq` ignores
`@PG`). Comparing record bodies via `samtools view` isolates it.

Structurally this holds because after the PR **no non-test code reads `cli.clip_r1` / `cli.clip_r2` /
`cli.three_prime_clip_r*` directly** except `effective_clips()` itself and the `--clump_only`
rejection block:

```
$ grep -rn "cli\.clip_r\|cli\.three_prime_clip_r" src/
src/cli.rs:1625:  assert_eq!(cli.clip_r1, Some(5));
…                                        # all 8 hits are inside #[cfg(test)] modules
```

and all four `report::TrimConfig` construction sites (`main.rs:1449`, `:2239`, and `:2618` inside
`write_paired_reports`, which four callers funnel through at `:1801`, `:2008`, `:2429`, `:2548`) set
`library: Some(cli.effective_clips())`. That single-funnel design is the right shape for this defect
class and is the main reason the trace comes out as clean as it does.

---

## Mutation testing — proving the load-bearing tests can fail

Each mutation was applied to the worktree, the suite run, then the file restored from a pre-mutation
snapshot. None of this is in the PR.

| # | Mutation | Tests that went RED |
|---|---|---|
| M1 | `Accel => (10, 15, 10, 10)` → `(10, 14, 10, 10)` | `library::preset_table_matches_the_agreed_numbers`, `library::preset_supplies_all_four_values_when_user_gave_none`, `cli::test_library_preset_supplies_all_four_clip_values`, `report::test_text_report_records_preset_name_and_expanded_values`, integration `preset_output_is_identical_to_the_flags_it_stands_for` (**5**) |
| M2a | `resolve()` lets the **preset** win over an explicit flag (record kept) | `library::explicit_clip_flag_wins_over_the_preset`, `cli::test_explicit_clip_flag_overrides_the_preset`, `report::test_json_report_records_library_block`, integration `explicit_clip_flag_overrides_the_preset_and_is_logged` (**4**) |
| M2b | `resolve()` stops recording the override (user value still wins) | `library::an_override_records_both_values`, `cli::test_explicit_clip_flag_overrides_the_preset`, `report::test_json_report_records_library_block`, `report::test_text_report_names_both_values_for_an_override`, integration `explicit_clip_flag_overrides_the_preset_and_is_logged` (**5**) |
| M3 | text-report library block deleted from `write_report_header` | `report::test_text_report_records_preset_name_and_expanded_values`, `report::test_text_report_names_both_values_for_an_override`, integration `both_reports_record_the_preset_and_its_expanded_values` (**3**) |
| M4 | JSON `library` key deleted entirely (both match arms) | `report::test_json_report_records_library_block`, integration `both_reports_record_the_preset_and_its_expanded_values` (**2**) — **but see F-5a**, two tests that should have caught this stayed green |
| M5 | `TrimConfig.three_prime_clip_r2` reads `cli.*` instead of `clips.*` — the exact silent-no-op shape, on one of the four values | integration `preset_output_is_identical_to_the_flags_it_stands_for` **only** (**1**); all 471 lib tests stayed green |
| M6 | SE JSON top-level `three_prime_clip_r2` set to `config value + 100`, disagreeing with the `library` block | **nothing** — 751/751 green. See F-5b |

Every load-bearing test the brief named does fail under an appropriate mutation, so none of them is
vacuous. M5 is the important structural result: **the entire silent-no-op class rests on one
integration test.** The unit tests all exercise `resolve()` / `effective_clips()`, which are correct;
the defect M5 simulates lives in the wiring from `clips` into `TrimConfig`, which only the
byte-differential test can see.

Also worth recording: `preset_changes_the_output_relative_to_no_clipping` (the SE non-vacuity test)
stayed **green** under M5. It only asserts "output differs from a no-clip run", so it is blind to
three of the four values being dropped, as long as one still applies. That is the
"an assertion is only sensitive to what it counts" pattern — the test is not wrong, but it is not
protection.

---

## Spec violations

None found. F-1 and F-2 are behavioural defects that no stated clause forbids; F-4 is a small
addition beyond the agreed alias list. Recorded separately from code quality below because they land
differently.

---

## Findings — behaviour

### F-1 (High) — `--library` is accepted and silently ignored on all four specialty modes

`src/main.rs:468`, `:499`, `:530`, `:552`

**What the brief requires.** Priority 2: "A preset that is accepted on the command line but never
applied is the worst failure mode here, because the user gets untrimmed ends and a report that may
still claim clipping."

**What the code does.** The four specialty modes `return Ok(())` before `setup_trimming()` — the only
place `effective_clips()` reaches `TrimConfig` — is ever called. Nothing in `Cli::validate()` rejects
`--library` alongside them. The run exits 0, prints no warning, and specialty modes write no trimming
report, so there is no record either. There is not even the `Library preset '…' selected:` line,
because that too lives in `setup_trimming` (`main.rs:1117-1132`).

**How I verified it.** One input, one `--hardtrim5 30`, three flag sets:

```
$ trim_galore --hardtrim5 30 --library pbat -o D1 test_files/BS-seq_10K_R1.fastq.gz
$ trim_galore --hardtrim5 30 --clip_R1 8 --clip_R2 8 \
      --three_prime_clip_R1 8 --three_prime_clip_R2 8 -o D2 …
$ trim_galore --hardtrim5 30 -o D3 …
$ md5 -q D{1,2,3}/BS-seq_10K_R1.30bp_5prime.fq.gz
5de6f7398f4f57025457eacaca1db1fa
5de6f7398f4f57025457eacaca1db1fa
5de6f7398f4f57025457eacaca1db1fa
```

All three identical. No `Library preset` line on stderr for `--hardtrim5`, `--hardtrim3`, `--clock`
or `--implicon=8`.

**Weighing it.** This is **not a regression** — the middle run shows the four clip flags are already
inert under `--hardtrim5`, so `--library` merely inherits that. Two things argue for handling it in
this PR anyway:

- This repo's own CHANGELOG carries "**Breaking:** `--fastqc` is now refused on the four specialty
  modes … instead of being silently ignored" (#421). Same defect shape; the maintainer already chose
  refusal over silence once, as a breaking change.
- A preset compresses four flags into one token, so a copy-pasted `--library emseq --hardtrim5 30`
  drops four values rather than one, and there is no report in which to notice.

**Suggested fix.** One `anyhow::bail!` per mode in `Cli::validate()`, following the existing
`--clump_only` block at `src/cli.rs:833-850` or the `--passthrough` family at `:955-967`. If the
maintainer would rather not widen the PR, the minimum is a line in
`docs/src/content/docs/guide/library-presets.md` §"What a preset does *not* set", which currently says
nothing about the specialty modes ignoring `--library`.

---

### F-2 (Medium) — on single-end runs, "Clipping in force" names two values that are not in force

`src/report.rs:276-277` (the `Clipping in force:` line) · `src/library.rs:121-134` (`flag_summary`)

**What the code does.** `flag_summary()` renders whichever of the four `Option`s are `Some`, with no
knowledge of `paired`. Every preset sets all four, so every SE run prints all four.

**Observed, from a real report** (`--dont_gzip --library pbat`, single-end):

```
Trimming mode: single-end
…
Library preset: pbat
Clipping in force: --clip_R1 8 --clip_R2 8 --three_prime_clip_R1 8 --three_prime_clip_R2 8
```

**Proof the R2 halves are inert.** A preset run and a run with only the R1 clips give identical bytes:

```
$ md5 -q se_pbat/BS-seq_10K_R1_trimmed.fq       # --library pbat
aacae45a41dd1c48ea127de248754500
$ md5 -q se_r1only/BS-seq_10K_R1_trimmed.fq     # --clip_R1 8 --three_prime_clip_R1 8
aacae45a41dd1c48ea127de248754500
```

**Why it matters here specifically.** The spec makes the report load-bearing: it is what licenses
changing preset values between releases, because "an old report stays self-describing and a past run
[is] reproducible from its own record". A line asserting four numbers are *in force* when two were
never applied weakens exactly that property — someone reconstructing the run from the report
reproduces it only by luck (the two extra flags happen to be inert on SE input too).

The pre-PR JSON already carried `clip_r2` on SE runs when a user typed `--clip_R2`, so that half is
inherited. The **text** line is new: pre-PR the text report named no clip values at all, so
"in force" is a new claim by this PR.

**Suggested fix.** `TrimConfig.paired` is already in scope at `src/report.rs:270`. Either pass it into
`flag_summary(paired: bool)` and skip the two R2 entries when `!paired`, or keep the full expansion
and reword to `Preset expands to:`, which is true either way. The former is better — it keeps the
line reproducible exactly as typed.

---

## Findings — tests

### F-5a (Medium) — the "no preset → `library` is null" tests cannot tell null from absent

`src/report.rs:1447-1462` (`test_json_report_library_is_null_without_a_preset`) ·
`tests/integration_library.rs:208-215` (`a_run_without_a_preset_reports_no_library_block`)

Both assert `json["parameters"]["library"].is_null()`. `serde_json` returns `Value::Null` when you
index a **missing** key, so the assertion holds whether the key is present-with-null or absent
entirely. Mutation M4 deleted the `library` key from the JSON writer completely and **both tests
stayed green**.

The documented contract is the stronger one — `CHANGELOG.md` says "`parameters.library` in JSON,
`null` when no preset was used" and `docs/…/library-presets.md` repeats it — so a consumer keying on
the field's presence is not protected by these tests. `test_json_report_records_library_block` does
catch outright removal, so the gap is narrow, but a change that emitted the block only when a preset
exists (and omitted the key otherwise) would pass the entire suite while breaking the documented
shape.

**Suggested fix.** Assert on presence, not on null-ness — either
`assert!(json["parameters"].as_object().unwrap().contains_key("library"))` alongside the existing
check, or assert the raw text contains `"library": null`.

### F-5b (Medium) — nothing detects the `library` block disagreeing with the top-level clip values

`tests/integration_library.rs:193-203`

The brief's priority 3 asks whether the two JSON locations can disagree. **Today they cannot** — the
only mutation between `clips` and `config.clip_*` is the RRBS auto-`clip_r2`, which is gated on
`clips.clip_r2.is_none()`, and every preset sets all four, so the guard can never fire under a
preset. I confirmed agreement on the SE, PE, uBAM-PE and interleaved paths, and on
`--rrbs --library emseq` (both read 10).

But **no test enforces it.** `both_reports_record_the_preset_and_its_expanded_values` asserts
`json["parameters"]["clip_r1"]` against the top level and `library["clip_r1"]` /
`library["three_prime_clip_r2"]` inside the block — it never compares the other three top-level
values. Mutation M6 made the SE top-level `three_prime_clip_r2` disagree with the block by +100 and
the whole 751-test suite stayed green. The real report the mutated binary produced:

```
library block : {"preset": "pbat", "clip_r1": 8, "clip_r2": 8, "three_prime_clip_r1": 8, "three_prime_clip_r2": 8, "overrides": []}
top level     : {'clip_r1': 8, 'clip_r2': 8, 'three_prime_clip_r1': 8, 'three_prime_clip_r2': 108}
```

Note also the latent version of this: the "they cannot disagree" invariant holds *because all five
presets happen to set all four values*, not by construction. A future preset that left one value
unset would let the RRBS auto-clip fire and put `2` at the top level while the `library` block
printed `null` — silently, since nothing checks.

**Suggested fix.** One loop in that test asserting all four keys agree between
`json["parameters"]` and `json["parameters"]["library"]`. Cheap, and it converts the current
coincidence into an enforced invariant.

### F-5c (Low) — no differential test on the uBAM or interleaved paths

`tests/integration_library.rs`

The byte-differential test that carries the whole silent-no-op class (M5 above) covers only the
paired FASTQ path; `preset_changes_the_output_relative_to_no_clipping` covers SE non-vacuously but is
blind to three of four values (shown above). The uBAM-output paths and the interleaved-uBAM path have
no `--library` test at all. They work today — I verified all four by record-body md5 — but the
structural argument that one funnel protects them is only as good as `setup_trimming` staying the
single source, which is exactly what a future refactor could break.

**Suggested fix.** Optional, low cost: one `--output-format ubam` case in
`tests/integration_library.rs` comparing `samtools`-free record bytes, or simply reusing the existing
`assert_ubam_eq` helper the repo already has for `@PG`-insensitive comparison.

---

## Findings — code quality

### F-3 (Low) — CHANGELOG entry is 19 lines across two paragraphs

`CHANGELOG.md:8-25`

The global rule is "one sentence per entry, no exceptions", and this entry is a long way past it — two
paragraphs covering the preset table, the membership criterion, override semantics, the RRBS
interaction, the reporting contract and the no-preset guarantee.

In fairness, **it matches this file's own house style exactly**: the adjacent #434 entry immediately
below runs to two paragraphs and ~26 lines, #428 and #398 are similar. So the entry is consistent
with the file it lives in and inconsistent with the stated rule, and only the maintainer can say
which wins here. If it is to be cut, the reduction that loses least is: keep the first sentence plus
the preset list, and move the membership criterion, override semantics and reporting contract to
`docs/guide/library-presets.md`, which already says all three.

### F-4 (Low) — three aliases exist in code but in no user-facing text

`src/library.rs:37` (`em_seq`, `em-seq`), `:48` (`single-cell`)

The agreed alias list is `swift`/`xgen` for accel and `single_cell` for scbs. The PR also ships
`em_seq`, `em-seq` and `single-cell`, all of which work:

```
--library em_seq       -> Library preset 'emseq' selected: --clip_R1 10 …
--library em-seq       -> Library preset 'emseq' selected: --clip_R1 10 …
--library single-cell  -> Library preset 'scbs' selected: --clip_R1 6 …
```

None appears in `--help`, in the docs table, or in the CHANGELOG (clap hides `PossibleValue` aliases
from help by design). Harmless behaviourally — they cannot change any of the five rows — but they are
an accepted-input surface nobody has documented and nobody will remember to keep. Either document
them next to the others or drop them.

### F-6 (Low) — preset names are case-sensitive, unlike every clip flag alias

`src/cli.rs:143`

`--library EMSEQ` is refused (`error: invalid value 'EMSEQ' for '--library <PRESET>'`). The repo goes
out of its way to accept both `--clip_R1` and `--clip_r1`, and `--library XGEN` is the kind of thing a
user copying from a vendor PDF will type. One attribute (`ignore_case = true` on the arg) fixes it.
Not a spec matter — noting it because the rest of this CLI is forgiving about case.

### F-7 (Low) — `json_null()` exists but the `library: null` line is hand-written

`src/report.rs:900` writes `writeln!(w, "{}\"library\": null,", i2)?` where
`json_null(w, i2, "library", true)` (`src/report.rs:1204`) does exactly that and is what every other
optional field in this writer uses. Same for the `overrides` array, which is emitted on one line with
manual comma placement while every other field in the file is one-key-per-line. Both are *correct* —
I parsed real reports with 0, 1 and 4 overrides and all are valid JSON with right comma placement —
so this is consistency, not a bug.

### F-8 (Low) — comment verbosity against CLAUDE.md's stated rule

CLAUDE.md asks for one line stating the fact, and specifically "never justify the change in the
source". The clearest cases:

- `src/library.rs:1-18` — an 18-line module doc. The membership criterion itself is well sited (it is
  an invariant on the enum someone will edit, which CLAUDE.md explicitly endorses), but the
  *justifications* around it — "That is why `--rrbs` stays its own flag …", "That is safe because both
  trimming reports record the expanded values" — are commit-message material.
- `src/library.rs:30-33` — "One enum rather than one boolean per kit: a prep is a single choice, and an
  enum makes 'two kits at once' unrepresentable instead of needing a conflict rule for every pair."
  Pure design justification.
- `src/library.rs:145-150` — "An error instead of an override would push the user back to writing all
  four values by hand …" Same.
- Test doc comments in all three test files are discursive in the same way ("Non-vacuity: … Without
  this, the equivalence test above would still pass if presets did nothing at all"). CLAUDE.md's
  comment rule covers committed tests.

Two things in the PR's favour, both checked: **the "current state only, never reference the state you
are fixing" sub-rule is fully satisfied** — I found no "no longer", "previously", or historical
narrative anywhere in the new comments. And the surrounding repo is written in exactly this register
(`src/main.rs` is dense with multi-paragraph comments citing plan sections and PR numbers), so the PR
is stylistically consistent with the file it joins. As with F-3, the honest reading is "matches the
codebase, exceeds the rule".

### F-9 (Low) — the preset log line repeats once per input file, after adapter detection

`src/main.rs:1117-1132`

`setup_trimming` runs per input file (SE) and per pair (PE), so a 10-file SE run prints
`Library preset 'pbat' selected: …` ten times. It also prints *after* adapter auto-detection, which
can scan up to 1 M reads, so the first thing the user sees about their preset arrives well into the
run. Both are consistent with how the `--rrbs` messages behave a few lines below, so this is a nit,
not a defect. Hoisting the announcement to `main()` before dispatch would print it once and early —
and would incidentally give F-1 a natural home.

---

## Docs check

`docs/src/content/docs/guide/library-presets.md` (new), `guide/flags.md`, `astro.config.mjs`.

Register is fine — plain scientific prose, no AI-slop phrasing, no banned constructions. Every factual
claim on the page I could test, I tested, and all hold:

- the preset table matches `src/library.rs` and the spec, all five rows;
- the "`--library emseq` is exactly these four flags" equivalence — verified byte-identical;
- the override log-line transcript — matches real stderr verbatim;
- the text-report transcript and the JSON `parameters.library` shape — match real reports;
- "a preset that supplies `--clip_R2` also suppresses the `--rrbs` auto-clip" (`flags.md`) — verified:
  `--paired --rrbs --library emseq` gives `clip_r2: 10` and no `Setting the option '--clip_r2 2'`;
- the cross-link `/modes/rrbs/#when-not-to-use---rrbs` resolves (`docs/…/modes/rrbs.md:40`,
  `## When NOT to use \`--rrbs\``);
- the sidebar entry is registered in `astro.config.mjs`.

One gap, tied to F-1: §"What a preset does *not* set" lists no adapter, no `--length`, no mode
changes — but not that the four specialty modes ignore `--library` outright.

`--help` output is accurate and complete for the documented five values; it does not print the
numbers, which is reasonable given they live one link away.

---

## Hazards observed

- Every trimming invocation used `-o <fresh dir>`; `test_files/` was checksummed before and after the
  session and is byte-unchanged (`806b111292184f6a50f8eb5fb6ef2a59`).
- File counts were checked per case, and "produced no output" was treated as a failure — that is how
  the `--implicon 8` case was caught as a harness error (`--implicon` needs `require_equals`, so
  `--implicon=8`) rather than being scored as a pass.
- No `gh` command was run. Nothing was committed, pushed, approved or commented.
- All mutations reverted; `git status --short` in the worktree is empty at `145432e`.
