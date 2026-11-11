# CODE REVIEW — implementation, Reviewer A

**Target:** `3161b8a541c09d76d77e4e140f7b2638b8333d5c` on `fix/440-clip-r2-se-warning`, parent `fdade7e`
**Plan:** `plans/08172026_clip-r2-se-warning/PLAN.md` v2 (incl. §12)
**Worktree:** `/private/tmp/claude-501/crA440` (detached at `3161b8a` — verified)
**Reviewer:** A (independent; Reviewer B's report not read)

```
$ git rev-parse HEAD
3161b8a541c09d76d77e4e140f7b2638b8333d5c
$ git log --oneline -2
3161b8a fix(cli): warn when a clip flag cannot be honoured, and stop announcing it (#440)
fdade7e perf(clump): refit the --clumpify memory model against measured peak RSS (#439) (#446)
```

---

## Verdict — **READY WITH CHANGES**

**No defect in shipped behaviour.** I could not find a single reachable invocation where the printed
reason is false, where an inert flag goes unwarned, or where the two gated surfaces disagree — 41 live
invocations, enumerated below: all four specialty modes with and without `--paired`, single-file and
two-file, `--clump_only`, `--clumpify`, `--demux`, `--passthrough`, uBAM in, uBAM out, the
interleaved-uBAM paired shape, multi-pair (4 files) and multi-file single-end. The gate is correct on
both surfaces at every report site, including the two uBAM report sites §12's differential did not cover.

What needs changing is **test coverage of two of the code's own branches**. I mutated each and watched it
ship green over the full 780-test suite:

| | Mutation | Suite | User-visible consequence |
|---|---|---|---|
| **H-1** | drop the `all_four` narrowing (`cli.rs:1155-1163`) | **780 pass** | a plain `--clip_R1 3 <se>` run prints a sentence false in both halves |
| **M-1** | delete the `--three_prime_clip_R1` entry from the warning loop | **780 pass** | that flag silently ignored again under all four specialty modes |

Plus four documentation-accuracy items, all of the same shape: the new user-facing prose ("the trimming
report", "either surface") and one pre-existing source comment do not survive the JSON asymmetry that
Q-2 deliberately kept.

Everything §12 claims about itself that I re-ran, held: **780 tests**, `fmt` clean, **release** clippy
clean (with a probe proving the harness can fail), V6, V8a, V8b and V12 each red on exactly the tests
listed, V16's triple, and V18's byte-identity — which I reproduced independently at **43 files, 39
identical, 4 different, all four being the same single line**.

---

## Method — what was actually run

- `cargo build --tests`, then a **41-invocation mode matrix** against the debug binary (script:
  `scratch_a/matrix.sh`, outputs under `scratch_a/runs/`). Every invocation passed `-o <scratch dir>`.
- `cargo test --no-fail-fast` (full, 18 binaries) as the baseline for every mutation.
- Six mutations, each applied with a `python3` in-place edit asserting a unique match, then reverted with
  `git checkout --` and confirmed with `git status --short`.
- A clean **base-vs-new differential**: both sides built with a pinned `SOURCE_DATE_EPOCH`, invoked as
  `./tg` from their own directory with identical absolute input and `-o` paths, so `command_line`, the
  provenance line and the logged output paths are identical by construction.
- `cargo fmt --all -- --check`, `cargo clippy --all-targets -- -D warnings`, and
  `cargo clippy --all-targets --release -- -D warnings` — the last from a cleaned target dir.

**One harness correction, stated because it changed a result.** My first differential was contaminated:
after reverting the V6 mutation I re-ran the binary without rebuilding, so the "new" side carried the
`effective_clips()` mutation and reported spurious warnings and a differing `.bam`. I rebuilt, re-ran
everything affected (the V16 triple and the whole differential), and only the clean numbers are reported
above and below. The `.bam` difference disappeared: it was the `@PG CL:` line recording the two `-o`
paths, not a content change. The matrix, `surfaces.sh` and the multi-pair probes all predate any
mutation and were unaffected.

---

## Findings

### Critical
None.

### High

#### H-1 — the `all_four` narrowing is the load-bearing half of the ladder and no test covers it

`cli.rs:1153-1173` decides *which* flags are candidates by threading `all_four` through the array:

```rust
// src/cli.rs:1154-1165
for (flag, value) in [
    ("--clip_R1", if all_four { self.clip_r1 } else { None }),
    (
        "--three_prime_clip_R1",
        if all_four { self.three_prime_clip_r1 } else { None },
    ),
    ("--clip_R2", self.clip_r2),
    ("--three_prime_clip_R2", self.three_prime_clip_r2),
]
```

That `else { None }` is what keeps a single-end run from warning about a flag it honours perfectly well.
I removed it (kept `all_four` alive with `let _ = all_four;` so the compile is clean) and ran the whole
suite:

```
M4 applied: all_four narrowing removed (R1 flags now warn on single-end too)
M4 TOTAL passed=780 failed=0
```

Then the same binary, on the most ordinary invocation there is:

```
$ trim_galore --dont_gzip --clip_R1 3 -o <dir> test_files/BS-seq_10K_R1.fastq.gz
WARNING: --clip_R1 3 was given but is not used in this mode (it applies to Read 2 of a pair,
         and this is a single-end run). Ignoring.
```

Both halves are false — `--clip_R1` **is** used, and it does not apply to Read 2 — and the run then
clips R1 by 3 anyway. This is a worse sentence than the one v1 shipped, on a far more common invocation,
and the suite is silent about it.

The plan's V4 tests the *opposite* direction only (hard-code `all_four = false` → the specialty rows go
red). The reverse mutation has no row and no test. Note the structural reason it slipped: `read_2` is
now encoded in **three** places — `flag_summary`'s array (`library.rs:124-128`, pinned by its own
string-equality tests), `pick`'s four call sites (`library.rs:187-200`, pinned by the new oracle test) —
and **positionally, in this loop**, where the R1 entries are gated and the R2 entries are not. The
commit built an oracle for the second list and left the third one bare.

**Fix (one test):** a single-end row asserting silence with `--clip_R1` set, e.g. in
`tests/integration_inert_clip_flags.rs` beside `preset_alone_on_single_end_does_not_warn`:

```rust
/// A single-end run honours the Read 1 pair, so naming those flags must stay silent.
#[test]
fn single_end_read_1_clip_flags_do_not_warn() {
    let dir = tempdir("se_r1_quiet");
    let (ok, err) = run(&[
        "--dont_gzip", "--clip_R1", "3", "--three_prime_clip_R1", "4",
        "-o", dir.to_str().unwrap(), R1,
    ]);
    assert!(ok, "run failed:\n{err}");
    assert!(!err.contains(INERT), "a Read 1 clip a single-end run performs must not warn:\n{err}");
}
```

That single test kills H-1 and M-1's single-end half at once, and it is the natural twin of
`paired_clip_flags_do_not_warn` (`tests/integration_inert_clip_flags.rs:126`).

### Medium

#### M-1 — `--three_prime_clip_R1` is never exercised by the warning loop

The mode table (`tests/integration_inert_clip_flags.rs:51`) uses `--clip_R1` on the four specialty rows
and `--clip_R2` on the single-end row; `each_inert_flag_gets_its_own_line` (`:166`) uses the two R2
flags. Nothing anywhere passes `--three_prime_clip_R1`. Deleting its array entry outright:

```
M5 applied: --three_prime_clip_R1 dropped from the warning loop
M5 TOTAL passed=780 failed=0

$ trim_galore --dont_gzip --hardtrim5 30 --three_prime_clip_R1 4 -o <dir> …R1.fastq.gz
0     # warning lines — the flag is silently ignored again
```

One quarter of the feature ships uncovered. This is precisely Reviewer A's plan-review I-6 against v1
("half the feature is never exercised, so deleting its `if let` block would keep all tests green"),
re-appearing one flag over. **Fix:** use `--three_prime_clip_R1` on one of the four specialty rows —
the `hardtrim3` row is the natural home, since its reason string is identical to `hardtrim5`'s and the
row is otherwise a duplicate.

#### M-2 — the new `flags.md` bullet says "the trimming report" where two reports disagree

`docs/src/content/docs/guide/flags.md:15`:

> The trimming report lists only the clipping that applied, so a single-end report names Read 1 alone.

Trim Galore writes **two** trimming reports, and this is true of one of them. Measured on
`--dont_gzip --clip_R2 5 <se>`:

```
-- text report clip lines: [0]
mode = single-end | clip_r2 = 5 | clip_r1 = None | library = None
```

and with a preset, the JSON keeps an `overrides` entry naming `--clip_R2` (`json_ovr_entries=1` for the
single-end shape in my surface table below). The same CHANGELOG entry two sentences later uses "**Both**
trimming reports" for exactly this txt/JSON pair, so "the trimming report" is ambiguous in a file that
elsewhere distinguishes them, and false for the JSON. Second, smaller half: "names Read 1 alone" implies
it names Read 1 — on a bare `--clip_R2 5` single-end run the text report names **no** clipping at all.

**Fix:** "The **text** trimming report lists only the clipping that applied, so a single-end report names
Read 1 clipping only; the JSON report records what was requested."

Register otherwise checked and clean: plain scientific prose, no second person, no conversational
heading, same bullet shape as its four neighbours, and the `--library` bullet above it is untouched.

#### M-3 — "nothing on either surface" has no antecedent for a changelog reader, and the clause reads as preset-only

`CHANGELOG.md:20-22`:

> — and a clip flag those modes cannot honour, or a Read 2 clip flag on a single-end run, now warns
> instead of passing silently, with nothing on either surface announcing a Read 2 clip the run did not
> perform.

Two problems, both from internal vocabulary leaking outward:

1. **"either surface"** names neither. In the plan they are the stderr banner and the text report; to a
   reader holding a JSON report that still says `clip_r2: 5`, the sentence is simply wrong. Naming them
   costs three words: "…and neither the log nor the text report announces a Read 2 clip the run did not
   perform."
2. The clause sits inside the `--library <preset>` bullet, but the warning fires **with no preset at
   all** — `--clip_R2 5 <se>` and `--hardtrim5 30 --clip_R1 3` both warn. A reader scanning for
   clip-flag behaviour will not look inside a preset entry. Folding rather than adding a second bullet
   was the deliberate call (§11, B's O-4) and I am not reopening it; the fix is to make the folded
   sentence say so — "a clip flag those modes cannot honour, preset or not, …".

Grammar itself is fine: the `A, or B, now warns` subject is singular and the verb agrees, which is the
I-1 defect from both plan reviews correctly retired.

#### M-4 — the JSON's own comment now contradicts the invariant this commit establishes

`src/report.rs:867-869`, untouched by the commit:

```rust
// #440 — `null` unless a `--library` preset was selected. The four values
// repeated here are the ones in force, overrides included, so a consumer
// never has to reconstruct them from the preset name.
```

After this commit, "the ones in force" is exactly what the four JSON values are **not** on a single-end
run: the text report is now the "in force" surface and the JSON is deliberately the "as requested" one
(§Integration, Q-2). The next person to ask "should the JSON be gated too?" reads this comment and is
told the JSON already claims to be the in-force view — i.e. the opposite of the recorded decision.
CLAUDE.md's rule puts the invariant on the thing it constrains, and this is that thing. One-clause fix:
"The four values repeated here are the ones **requested**, overrides included; the text report carries
the ones in force."

### Low

- **L-1 — one CHANGELOG line left unwrapped by the re-wrap.** `CHANGELOG.md:23` is 111 characters;
  every other line in the paragraph is 73–97 (`awk` measurement, lines 8–30). Cosmetic only — there is
  no prettier or markdownlint config in the repo and no markdown step in `.github/workflows/`. Rewrap
  the paragraph.
- **L-2 — `single_end_announces_no_read_2_clip` is two negatives and one positive.**
  `tests/integration_library.rs:339`. It is **not** vacuous today (V8a below turns it red), and
  `report_text` panics if the report is missing, so it cannot pass on an absent file. But both
  `--clip_R2`-shaped assertions would survive someone gating the whole preset block away on single-end.
  One line closes it: `assert!(text.contains("Clipping in force: --clip_R1 10 --three_prime_clip_R1 10"))`.
- **L-3 — the #440 CHANGELOG entry is now ~270 words for one issue.** Against CLAUDE.md's "one sentence
  per entry" this is a pre-existing violation that the fold makes larger; noting it because the brief
  asks whether the resulting sentence is "one sentence's worth of claim" — the clause is, the entry
  containing it is not. Trimming the whole entry is a separate call.
- **L-4 — the new test file's module header argues against a rejected alternative.**
  `tests/integration_inert_clip_flags.rs:6-9`: "A single `!paired` test cannot express that — it reports
  … which sends the user to add `--paired`, silencing the warning without making the flag work." That is
  design justification in committed source, which CLAUDE.md discourages — but it is a faithful copy of
  the precedent it cites (`tests/integration_adapter2.rs:414-422` records an *earlier version of the
  fix* in the same way), so this is an observation, not a request.
- **L-5 — the hardtrim reason names both flags when only one was given.** `--hardtrim3 30 --clip_R1 3`
  says "`--hardtrim5/--hardtrim3` do not honour the clip flags". Verbatim the shape `-a2` already ships
  (`cli.rs:1053`), so this is consistency, not a defect. Recorded so the absence of a finding is a
  decision.

---

## Priority 1 — is the reason true in every reachable mode?

41 invocations, debug binary, every one with `-o <scratch dir>`. `rc` is the process exit code,
`warnlines` counts `is not used in this mode`.

```
=== A. trimming path ===
se_clip_r2            rc=0 warnlines=1  --clip_R2 5 … (it applies to Read 2 of a pair, and this is a single-end run)
se_clip_r1            rc=0 warnlines=0
se_tp_clip_r2         rc=0 warnlines=1  --three_prime_clip_R2 7 … (single-end run)
se_clip_r2_zero       rc=0 warnlines=1  --clip_R2 0 … (single-end run)
pe_clip_r2            rc=0 warnlines=0
pe_clip_all_four      rc=0 warnlines=0
se_multifile_clip_r2  rc=0 warnlines=1   ← two files, no --paired: one line for the invocation
=== B. specialty modes, no --paired (two files) ===
ht5_clip_r1           rc=0 warnlines=1  (--hardtrim5/--hardtrim3 do not honour the clip flags)
ht3_clip_r1           rc=0 warnlines=1  (--hardtrim5/--hardtrim3 do not honour the clip flags)
clock_clip_r1         rc=0 warnlines=1  (--clock does not honour the clip flags)
implicon_clip_r1      rc=0 warnlines=1  (--implicon does not honour the clip flags)
clock_all_four        rc=0 warnlines=4  all four named, all with the clock reason
=== C. specialty modes WITH --paired ===
ht5_paired_clip_r1    rc=0 warnlines=1  (hardtrim reason — not the single-end one)
clock_paired_clip_r2  rc=0 warnlines=1  (clock reason)
implicon_paired_tp2   rc=0 warnlines=1  (implicon reason)
=== D. specialty modes, single file ===
ht5_single_clip_r1    rc=0 warnlines=1  (hardtrim reason)
ht3_single_tp_r1      rc=0 warnlines=1  (hardtrim reason)
=== E. clump_only / clumpify ===
clumponly_clip_r2     rc=1 warnlines=0  Error: --clump_only does not clip; --clip_r2 is not compatible
clumponly_clip_r1     rc=1 warnlines=0  Error: --clump_only does not clip; --clip_r1 is not compatible
clumpify_se_clip_r2   rc=0 warnlines=1  (single-end reason) — with --cores 2 and gzip output
=== F. passthrough / demux ===
passthrough_clip_r2   rc=0 warnlines=0   ← --paired required, R2 clipping genuinely applies
demux_clip_r2         rc=0 warnlines=1  (single-end reason)
demux_clip_r1         rc=0 warnlines=0   ← demux runs after the trim, so R1 clipping applies
=== G. uBAM ===
ubam_in_se_clip_r2    rc=0 warnlines=1  (single-end reason)
ubam_in_pe_clip_r2    rc=0 warnlines=0   ← --paired + one interleaved uBAM
ubam_out_se_clip_r2   rc=0 warnlines=1  (single-end reason)
ubam_out_pe_clip_r2   rc=0 warnlines=0
ubam_out_interleaved  rc=0 warnlines=0
=== H. preset shapes ===
lib_se_only           rc=0 warnlines=0   ← the V6 guard, on the clean binary
lib_se_clip_r2        rc=0 warnlines=1  (single-end reason)
lib_pe_clip_r2        rc=0 warnlines=0
lib_ubam_se_clip_r2   rc=0 warnlines=1
lib_ubam_pe_clip_r2   rc=0 warnlines=0
lib_ubam_interleaved  rc=0 warnlines=0
```

Plus `multipair_pe` (4 files, `--paired --library emseq --clip_R2 5`) → 0 warnings, and `multifile_se`
(2 files, no `--paired`) → 1 warning.

**Every reason is true for its mode. No inert flag went unwarned. No honoured flag was warned about**
(subject to H-1, which is about the *test* for that property, not the property).

The specific thing the commit exists to fix is confirmed dead: `clock_paired_clip_r2` and
`ht5_paired_clip_r1` are genuinely paired runs and neither is told "this is a single-end run"; the
specialty modes name all four flags (`clock_all_four` → 4 lines).

Two supporting facts I verified rather than assumed:

- **`specialty.rs`, `demux.rs` and `clump_only.rs` read no clip flag.**
  `grep -rn "clip_r1\|clip_r2\|clip_R1\|clip_R2" src/specialty.rs src/demux.rs src/clump_only.rs` →
  zero matches. Assumption 5 holds.
- **`--clump_only` bails on all four**, not just the two the plan cites: `cli.rs:836` `clip_r1`,
  `:839` `clip_r2`, `:842` `three_prime_clip_r1`, `:847` `three_prime_clip_r2`, `:852` `library`.
- **`--clumpify` (not `--clump_only`) does reach the ladder** and warns correctly; two `validate()`-level
  refusals sit in front of it (`--cores >= 2`, and `--dont_gzip` mutual exclusion), and both stay quiet,
  which is assumption 4 behaving as designed.
- **`--implicon`'s bare form works in the test table.** `cli.rs:421-423` gives it
  `default_missing_value = "8", num_args = 0..=1, require_equals = true`, so the table's bare
  `--implicon` is `Some(8)` and does not leak a positional. (My own probes used `--implicon=8`.)

---

## Priority 2 — the gate, on both surfaces

`main.rs:1128` keys on `cli.paired`; `report.rs:287` on `config.paired`. These are different values, so
I checked whether they can ever disagree. `report::TrimConfig` is constructed at exactly three sites:

| Site | `paired` | Reached when |
|---|---|---|
| `main.rs:1464` (`run_single_file`) | `false` (literal) | the single-end loop, which `main.rs:893` guards with `if cli.paired { … } else` |
| `main.rs:2254` (`run_ubam_output_single`) | `false` (literal) | the uBAM-output single-end loop |
| `main.rs:2633` (`write_paired_reports`, shared) | `true` (literal) | every paired path, FASTQ and uBAM, including the interleaved single-file one |

So `config.paired` is not a copy of `cli.paired` but a property of the writer, and the two can only
disagree if a paired invocation reached a single-end writer or vice versa — which the `if cli.paired`
dispatch prevents. Measured across six preset shapes (`stderr_ovr` = the override sentence on stderr,
`txt_ovr` = the same sentence summed over that run's report files):

```
lib_se_only          stderr_ovr=0  txtfiles=1 txt_ovr=0 clipforce_lines=1  jsonfiles=1 json_ovr_entries=0
lib_se_clip_r2       stderr_ovr=0  txtfiles=1 txt_ovr=0 clipforce_lines=1  jsonfiles=1 json_ovr_entries=1
lib_pe_clip_r2       stderr_ovr=1  txtfiles=2 txt_ovr=2 clipforce_lines=2  jsonfiles=2 json_ovr_entries=2
lib_ubam_se_clip_r2  stderr_ovr=0  txtfiles=1 txt_ovr=0 clipforce_lines=1  jsonfiles=1 json_ovr_entries=1
lib_ubam_pe_clip_r2  stderr_ovr=1  txtfiles=2 txt_ovr=2 clipforce_lines=2  jsonfiles=2 json_ovr_entries=2
lib_ubam_interleaved stderr_ovr=1  txtfiles=2 txt_ovr=2 clipforce_lines=2  jsonfiles=2 json_ovr_entries=2
```

and the `Clipping in force:` lines, verbatim:

```
[lib_se_clip_r2]       Clipping in force: --clip_R1 10 --three_prime_clip_R1 10
[lib_ubam_se_clip_r2]  Clipping in force: --clip_R1 10 --three_prime_clip_R1 10
[lib_pe_clip_r2]       Clipping in force: --clip_R1 10 --clip_R2 5 --three_prime_clip_R1 10 --three_prime_clip_R2 10
                       --clip_R2 5 was given on the command line and overrides the emseq preset value 10
[lib_ubam_interleaved] (both R1 and R2 reports: same as pe_lib_ovr)
```

Multi-pair (4 files) gives 2 stderr banners and 4 reports, all four carrying the override; multi-file
single-end gives 0 and 0. **The two surfaces agree on every shape**, and the `Clipping in force:` line
(gated since `2142577`) agrees with both — the block is internally consistent for the first time.

`json_ovr_entries=1` on the two single-end preset shapes is the deliberate Q-2 asymmetry, and it is what
M-2/M-3/M-4 are about.

I also checked there is no fourth text surface: the only clip-naming writes in `report.rs` are the
`Clipping in force:` line (`:279-283`, via `flag_summary(config.paired)` at `:282`) and the override
loop (`:284-296`).
The JSON has three ungated clip blocks (`:877-890` preset-resolved, `:894-902` overrides, `:911-924`
raw request).

---

## Priority 3 — the `read_2` two-list hazard, and whether a third place needs it

There are now **three** encodings of "this flag belongs to Read 2":

1. `library.rs:124-128` — `flag_summary`'s tuple array. Pinned independently by
   `flag_summary_lists_all_four_for_paired` and `flag_summary_omits_read_2_for_single_end`
   (`library.rs:329-335`, asserting `"--clip_R1 10 --three_prime_clip_R1 10"`), both exact-string.
2. `library.rs:187-200` — `pick`'s four call sites. Pinned by the new oracle test
   (`every_override_knows_which_read_it_belongs_to`, `library.rs:297`), which uses `flag.ends_with("_R2")`
   as the oracle while production never inspects the string. **This guard works** — V12 below.
3. `cli.rs:1154-1165` — the warning loop, where the same fact is encoded *positionally* (R1 entries
   gated by `all_four`, R2 entries ungated). **No oracle, no test.** That is H-1.

The fourth candidate — the JSON `overrides` writer at `report.rs:894` — is a genuine consumer of
`ClipOverride` that ignores `read_2`, and that is correct per Q-2. It is the reason M-4 matters: the
decision is right, but the comment above it now says the opposite.

---

## Priority 5 — mutation spot-checks (I did not trust the green run)

All six re-run by me. Every one reverted with `git checkout --` and confirmed clean before the next.

**V8a — only the `main.rs` gate reverted** (`for o in &clips.overrides {`):

```
--lib                             491 passed; 0 failed        ← report.rs test correctly unaffected
--test integration_inert_clip_flags 4 passed; 0 failed
--test integration_library        9 passed; 1 failed
---- single_end_announces_no_read_2_clip stdout ----
panicked at tests/integration_library.rs:353:5:
stderr announces a Read 2 override on a single-end run:
WARNING: --clip_R2 5 was given but is not used in this mode (…single-end run). Ignoring.
…
Library preset 'emseq' selected: --clip_R1 10 --three_prime_clip_R1 10
--clip_R2 5 was given on the command line and overrides the emseq preset value 10
```

That panic body is the exact state `integration_adapter2.rs:414-422` calls the bug — warned it was
unused, then announced it had displaced the preset. The commit's central claim is real.

**V8b — only the `report.rs` gate reverted:** exactly 2 red, and *different* ones —
`report::tests::test_text_report_omits_r2_override_lines_on_single_end` **and**
`single_end_announces_no_read_2_clip`; `integration_inert_clip_flags` stayed 4/4 green. **The two
surfaces are independently covered**; neither gate can be reverted silently. This was the commit's
central claim and it holds.

**V12 — two `read_2` arguments swapped at `pick`'s call sites:** 6 red, including the oracle:
`library::tests::every_override_knows_which_read_it_belongs_to`, `library::tests::an_override_records_both_values`,
`report::tests::test_text_report_names_both_values_for_an_override`,
`report::tests::test_text_report_omits_r2_override_lines_on_single_end`,
`explicit_clip_flag_overrides_the_preset_and_is_logged`, `single_end_announces_no_read_2_clip`.

**V6 — the warning reads `effective_clips()`:** exactly the 2 §12 names —
`preset_alone_on_single_end_does_not_warn` and `explicit_clip_flag_overrides_the_preset_and_is_logged`.
The plan's central design decision is load-bearing and guarded.

**M4 (mine) — `all_four` narrowing removed:** 780 pass. See H-1.
**M5 (mine) — `--three_prime_clip_R1` entry deleted:** 780 pass. See M-1.

Seams I looked for and did **not** find: no new test passes without the code beneath it except the two
above; `report_text` panics on a missing report so the negative report assertions cannot pass on an
absent file; the new file's `tempdir(tag)` prefix (`tg_inert_`) cannot collide with
`integration_library.rs`'s (`tg_library_`); and the mode table's `err.contains(clip_flag)` cannot be
satisfied without the warning, since the `INERT` assertion runs first.

---

## Priority 4 — what the change should not have touched

Independent differential, base emulated by reverting **both** gates in the working tree (nothing else in
the commit can reach a file — verified by reading the diff: `trimmer.rs`, `parallel.rs`, `fastq.rs`,
`io.rs` are untouched). Both binaries built with `SOURCE_DATE_EPOCH=1700000000`, invoked as `./tg` from
their own directory with identical absolute paths, 7 shapes:

```
FILES COMPARED = 43   identical = 39   different = 4

--- DIFFERS: se_lib_ovr/BS-seq_10K_R1.fastq.gz_trimming_report.txt
14d13
< --clip_R2 5 was given on the command line and overrides the emseq preset value 10

--- DIFFERS: se_lib_ovr/_stderr
12d11
< --clip_R2 5 was given on the command line and overrides the emseq preset value 10

--- DIFFERS: ubam_se_lib_ovr/BS-seq_10K_R1.fastq.gz_trimming_report.txt
14d13
< --clip_R2 5 was given on the command line and overrides the emseq preset value 10

--- DIFFERS: ubam_se_lib_ovr/_stderr
12d11
< --clip_R2 5 was given on the command line and overrides the emseq preset value 10
```

**Four differences, all the same single line.** Every FASTQ, every `.bam`, every JSON, and every file on
all three paired shapes: byte-identical, including `ubam_se_lib_ovr/BS-seq_10K_R1_trimmed.bam`.
§12's V18 claim is confirmed and **extended**: my 7th shape (uBAM-output single-end) shows the gate also
reaches the second report site (`main.rs:2244`), which §12's five shapes did not cover. That is the
right behaviour, and it is the only respect in which my differential differs from §12's.

The JSON is unchanged on every shape, as intended.

---

## Verified correct — so a missing finding means a check, not an omission

- **`git rev-parse HEAD` = `3161b8a…`**, 8 files, +437/−18, and the diff touches no trimming, I/O or
  parallel code.
- **780 passed, 0 failed, 0 ignored, 18 binaries** (`cargo test --no-fail-fast`, macOS, debug) — 768 base
  + 12 new, matching §12. All 12 new test names present in the log. The new integration file raised the
  binary count 17 → 18, so cargo does pick it up.
- **`cargo fmt --all -- --check` → exit 0.**
- **`cargo clippy --all-targets -- -D warnings` → exit 0**, zero diagnostics.
- **`cargo clippy --all-targets --release -- -D warnings` → exit 0**, zero diagnostics, run from a
  cleaned target dir (207-line log, full dependency set, `trim-galore` checked at line 189). §12's
  release-clippy claim — the one thing both plan reviewers left owed — is confirmed. **And the harness
  can fail:** injecting a `len_zero` violation gives `PROBE_EXIT=101`, `error: length comparison to
  zero`, in both profiles.
- **V16's triple, on a clean binary:** stderr warns and names the flag; the text report contains **zero**
  occurrences of the string `clip`; the JSON has `clip_r2 = 5`, `clip_r1 = None`, `mode = single-end`.
- **The ladder's placement.** `cli.rs:1138` is after the last `bail!` in `validate()` (`:1132`), so no
  `validate()`-level refusal is preceded by the warning; `--clump_only`, `--clumpify --cores 1` and
  `--clumpify --dont_gzip` all confirmed silent-then-refuse.
- **`--paired` + one interleaved uBAM stays silent** and its two reports both name the R2 override —
  the one shape where `self.paired` is true but only one file exists.
- **`--rrbs` cannot interact:** `main.rs:1141` is paired-gated and binds a local, never `cli.clip_r2`.
- **`--clip_R2 0` warns** (assumption 2), and a multi-file single-end run warns **once per flag for the
  invocation**, not per file.
- **The ordering trio.** `hardtrim_conflict_precedes_the_library_refusal`,
  `clock_arity_error_precedes_the_library_refusal` and `missing_input_precedes_the_library_refusal` all
  pass, each asserting the expected substring **and** `!err.contains("--library")`. The `--clock --paired
  --library pbat R1` case does not slip `validate_paired_input`'s early `Ok`, because that path needs
  `mode_label == "Paired-end"` and `:1029` passes `"--clock"`.
- **D3's narrowing is sound and non-vacuous.** The override sentence is `"--clip_R2 5 was given on the
  command line …"`; the warning is `"--clip_R2 5 was given but is not used …"`. They share no substring
  long enough to confuse the assertion, so `!err.contains("--clip_R2 5 was given on the command line")`
  can distinguish the two, and V8a proves it does. Narrowing it was correct: the plan's original
  assertion would have contradicted the feature.
- **Comment convention.** `cli.rs:1138` and `library.rs:92` are one line each, present tense, state the
  fact, no history and no self-justification. Compliant.
- **`flags.md` register**: plain prose, no second person, no bold pseudo-heading scaffold, bullet shape
  matches its neighbours; the pre-existing `--library` bullet is unmodified. Only the "trimming report"
  wording (M-2) is wrong.
- **CHANGELOG grammar**: the new clause's singular subject and verb agree — the plural-verb defect both
  plan reviewers reported against v1 ("`--clip_R2 5` and `--three_prime_clip_R2 7` **was** given") is
  gone, replaced by one line per flag (`clock_all_four` → 4 separate lines).

**Not checked:** the `validation` CI job against Perl 0.6.11 (it passes no clip flag on any shape, and
warnings are stderr-only, so it cannot be affected); Linux test counts; and whether every one of the
plan's 23 outline steps landed — that is `plan-manager`'s audit, not this one.

---

## Worktree state

```
$ cd /private/tmp/claude-501/crA440 && git status --short
?? scratch_a/
$ git diff --stat HEAD
$ git rev-parse HEAD
3161b8a541c09d76d77e4e140f7b2638b8333d5c
```

No tracked file modified; `scratch_a/` (324 K) is my own untracked scratch inside a throwaway worktree —
the four driver scripts (`matrix.sh`, `surfaces.sh`, `diff2.py`, `differential.sh`) and the test/clippy
logs are kept, the regenerable run outputs and the two comparison binaries were deleted to give the
volume its space back. `test_files/` was never written to — every invocation passed `-o`. `cargo clean`
run on completion; the volume went from 8.4 GiB to 7.5 GiB free over the review. No GitHub writes;
nothing pushed, commented or opened.
