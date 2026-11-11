# CODE REVIEW — Reviewer B — commit `3161b8a`

**Target:** `3161b8a` on `fix/440-clip-r2-se-warning`, parent `fdade7e`
**Worktree:** `/private/tmp/claude-501/crB440`
**Plan:** `plans/08172026_clip-r2-se-warning/PLAN.md` v2 (+ §12 implementation notes)
**Reviewer:** B (independent; Reviewer A's report not read)

```
$ git rev-parse HEAD
3161b8a541c09d76d77e4e140f7b2638b8333d5c
$ git log --oneline -2
3161b8a fix(cli): warn when a clip flag cannot be honoured, and stop announcing it (#440)
fdade7e perf(clump): refit the --clumpify memory model against measured peak RSS (#439) (#446)
$ git diff --stat fdade7e 3161b8a
 8 files changed, 437 insertions(+), 18 deletions(-)
```

---

## Verdict — **READY WITH CHANGES**

The central claim holds and I could not break it. The reason ladder is mode-accurate on all 32
invocations I enumerated, including every shape the previous design got wrong; both override-line
surfaces are gated and each is covered by a test that fails when *only that gate* is reverted; paired
output is byte-identical in a real base-vs-new differential of 39 files; `clippy --release -D warnings`
and `cargo fmt --check` are clean, and I proved the clippy harness can fail before believing that it
passed. 780 tests pass, matching §12 exactly.

Three things to change, none of them behavioural:

1. **M-1** — `--three_prime_clip_R1` under the specialty modes has **no test**. I deleted its entire row
   from the warning array and all 780 tests stayed green. This is exactly the class Reviewer A's
   plan-review I-6 raised; v2 closed it for `--three_prime_clip_R2` and left the R1 3' flag open. One
   extra flag in one existing table row fixes it.
2. **M-2** — `flags.md`'s new bullet says "The trimming report lists only the clipping that applied".
   That is false of the JSON report, which records all four values *and* both R2 override entries on a
   single-end run (measured below). `library-presets.md` gets the txt/JSON distinction right; the
   terser restatement loses it.
3. **M-3** — the **third** place that encodes "which read a clip flag belongs to" is
   `main.rs:398`, and the commit did not teach it. A single-end `--rename --output-format ubam
   --clip_R2 5` run now prints *"…is not used in this mode … Ignoring."* and is then **refused because
   `--clip_R2` is set**. Pre-existing over-refusal; the commit's invariant is what makes it read as a
   contradiction.

Everything else is Low or a style observation.

---

## Defects

### Medium

#### M-1 — one of the four flags the commit newly warns about has no test at all

`src/cli.rs:1156-1164` (the `--three_prime_clip_R1` entry of the warning array) and
`tests/integration_inert_clip_flags.rs:51-124` (the mode table).

The `--three_prime_clip_R1` warning is reachable **only** through the `all_four` branch, i.e. only under
`--hardtrim5/3`, `--clock`, `--implicon` — the newest behaviour in the commit. The mode table's four
specialty rows all pass `--clip_R1`, and the single-end row passes `--clip_R2`;
`each_inert_flag_gets_its_own_line` covers `--three_prime_clip_R2`. Nothing anywhere sets
`three_prime_clip_r1` under a specialty mode.

Mutation, run in full with `--no-fail-fast` — the array row deleted entirely:

```
  PATCHED: --three_prime_clip_R1 row removed from the warning array
  cargo exit=0
  passed=780 failed=0
  FAILCOUNT=0
```

Contrast every other mutation I ran, each of which turned something red. The **behaviour is correct** —
I verified it live:

```
### implicon_paired_all4               rc=0 warnlines=4
      W| WARNING: --clip_R1 1 ... (--implicon does not honour the clip flags). Ignoring.
      W| WARNING: --three_prime_clip_R1 3 ... (--implicon does not honour the clip flags). Ignoring.
      W| WARNING: --clip_R2 2 ... (--implicon does not honour the clip flags). Ignoring.
      W| WARNING: --three_prime_clip_R2 4 ... (--implicon does not honour the clip flags). Ignoring.
```

so this is a regression-protection gap, not a bug. Reviewer A's plan-review I-6 predicted it in these
words: *"Half the feature (`--three_prime_clip_R2`) is never exercised by any proposed test, so deleting
its `if let` block would keep all 776 tests green."* v2 answered that for the flag named and not for its
sibling.

**Fix:** add `--three_prime_clip_R1` to one specialty row of the table (e.g. give the `implicon` row two
flags and assert two warning lines), or add a sixth row. One row also closes the neighbouring gap: **no
test passes a specialty mode together with `--paired`**, so a future edit that nested the ladder inside
`if !self.paired` would silently stop warning on `--clock --paired` with nothing going red. The table's
specialty rows deliberately omit `--paired`, which is right for reason-accuracy but leaves that arm
untested.

#### M-2 — `flags.md`'s report sentence is false of the JSON report

`docs/src/content/docs/guide/flags.md:15`, final sentence:

> The trimming report lists only the clipping that applied, so a single-end report names Read 1 alone.

Measured on the single-end run that exercises all four overrides
(`--library accel --clip_R1 1 --clip_R2 2 --three_prime_clip_R1 3 --three_prime_clip_R2 4`):

```
BS-seq_10K_R1.fastq.gz_trimming_report.json
  mode = single-end
  library = {"preset": "accel", "clip_r1": 1, "clip_r2": 2, "three_prime_clip_r1": 3,
             "three_prime_clip_r2": 4, "overrides": [{"flag": "--clip_R1", ...},
             {"flag": "--clip_R2", ...}, {"flag": "--three_prime_clip_R1", ...},
             {"flag": "--three_prime_clip_R2", ...}]}
  params clip: {'clip_r1': 1, 'clip_r2': 2, 'three_prime_clip_r1': 3, 'three_prime_clip_r2': 4}
```

There are two trimming reports, and the JSON lists all four values plus both R2 override entries. That
is the plan's deliberate decision (Q-2) and `docs/.../library-presets.md:86-88` states it exactly:
*"The JSON block carries all four resolved values whatever the run type … That is why the same run's
text and JSON reports differ above."* The new bullet, four lines below the one that links to that page,
contradicts it.

Second, weaker half: for a run **without** `--library`, the text report lists no clipping at all —
`report.rs:275-276` gates the whole block on `config.library` *and* `library.preset`. Verified:

```
$ grep -n -i "clip" <se_trim_clipR2>/BS-seq_10K_R1.fastq.gz_trimming_report.txt   → no matches
$ grep -n -i "clip" <pe_trim_clipR2>/BS-seq_10K_R1.fastq.gz_trimming_report.txt   → no matches
```

So the sentence promises a listing that most runs do not get.

**Fix:** name the surface and the condition, e.g. *"The text report lists only the clipping a
`--library` preset put in force, so a single-end report names Read 1 alone; the JSON report records all
four values whatever the run type."*

#### M-3 — the third place that knows which read a clip flag belongs to, and it disagrees with the warning (pre-existing)

`src/main.rs:398`:

```rust
&& (cli.effective_clips().any_set() || cli.hardtrim5.is_some() || cli.hardtrim3.is_some())
```

`any_set()` is read-agnostic, so on a single-end run an R2 flag alone satisfies it. Measured:

```
### rename_ubamout_se_clipR2           rc=1 warnlines=1
      W| WARNING: --clip_R2 5 was given but is not used in this mode (it applies to Read 2 of a
         pair, and this is a single-end run). Ignoring.
      E| Error: --rename is refused with --output-format ubam when any input is FASTQ. ...
### rename_ubamout_se_noclip           rc=0 warnlines=0
```

The second row is the proof: the identical invocation without `--clip_R2 5` succeeds. So the tool says
it is ignoring the flag and then refuses the whole run because of it — and the stated remedy ("Use
FASTQ output, or drop --rename") does not mention the flag that actually triggered it.

`main.rs:391-393`'s own comment says the list "must match every `--rename`-driven `append_to_id` site
… without one of them set, `--rename` appends nothing and there is nothing to lose". On a single-end run
an R2 clip flag drives no `append_to_id` site, so by that comment's own rule it should not count. Not
introduced here — `main.rs:398` is untouched by the diff — but the brief asked whether any third place
now needs to know, and this is it. The narrow fix is a read-aware predicate at that site
(`clip_r1.is_some() || three_prime_clip_r1.is_some() || (cli.paired && …)`), which is a separate change;
the minimum is to record it, since the commit's invariant is what makes the pairing read as a
contradiction.

### Low

#### L-1 — the CHANGELOG's "either surface" names neither, in a paragraph that describes three

`CHANGELOG.md:20-22`:

> — and a clip flag those modes cannot honour, or a Read 2 clip flag on a single-end run, now warns
> instead of passing silently, **with nothing on either surface** announcing a Read 2 clip the run did
> not perform.

Six lines later the same paragraph describes *three* report surfaces (stderr log, text report, JSON),
and the third does record the value. A reader deciding whether an upgrade affects them cannot tell
which two "either" refers to. Naming them ("…on the log or in the text report…") costs three words and
removes the ambiguity. The clause is otherwise accurate: I verified both text surfaces on eight shapes.

Minor register point on the same clause: *"a clip flag … now warns"* attributes the warning to the flag
rather than the tool.

#### L-2 — a re-wrap leftover in the CHANGELOG

`CHANGELOG.md:23` is **111 characters** against 84–97 for every other line of the paragraph:

```
  8  74 | 16  88 | 20  97 | 23 111  <-- outlier
  9  92 | 17  85 | 21  96 | 24  93
 10  84 | 18  94 | 22  94 | 25  73
```

The inserted clause was wrapped at :19-22 and the overflow was pushed into :23 without re-wrapping it.
Cosmetic; the file is otherwise consistent.

#### L-3 — one of the 12 new tests passes with the feature entirely absent

`src/cli.rs:2589` `inert_clip_flag_is_a_warning_not_an_error` asserts only that `validate()` returns
`Ok`. Delete the whole warning block and it still passes. That is legitimate for a negative test — the
plan says as much at step 5, since `eprintln!` is invisible to the unit harness — but it means the
"12 new tests" figure contains one row with zero coupling to the feature. No change needed; worth
knowing when the number is quoted.

The three ordering tests have the same shape in a milder form: they pin the refusal's *precedence*, not
its existence, and would pass if the `--library` refusal block were deleted. That is correct division of
labour — `library_plus_{hardtrim5,hardtrim3,clock,implicon}_rejected` at `cli.rs:2457-2516` pin the
refusal itself — but the two sets are only jointly sufficient.

---

## Style observations (not defects)

#### S-1 — the new test file's module doc argues against a design that was never shipped

`tests/integration_inert_clip_flags.rs:4-15`:

```
//! Two inert sets, and they differ: … A single `!paired` test cannot express that — it
//! reports "this is a single-end run" for a two-file `--clock` run, which sends the
//! user to add `--paired`, silencing the warning without making the flag work.
//!
//! Modelled on `integration_adapter2.rs`, … The specialty rows
//! deliberately pass two files and no `--paired`, since that is the shape a bare
//! `!paired` predicate gets wrong.
```

CLAUDE.md, twice over: *"Current state only — never reference the state you are fixing"* and *"Never
justify the change in the source. If you find yourself pre-empting 'but doesn't this break X?', that
argument belongs in the commit message or the PR thread."* The rule is scoped to "code, tests and
config that gets committed", so tests are in scope. Roughly half of these twelve lines defends the
ladder against an alternative the repo never shipped.

I am ranking this an observation rather than a finding because the repo's own precedent cuts the other
way: `tests/integration_adapter2.rs:414-422` does exactly this, the plan cites it approvingly as where
"hard-won reasoning" lives, and the plan's Self-Review makes that a stated principle. That is a
convention question for Felix, not something a reviewer should decide. If the rule wins, the
load-bearing sentence to keep is the first — *the two inert sets differ, and why* — and the two
`!paired`-predicate paragraphs are what the commit message already carries.

The other new doc comments are compliant: `library.rs:92`, `cli.rs:1138`, `report.rs:1428-1429`,
`library.rs:294-295` and `tests/integration_inert_clip_flags.rs:47-50` each state a current fact in one
or two lines. `cli.rs:2585-2586`'s *"which is refused because it is new"* is a rationale in source, but
it is one clause.

#### S-2 — the two-list invariant is written down nowhere in the source

`library.rs:93`'s field doc (*"A Read 2 flag, which a single-end run cannot honour"*) is exactly right
and one line. But the hazard the plan identifies — `pick`'s four call sites (`library.rs:186-201`) and
`flag_summary`'s triple array (`library.rs:124-128`) encode the same fact and can drift — is recorded
only in the plan and in the oracle test's doc comment. CLAUDE.md's *"put an invariant on the thing it
constrains, not in the function that relies on it"* would put a clause at one of the two lists.

Optional, and arguably already answered: a mechanical guard beats a comment, and I confirmed the guard
is tight (see V12b below). Noting it only because the plan flags the hazard as its main new risk and the
source does not.

#### S-3 — the CHANGELOG entry is now ~290 words

Folding into #440's existing bullet was the right call (both plan reviewers' O-4/O-5 disagreed, and
folding is what avoids two bullets citing one unreleased issue). The result is a four-sentence,
~290-word entry against CLAUDE.md's *"One sentence per entry. No exceptions."* Pre-existing debt that
this commit adds ~45 words to. If the entry is ever trimmed, the added clause is one sentence's worth
of claim and survives on its own.

---

## Verified as correct — so a missing finding means a check, not an omission

### 1. The reason is true in every reachable mode (brief priority 1)

32 invocations, each with its own `-o` scratch dir, run against the debug binary from the crate root.
Full log: `/private/tmp/claude-501/crB440_probe.log`. Every row's reason names the mode that will
actually run; `main.rs`'s dispatch order (473 hardtrim5, 504 hardtrim3, 535 clock, 557 implicon) is the
same order as the ladder, so the first-match arm is the arm that executes.

| Invocation | warn lines | Reason given | Correct? |
|---|---|---|---|
| `--paired --clip_R2 5 R1 R2` | 0 | — | ✅ R2 clipping applies |
| `--clip_R2 5 R1` | 1 | single-end | ✅ |
| `--clip_R1 5 R1` | 0 | — | ✅ R1 clipping applies |
| `--clip_R2 0 R1` | 1 | single-end | ✅ explicit zero still an instruction |
| `--clip_R2 5 --three_prime_clip_R2 7 R1` | 2 | single-end ×2 | ✅ one line per flag |
| `--clip_R2 5 R1 R2` (no `--paired`, 2 SE runs) | 1 | single-end | ✅ once per invocation |
| `--hardtrim5 30 --clip_R1 3 R1` | 1 | hardtrim5/3 | ✅ |
| `--hardtrim5 30 --clip_R1 3 R1 R2` | 1 | hardtrim5/3 | ✅ **the shape the old design got wrong** |
| `--paired --hardtrim5 30 --clip_R2 3 R1 R2` | 1 | hardtrim5/3 | ✅ |
| `--hardtrim3 30 --clip_R1 3 R1 R2` | 1 | hardtrim5/3 | ✅ |
| `--clock --clip_R1 3 CR1 CR2` | 1 | clock | ✅ |
| `--clock --paired --clip_R2 3 CR1 CR2` | 1 | clock | ✅ |
| `--implicon=8 --clip_R1 3 CR1 CR2` | 1 | implicon | ✅ |
| `--implicon=8 --paired` + all four | 4 | implicon ×4 | ✅ all four named |
| `--clump_only --clip_R2 5 R1` | — | refused at `cli.rs:839` | ✅ unreachable |
| `--paired --passthrough I1 --clip_R2 5 R1 R2` | 0 | — | ✅ R2 clipping applies |
| `--demux <sheet> --clip_R2 5 <se>` | 1 | single-end | ✅ demux is post-trim (`main.rs:1516`) |
| `--demux <sheet> --clip_R1 5 <se>` | 0 | — | ✅ R1 clipping applies |
| `--output-format ubam --clip_R2 5 R1` | 1 | single-end | ✅ |
| `--paired --output-format ubam --clip_R2 5 R1 R2` | 0 | — | ✅ |
| `--clip_R2 5 ubam_test.bam` | 1 | single-end | ✅ |
| `--paired --clip_R2 5 ubam_paired_test.bam` (interleaved) | 0 | — | ✅ R2 exists after de-interleaving |
| `--paired --output-format ubam --clip_R2 5 <interleaved>` | 0 | — | ✅ |
| `--paired --clip_R2 5 <one FASTQ>` | 0 | refused: "only legal if that file is a uBAM" | ✅ no false claim |
| 8 `--library` shapes | as tabulated below | — | ✅ |

- **No mode is missing a warning.** The only ways a clip flag can be inert are the four specialty modes
  (all covered), the two R2 flags on a single-end run (covered), and `--clump_only`, which bails on all
  four (`cli.rs:835-849`) plus `--library` (`:851`) long before the ladder.
- **No mode gets a false warning.** `specialty.rs` reads none of the four flags —
  `grep -n "clip_r1\|clip_r2\|three_prime" src/specialty.rs` → **zero matches** — confirming assumption
  5. And all four modes do clip at fixed positions of their own (`clip_5prime` at `specialty.rs:107,219`;
  clock clips R1@13/R2@15 per `:241-243`; implicon clips R2 by `umi_length` per `:350-352`), so the
  chosen wording *"do not honour the clip flags"* is accurate where *"performs no clipping"* would have
  been false.
- **`--rrbs` cannot interact**: `main.rs:1141` is paired-gated and binds a **local** `clip_r2`, never
  `cli.clip_r2`.
- **Reads the raw fields, not `effective_clips()`**: `--library emseq <se>` with no clip flag typed is
  silent (`preset_alone_on_single_end_does_not_warn`'s subject). §12's V6 records the mutation; I did
  not re-run it, but both plan reviewers reproduced it independently on the identical snippet and I
  confirmed the current behaviour is silent.

### 2. The gate, on both surfaces (brief priority 2)

- **Only two surfaces print an override line.** `grep -rn "was given on the command line and overrides"
  src/ tests/` → `report.rs:291` and `main.rs:1130` only (plus two doc examples, both `--clip_R1`, so
  unaffected). Both are gated.
- **`cli.paired` vs `config.paired` is right at every site.** `report::TrimConfig.paired` is a
  **literal**, not `cli.paired`: `false` at `main.rs:1464` (`run_single_file`) and `main.rs:2254`
  (`run_ubam_output_single`), `true` at `main.rs:2633` (`write_paired_reports`). Tracing the dispatch:
  `run_single_file` is reachable only from the `else` of `if cli.paired` (`main.rs:1039`);
  `run_ubam_output_single` only after both paired arms have returned (`main.rs:2168`); every paired
  arm — two-file FASTQ, interleaved-uBAM FASTQ-out, two-file uBAM-out, interleaved uBAM-out — routes
  through `write_paired_reports`. So the banner and the report **cannot disagree**, and I confirmed it on
  8 shapes including multi-file SE, uBAM-out SE, and interleaved-uBAM paired.
- **Every other R2-announcing surface was already gated.** `Clipping in force` uses
  `flag_summary(paired)` on both surfaces (`main.rs:1126`, `report.rs:282`); the `--rrbs` auto-clip note
  (`main.rs:1141-1143`) is paired-gated. `grep` for R2 clip values in any `writeln!`/`eprintln!` across
  `main.rs`/`report.rs`/`trimmer.rs`/`parallel.rs` finds nothing else.
- **The invariant holds, measured.** On `--library accel` single-end with all four flags:
  `Clipping in force: --clip_R1 1 --three_prime_clip_R1 3`, two override lines (R1 only), and
  `grep -c -e "--clip_R2" -e "--three_prime_clip_R2"` over the text report = **0**. On the paired twin
  and on the interleaved-uBAM paired shape, both R2 lines are present on both surfaces.

### 3. The `read_2` two-list hazard (brief priority 3)

- **One production constructor.** `grep -rn "ClipOverride {" src/` → `library.rs:174` (inside `pick`)
  and `library.rs:270` (a test). No third builder.
- **The oracle guard is tight.** I ran the narrowest possible drift — flipping **only**
  `three_prime_clip_r1`'s `read_2` from `false` to `true`:

  ```
    PATCHED 193:            true,      # library.rs, the --three_prime_clip_R1 pick call
    passed=779 failed=1
      test library::tests::every_override_knows_which_read_it_belongs_to ... FAILED
  ```

  Exactly one test catches it, and it is the oracle. Nothing else in the suite would have.
- **The other list is guarded too.** `flag_summary`'s triple array (`library.rs:124-128`) is pinned by
  two exact-string tests (`flag_summary_lists_all_four_for_paired`,
  `flag_summary_omits_read_2_for_single_end`), so a flip there fails on the string, not on a bool.
- **The new field does not leak into the JSON schema.** `report.rs:894-903` writes
  `flag`/`preset_value`/`user_value` only; the emitted JSON carries no `read_2` key (verified on a
  four-override run). No consumer-visible schema change.
- **Third place:** see M-3. That is the only one I found.

### 4. What the change should not have touched (brief priority 4)

I built the **base binary from `fdade7e`** in the same worktree (`git checkout fdade7e -- src/`,
`cargo build --bin trim_galore`, source restored) and ran both binaries **through the same path**, so
the report's `command_line` field is identical except for the `-o` directory, which I excluded by line.
Harness proved before use: base emits 0 warning lines and 1 override line on `se_emseq_clipR2`; new
emits 1 and 0.

```
pe_accel_all4              files= 6 differing=0
pe_emseq_clipR2            files= 6 differing=0
pe_emseq_only              files= 6 differing=0
pe_noprese_clipR2          files= 6 differing=0
pe_ubam_il_emseq_r2        files= 6 differing=0
se_accel_all4              files= 3 differing=1   (two lines deleted: --clip_R2, --three_prime_clip_R2)
se_emseq_clipR1            files= 3 differing=0
se_emseq_clipR2            files= 3 differing=1   (one line deleted: --clip_R2)
===== compared=39 identical=37 differing=2 (command_line excluded) =====
```

- **All five paired shapes byte-identical**, including one the plan's own differential did not have:
  `--paired --library emseq --clip_R2 5` on an interleaved uBAM.
- **Every FASTQ output identical**, on every shape. No shape produced zero files (checked, not assumed).
- **Every JSON identical** — the deliberate no-change is confirmed, not asserted.
- **The only behavioural change anywhere is line deletion in the single-end text report**, on exactly
  the two shapes that override an R2 flag. §12's claim is verified and extended from 24 to 39 files.
- The `validation` CI job is unaffected: it passes no clip flag on any shape, and the warning is
  stderr-only.

### 5. Test quality — mutations run, not trusted (brief priority 5)

Baseline: `cargo test --no-fail-fast`, 18 targets, **780 passed / 0 failed** — matching §12 exactly.
Every mutation below was applied to a `cat`-snapshot, run as a **full** `--no-fail-fast` suite, then
restored and `cmp`-verified byte-identical against the snapshot (all four files, every run).

| Mutation | Result | §12 says |
|---|---|---|
| **V8a** — revert **only** the `main.rs` gate | **1 red**: `single_end_announces_no_read_2_clip`, at `integration_library.rs:353` (the *stderr* assertion) | 1, same test ✅ |
| **V8b** — revert **only** the `report.rs` gate | **2 red**: `test_text_report_omits_r2_override_lines_on_single_end` at `report.rs:1447`, and `single_end_announces_no_read_2_clip` at `integration_library.rs:363` (the *report* assertion) | 2 ✅ |
| **V2** — ladder loses its three specialty arms | **1 red**: `inert_clip_flag_warns_with_a_mode_accurate_reason`, first failing row `[hardtrim5]` | table red ✅ |
| **V12** — swap `clip_R1`/`clip_R2`'s `read_2` at `pick` | **6 red**, incl. `every_override_knows_which_read_it_belongs_to` | 6 incl. the oracle ✅ |
| **V12b** (mine) — flip **only** `three_prime_clip_r1`'s `read_2` | **1 red**, the oracle alone | — |
| **V19** (mine) — delete the `--three_prime_clip_R1` array row | **0 red, 780 green** | — → **M-1** |

The commit's central claim — that the two surfaces are covered **independently** — holds, and the
failure *lines* prove it rather than the test names: V8a fails on the stderr assert, V8b on the report
assert, in the same test.

### 6. Deviation D3 (brief priority 6)

The narrowing is correct and the assertion is not vacuous. `single_end_announces_no_read_2_clip`
(`integration_library.rs:339-366`) asserts three things: the override *sentence*
`"--clip_R2 5 was given on the command line"` is absent, `"is not used in this mode"` **is** present,
and the report contains no `--clip_R2` at all. V8a and V8b each turn it red, from different halves, so
it cannot pass vacuously. The plan's literal wording (stderr contains no `--clip_R2`) would indeed have
contradicted the feature, since the warning names the flag on purpose — D3's reasoning is right.

The one residual: the negative assertion is keyed on the override sentence's wording, so a reword of
`report.rs:291`/`main.rs:1130` would make it pass vacuously. `paired_still_announces_r2_overrides`
(`:371-389`) asserts the **presence** of the same substring on the paired twin, so the pair pins the
wording between them. Adequate.

### 7. Conventions (brief priority 7)

```
$ cargo clippy --all-targets --release -- -D warnings
CLIPPY_RELEASE_EXIT=0
grep -cE "^(warning|error)" → 0
   Compiling trim-galore v2.3.0 (/private/tmp/claude-501/crB440)
    Finished `release` profile [optimized] target(s) in 2.26s
$ cargo fmt --all -- --check
FMT_EXIT=0
```

Release **confirmed** (both plan reviewers owed it). I did not trust the clean run: I appended a
`return true;` to `src/lib.rs` and re-ran, which failed —

```
error: unneeded `return` statement
error: could not compile `trim-galore` (lib) due to 1 previous error
error: could not compile `trim-galore` (lib test) due to 1 previous error
```

— then restored the file and confirmed `git status --short` empty. So the harness can fail, and the
clean result means something.

Also verified:

- **`--help` agrees with the new warning.** `--clip_R2` / `--three_prime_clip_R2` both read
  "(paired-end only)", so the flag's own documentation and the runtime message say the same thing.
- **`flags.md`'s new bullet matches the docs register** — plain prose, bold only on the flag names
  exactly as the nine sibling bullets do, no conversational heading, no second person, no
  bold-pseudo-heading scaffold. Its factual problem is M-2, not its register.
- **`flags.md`'s mode claims are accurate**: "all four under `--hardtrim5`, `--hardtrim3`, `--clock` and
  `--implicon`, which do their own fixed clipping" — verified against `specialty.rs` and live.
- **`library-presets.md` needs no update**: both of its worked examples use `--clip_R1 12` on single-end,
  which the gate leaves untouched, and its `Overrides` section already explains that R2 values are
  reported for paired runs only.
- **CHANGELOG**: fold-into-#440 was the right structural call; the resulting clause is true of both text
  surfaces and grammatical, with the ambiguity at L-1 and the wrap at L-2.
- **Efficiency**: the ladder is ≤4 `Option::is_some()` checks and one stack-allocated 4-element array
  walk, once per invocation (`validate()` has a single call site, `main.rs:242`). The two gated loops add
  one `bool` test per override, ≤4 per run. Nothing on a hot path; no allocation added.
- **Test placement** matches the plan: the three ordering tests and
  `inert_clip_flag_is_a_warning_not_an_error` sit in the `── #440` section (`cli.rs:2454`), directly
  after `hardtrim5_range_error_precedes_the_library_refusal` (`:2517`).
- **No test writes into `test_files/`**: every new invocation passes `-o <tempdir>`, and
  `git status --short` is empty after a full suite run plus 6 mutation runs.

### 8. Checked and found nothing

- The warning printing **above** the version banner (`cli.validate()?` is `main.rs:242`, the banner
  `:249-251`) is unchanged `-a2` behaviour and a design point Felix already ruled on; not re-litigated.
- `--paired` with an odd file count, `--paired` with one FASTQ, `--clump_only` with any clip flag, and
  `--library` with any specialty mode all refuse **before** the ladder, so none of them can print a
  claim about a flag.
- `--basename`, `--retain_unpaired`, `--no_report_file`, `--cores` touch none of the changed code.
- The three new ordering tests exercise the paths the plan claims (`"cannot be combined in one
  invocation"`, `"requires an even number of input files"`, `"Input file not found"`), each also
  asserting `!err.contains("--library")`; all pass, and V15's equivalent is recorded in §12.

**Not checked:** Linux (macOS only; the plan's "never gate on one test number" caveat still applies).
The `validation` CI job was reasoned about from the workflow rather than executed. V6 and V10 and V15
were not re-run by me — V6's behaviour I confirmed live and both plan reviewers reproduced its mutation;
V10 and V15 are recorded in §12 and I verified V10's guard exists (`paired_still_announces_r2_overrides`
asserts presence on both surfaces).

---

## Worktree state

```
$ cargo clean
     Removed 4976 files, 1.3GiB total
$ git status --short          # run, not asserted — empty
$ git rev-parse HEAD
3161b8a541c09d76d77e4e140f7b2638b8333d5c
$ ls -d target
ls: target: No such file or directory
$ df -h . | tail -1
/dev/disk3s5   926Gi   891Gi   7.1Gi   100%    2.9M   75M    4%   /System/Volumes/Data
```

Clean: no tracked file modified, no untracked file left in the worktree, HEAD still `3161b8a`,
`target/` removed. Scratch output, logs and the base binary live outside the worktree under
`/private/tmp/claude-501/crB440_{probe.log,scratch,diff,mut_*.log,clippy_rel.log}`.

No GitHub writes were made; nothing was pushed, commented or opened.
