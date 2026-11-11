# PLAN REVIEW — Reviewer A (generation g1)

**Plan:** `plans/08162026_pr441-library-presets/PLAN.md`
**Base verified:** `/private/tmp/claude-501/pr441-A` at `145432e179ae614814dd9754846b76c77137127b`.
Tree clean at review start and at review end (`git status --short` → empty output, quoted in §5).
**Reviewer:** A, independent. Reviewer B's plan report not read.
**Settled decisions:** not relitigated. None of the four is unimplementable as specified.

---

## Verdict

**READY WITH CHANGES.**

The plan is unusually well-grounded: **all 14 line-number anchors it cites resolve to what it says
they are** (12 exact; 2 off by one or two lines at the end of a range, neither wrong about what it
points at), and the three claims most likely to be wrong — the two `flag_summary` call
sites, `both_reports_record_the_preset_and_its_expanded_values` being single-end, and T-1/F-5a being
genuinely blind — I confirmed empirically, not by reading. I applied the plan's own B2 patch to my
worktree and ran the whole suite with `--no-fail-fast`: **751 tests, exactly one failure**, the one the
plan predicts. Its self-review's two "caught" items are real catches.

Four things should change before implementation:

1. **`CHANGELOG.md:21-22` — a sentence in #441's own entry is falsified by B2 and is not on the fix
   list.** The plan's §7 rule ("docs changes required for factual correctness by our behaviour change
   are in scope") covers it; the plan applies that rule to `library-presets.md` and not to the
   changelog. Both entries land in the same `### Unreleased` block.
2. **The stderr half of B2 is not covered by any test — it is the one fix of the seven that is
   silently revertible.** The plan calls the two call sites its "load-bearing detail" and then tests
   only one of them.
3. **`library_plus_implicon_rejected` as written tests the wrong thing** (`--implicon 8` puts `8` in
   `input`). Proven; the fix is to copy the existing `fastqc_plus_implicon_rejected` shape.
4. **V5 is vacuous:** `cargo test flag_summary` matches zero tests today and would still match zero
   after implementation.

None of these is architectural. With them addressed the plan is implementable as written.

---

## 1. Claim-verification table

Every claim the lead asked about, plus each line anchor. Commands and outputs are quoted in §2-§5.

### 1a. Line-number anchors (all from `pr441-A` at `145432e`)

| Plan anchor | Claimed to be | Result |
|---|---|---|
| `cli.rs:1071-1097` | #421 `--fastqc` refusal on the four modes | **VERIFIED**, precise. `1071` = `// #421 — these arms reach no \`fastqc::run\` call site…`, `1073` = `if self.fastqc_requested() {`, `1097` = closing brace of the `implicon` bail. The four per-mode bails are exactly `1074-1097`. |
| `cli.rs:2309-2350` | that block's four `fastqc_plus_<mode>_rejected` tests | **VERIFIED** (off-by-one on the end). `2309` = the `── #421: --fastqc refused …` section comment; the four tests run `2313-2351`, not `-2350`. Substantively right. |
| `cli.rs:2466` | `hardtrim5_range_error_precedes_the_fastqc_refusal` | **VERIFIED**, exact. `2466` is the `fn` line. |
| `cli.rs:833-850` | the `--clump_only` clip family whose wording was selected | **VERIFIED**, exact. `833` = `if self.clip_r1.is_some() {`, `850` = `anyhow::bail!("--clump_only does not clip; --library is not compatible");`. |
| `cli.rs:1009` | #386 `--hardtrim5` + `--hardtrim3` bail | **VERIFIED**. `1009` = `if self.hardtrim5.is_some() && self.hardtrim3.is_some() {`. |
| `cli.rs:1016-1024` | the two hardtrim range checks | **VERIFIED**. `1016` = `if let Some(n) = self.hardtrim5`, `1019`/`1024` = the two `1 and 999` bails. |
| `cli.rs:1026-1031` | `validate_paired_input` for `--clock` / `--implicon` | **VERIFIED**, exact (`1027` and `1030` are the two calls). |
| `main.rs:273-277` | the `--phred64` uniform-refusal precedent | **VERIFIED**, exact. `273-277` is the passage `// Rejected uniformly rather than per-mode … one refactor away from being wrong.` |
| `main.rs:393` | the #408 guard widened to read `effective_clips()` | **VERIFIED**, exact. `393` = `&& (cli.effective_clips().any_set() \|\| cli.hardtrim5.is_some() \|\| cli.hardtrim3.is_some())`. |
| `main.rs:1121` | the stderr log call site | **VERIFIED**, exact (`clips.flag_summary()`). |
| `report.rs:279` | the text-report call site | **VERIFIED**, exact (`writeln!(w, "Clipping in force: {}", library.flag_summary())?;`). |
| `library.rs:121-134` | the current `flag_summary` body | **VERIFIED**, exact — `121` is the `fn` line, `134` the closing brace. |
| `library-presets.md:48-54` | the single-end override transcript | **VERIFIED**, exact (`48` command, `51-54` fence). |
| `library-presets.md:62-66` | the single-end text-report transcript | **VERIFIED**, exact. |

### 1b. Substantive claims

| # | Claim | Result |
|---|---|---|
| 2 | `flag_summary` has exactly **two** call sites | **VERIFIED** — grep over `src tests docs CHANGELOG.md` returns 3 hits total: the definition plus `report.rs:279` and `main.rs:1121`. No test calls it. No third site. |
| 3 | `both_reports_record_the_preset_and_its_expanded_values` is **single-end**, asserts all four, must be updated | **VERIFIED empirically.** It is `--dont_gzip --library pbat -o dir R1` (no `--paired`). Applying the plan's gate made it the **only** failure in 751 tests. |
| 4 | The refusal's placement achieves the claimed ordering | **VERIFIED.** The `fastqc_requested` block closes at `cli.rs:1110`; between there and `Ok(())` at `1140` there are only `eprintln!` deprecation warnings — nothing fallible. So the new block genuinely is last. `--hardtrim5 0` (range, `:1019`), `--clock`/`--implicon` arity (`:1027`/`:1030`) and #386 (`:1009`) all precede it. |
| 4b | `--clump_only --hardtrim5 --library emseq` | **ANSWERED** (plan does not state it): the user gets `--clump_only does not clip; --library is not compatible`, because `:850` precedes the specialty mutual-exclusion bails at `:893-903`, and the whole `clump_only` block at `:780` precedes the new block. See finding I-6. |
| 5 | "does not clip" is accurate for all four specialty modes | **VERIFIED structurally** — `src/specialty.rs` contains **zero** occurrences of `clip_r1`, `clip_r2`, `three_prime_clip*`, `effective_clips` or `clips.`. See finding O-3 on the wording. |
| 6 | The four-distinct-values lever catches every key permutation and fixes T-1 + T-2 + the ≥2 comma path | **VERIFIED.** `resolve()` evaluates `pick` in struct-literal field order, so `overrides` is deterministically `[--clip_R1(10→1), --clip_R2(15→2), --three_prime_clip_R1(10→3), --three_prime_clip_R2(10→4)]`. Resolved values `1/2/3/4` and displaced values `10/15/10/10` are independently assertable, exactly as claimed. **T-1 and F-5a confirmed real by mutation** (§3). |
| 7 | Two docs transcripts are single-end and would be falsified | **VERIFIED**, and I captured the real replacement output (§4). Incomplete — see I-3. |
| 8 | Refusing breaks no existing test | **VERIFIED.** Suite-wide, `--library` appears only in `src/cli.rs` (7 sites, all `--library`-only or `--clump_only`) and `tests/integration_library.rs`. No test anywhere combines `--library` with `--hardtrim5/3`, `--clock` or `--implicon`. |
| 9 | Specialty modes write no trimming report; exactly three `report::TrimConfig` sites; the gate reaches the text line on all three | **VERIFIED.** `report::TrimConfig {` occurs at `main.rs:1449` (`run_single_file`, `paired: false`), `main.rs:2239` (`run_ubam_output_single`, `paired: false`), `main.rs:2618` (`write_paired_reports`, `paired: true`), plus the `report.rs:1345` test helper. All three set `library: Some(cli.effective_clips())` (`main.rs:1471/2261/2640`), so `config.paired` reaches the gate on every one. |
| — | F-5b's agreement is a genuine invariant, not a coincidence | **VERIFIED, and stronger than the plan claims.** Top-level JSON `clip_r*` come from `config.clip_r*` (`main.rs:1485-1488`), which come from `clips.*` (`main.rs:1223-1226`). The one place they could diverge is the RRBS auto-clip at `main.rs:1136`, which is gated on `clips.clip_r2.is_none()` — and every preset supplies `clip_r2`, so the auto-clip is unreachable whenever a preset is in force. The invariant therefore holds under `--rrbs` too. |
| — | Suite is 751 tests at `145432e` | **VERIFIED** by summing the 17 `test result` lines from `cargo test --no-fail-fast` (§2). |
| — | T-3's premise: the widened guard refuses `--rename --output-format ubam --library emseq` today | **VERIFIED empirically** (§4), and `main.rs:391` is the **only** `--rename`+uBAM guard, so reverting it makes the run succeed — V10 is not vacuous. |
| — | The proposed `flag_summary` passes clippy `-D warnings` | **VERIFIED** — clean under the applied patch. |
| — | V5's `cargo test flag_summary` | **REFUTED** — matches zero tests (§3). |
| — | O-1's "a machine consumer … can read `paired` itself" | **REFUTED (key name).** There is no `paired` key. It is top-level `json["mode"]` = `"single-end"` / `"paired-end"` (`report.rs:803-813`). The argument survives; the key name is wrong. |
| — | V13 (`--rrbs` byte-identity) will pass | **VERIFIED by construction** — no `--library` means `config.library.preset` is `None`, so the `Clipping in force:` line is never written. See O-4: V13 is near-tautological. |

---

## 2. The decisive experiment: the plan's own B2 patch, whole suite

I applied step 1 and step 2 verbatim (3-tuple `flag_summary(paired)`, `clips.flag_summary(cli.paired)`,
`library.flag_summary(config.paired)`) and ran `cargo test --no-fail-fast`:

```
test result: ok. 471 passed; ...   (lib)
test result: ok. 15 passed; ...    test result: ok. 13 passed; ...
test result: ok. 16 passed; ...    test result: ok. 11 passed; ...
test both_reports_record_the_preset_and_its_expanded_values ... FAILED
test result: FAILED. 5 passed; 1 failed; ...   (integration_library)
test result: ok. 1 / 8 / 57 / 14 / 21 / 3 / 57 / 19 / 39 passed ...
```

471+15+13+16+11+6+1+8+57+14+21+3+57+19+39 = **751**. Exactly one failure, and it is step 10's test.
The produced report line is byte-for-byte what step 10 specifies:

```
Library preset: pbat
Clipping in force: --clip_R1 8 --three_prime_clip_R1 8
```

`cargo clippy --all-targets -- -D warnings` was clean on the patched tree (it did recompile:
`Compiling trim-galore v2.3.0` … `Finished`). Patch reverted; see §5.

**This is the strongest thing in the plan's favour: the blast radius is exactly one test, as claimed.**

---

## 3. Mutation checks: are T-1 and F-5a real, and is V5 real?

**T-1 and F-5a: both real, both confirmed blind.** I applied two mutations to `report.rs` at once —
swapped the `"clip_r2"` / `"three_prime_clip_r1"` JSON key labels inside the `library` block (the plan's
V8), and replaced `None => writeln!(w, "{}\"library\": null,", i2)?,` with `None => {}` (the plan's
V11):

```
running 3 tests
test io::tests::test_json_report_name ... ok
test report::tests::test_json_report_records_library_block ... ok
test report::tests::test_json_report_library_is_null_without_a_preset ... ok
test result: ok. 3 passed; 0 failed
```

and the whole integration file too:

```
running 6 tests
test a_run_without_a_preset_reports_no_library_block ... ok
test both_reports_record_the_preset_and_its_expanded_values ... ok
... test result: ok. 6 passed; 0 failed
```

Both mutations survive at unit **and** integration level. The plan's premise for steps 6, 9 and 12 is
correct, and the `contains_key` fix is the right one (`serde_json`'s `Index` returns `&Value::Null` for
a missing key, so `is_null()` cannot distinguish absent from present-and-null).

**V5 is vacuous.** Verbatim:

```
$ cargo test flag_summary
running 0 tests
test result: ok. 0 passed; 0 failed; 0 ignored; 0 measured; 471 filtered out
... (nine more binaries, all "running 0 tests ... ok")
```

Zero tests, exit ok — the repo's own "green over zero" hazard. And it stays zero after implementation:
the plan adds no test whose name contains `flag_summary` (they are
`test_text_report_omits_r2_clipping_on_single_end` and `paired_report_lists_all_four_values`).

**`--implicon 8` does not parse the way step 4 item 4 assumes.** `cli.rs:419-421` is
`#[clap(long = "implicon", alias = "umi_from_r2", default_missing_value = "8", num_args = 0..=1, require_equals = true)]`.
A probe test (inserted, run, reverted) printed:

```
PROBE  implicon=Some(8) input=["8", "test_files/BS-seq_10K_R1.fastq.gz", "test_files/BS-seq_10K_R2.fastq.gz"]
PROBE  validate=Some("Paired-end mode requires an even number of input files (R1/R2 pairs), got 3")
PROBE2 implicon=Some(8) input=["test_files/BS-seq_10K_R1.fastq.gz", "test_files/BS-seq_10K_R2.fastq.gz"]
PROBE2 validate=None
```

`require_equals` makes the bare `8` a positional input, so `validate()` reports an arity error that
contains neither `--library` nor `--implicon`. PROBE2 also independently confirms the premise of the
whole refusal: `--implicon --paired --library pbat R1 R2` is **accepted today** (`validate=None`).

---

## 4. Measured outputs for the docs steps

Under the applied gate, the real single-end override run:

```
$ trim_galore --dont_gzip --library emseq --clip_R1 12 -o $D test_files/BS-seq_10K_R1.fastq.gz
Library preset 'emseq' selected: --clip_R1 12 --three_prime_clip_R1 10
--clip_R1 12 was given on the command line and overrides the emseq preset value 10
```

(3 files in the output dir; no fixture was touched — see §5.) This is the exact replacement for
`library-presets.md:52`, and it is also the evidence for finding I-1.

T-3's premise, at unmodified `145432e`:

```
$ trim_galore --rename --output-format ubam --library emseq -o $D test_files/BS-seq_10K_R1.fastq.gz
Error: --rename is refused with --output-format ubam when any input is FASTQ. ...
=== files=       0
```

Refused today, writes nothing — step 15 is writable exactly as specified.

`cargo fmt --all -- --check` on the patched tree is **not** clean:

```
Diff in .../src/report.rs:276:
-        writeln!(w, "Clipping in force: {}", library.flag_summary(config.paired))?;
+        writeln!(
+            w,
+            "Clipping in force: {}",
+            library.flag_summary(config.paired)
+        )?;
```

---

## 5. Worktree hygiene

Every mutation was reverted and the revert **verified, not asserted**:

```
$ git status --short
$ git rev-parse HEAD
145432e179ae614814dd9754846b76c77137127b
$ git status --short test_files/
```

Both `git status --short` invocations produced no output. `test_files/` untouched (all runs used
`-o $TMPDIR/...`, never `-o .`). `pr441-B` and the `cmp436-*` worktrees were not touched. No `gh`
command was run. `PLAN.md` was not edited.

Re-verified after the report was written, on request:

```
$ git status --short
$ git diff --stat HEAD
$ git rev-parse HEAD
145432e179ae614814dd9754846b76c77137127b
```

Both `git status --short` and `git diff --stat HEAD` produced no output — no tracked modifications.

---

## 6. Findings

### Critical

**C-1 · `CHANGELOG.md:21-22` is falsified by B2 and is not on the plan's fix list.**

#441's own entry, in the same `### Unreleased` / `#### Changes` block this plan's entry will join, says:

```
`--clip_R2` already did. Both trimming reports record the preset name and all four values in
force (`parameters.library` in JSON, `null` when no preset was used), so preset values can
change in a future release ...
```

After B2 that is false for single-end text reports. The plan's §7 states its own rule precisely —
*"Docs changes required for factual correctness by our behaviour change are in scope — that is a
different thing from register editing"* — and then applies it only to `library-presets.md`. The
changelog is the artifact a reader consults to decide whether an upgrade affects them, and here the
two entries sit adjacent: one says all four values are always recorded, the next says Read 2 is no
longer listed on single-end.

This is not the excluded item. The exclusion is *"trimming #441's CHANGELOG entry"* — a
length/register judgement. Correcting a claim your own change falsifies is the §7 case.

*Recommendation:* because `--library` is unreleased, the cleanest fix is to edit #441's entry in place
(`all four values in force` → `the values in force`) and **not** add a second entry describing the
retraction. A user reading the release notes has never seen the un-gated behaviour, so "no longer
lists Read 2 clipping" describes a transition that never existed for them. The specialty refusal
still deserves a clause, which fits inside #441's existing entry ("…and the four specialty modes,
which perform no end clipping, refuse it"). That also removes the awkwardness of the plan's proposed
entry, which is otherwise format-compliant and genuinely one sentence.

**C-2 · The stderr half of B2 is silently revertible — no test covers it.**

The plan names the two call sites as "the load-bearing detail" and adds tests for `report.rs`
(step 7 unit + step 10/13 integration) but **nothing** for `main.rs:1121`. I checked every candidate:
`explicit_clip_flag_overrides_the_preset_and_is_logged` is the only test that inspects stderr, and it
asserts `err.contains("--clip_R1 12") && err.contains("preset value 10")` — both strings come from the
*override* line, not the `selected:` line. So changing `clips.flag_summary(cli.paired)` to
`clips.flag_summary(true)` at `main.rs:1121` leaves all ~763 tests green, and V7 (a manual one-off
comparison) is the only thing standing between that and a merge.

Of the seven findings being fixed, this is **the only one whose fix is silently revertible**. The
other six are all covered: V2 kills the four refusal tests, V3 the ordering test, V6 the report gate
(proven in §2), V8/V9 the T-1/T-2 test, V10 the uBAM test, V11 the two presence assertions, V12 the
agreement loop.

*Recommendation:* one line, in the test that already captures the right stderr. That test is a
single-end `--library emseq --clip_R1 12` run, and emseq's `--clip_R2` never appears in its expected
output, so:

```rust
assert!(!err.contains("--clip_R2"), "single-end log must not list R2 clipping:\n{err}");
```

Add the corresponding positive to step 13's new paired test (`err.contains("--clip_R2 15")`) and both
halves of the "load-bearing detail" are pinned. Then V7 can stay as documentation rather than being
the only guard.

**C-3 · `library_plus_implicon_rejected` (step 4 item 4) tests the wrong thing.**

Proven in §3: `["trim_galore", "--implicon", "8", "--paired", "--library", "pbat", R1, R2]` yields
`input = ["8", R1, R2]` and `validate()` returns *"Paired-end mode requires an even number of input
files (R1/R2 pairs), got 3"*. The plan's assumption 6 hedges (*"verify at implementation time and use
`--implicon=8` if not"*), which leaves a trap: with the plan's specified assertion shape
(`err.contains("--library") && err.contains("--implicon")`) the test fails loudly, but written as a
bare `assert!(cli.validate().is_err())` it would pass **for the arity error**, and the `--implicon`
arm of the refusal would ship untested.

*Recommendation:* drop assumption 6 and specify the proven shape, which is the one the sibling test
at `cli.rs:2343` already uses: `["trim_galore", "--implicon", "--paired", "--library", "pbat", R1, R2]`
— no value at all, relying on `default_missing_value = "8"`. PROBE2 confirms this parses to
`implicon=Some(8)` with the two files as input.

### Important

**I-1 · The stderr line's wording is now inaccurate, and contradicts the docs table two sections up.**

Measured (§4): `Library preset 'emseq' selected: --clip_R1 12 --three_prime_clip_R1 10`. But emseq
*is* 10/10/10/10, and `library-presets.md`'s Presets table says so. "selected" is a statement about
which preset was chosen and what it contains — not about what was applied.

That is precisely the distinction the plan draws in its own B2 edge cases to justify leaving the
override lines ungated: *"That sentence is a true statement about what the user typed; it makes no
claim that the value was applied. `Clipping in force:` does make that claim, which is the whole reason
it is being gated."* By that reasoning the `selected:` line sits on the **ungated** side, yet the plan
gates it. The plan is between two coherent positions:

- **(a)** Gate both, and reword `main.rs:1119` so the line is about application, e.g.
  `"Library preset '{}' selected; clipping in force: {}"`. Keeps V7's "log and report agree" exactly
  as written, and both lines become true. **Recommended** — it is one string, and the log is what a
  user actually reads.
- **(b)** Leave `main.rs:1121` as `flag_summary(true)` so the log shows the preset's full expansion
  (accurate to "selected"), and gate only the report. Then V7 must be rewritten, C-2's test becomes
  the opposite assertion, and the plan's headline "both call sites or they disagree" framing no longer
  applies.

Either is defensible; the plan should pick one deliberately rather than inherit (a) without noticing
the wording drifted. Note also A's F-9: this line prints once per input file, so a 10-file run repeats
the inaccuracy ten times.

**I-2 · The text/JSON asymmetry (assumption 3 / O-1) is defensible, but the plan under-argues it and
under-documents it.** The lead asked for the argument both ways.

*For the asymmetry.* The two artifacts have different readers and different contracts. The text
report is prose a human reads once; `Clipping in force: --clip_R2 8` on a single-end run is a false
statement to that reader, full stop. The JSON is a record a machine reads, it carries
`json["mode"] == "single-end"` at top level (verified, `report.rs:803-813`), and its top-level
`clip_r2` key has echoed inapplicable single-end values since **before** #441 — so gating the JSON
would make the nested `library.clip_r2` disagree with the top-level `clip_r2` on single-end runs,
which is exactly the invariant F-5b's new agreement loop is being added to enforce. **Gating the JSON
would break F-5b.** That is a decisive argument the plan does not make, and it is stronger than the
schema-change argument it does make. Symmetry here would buy consistency between two files at the
cost of consistency *within* one.

*Against.* The same run now produces two different answers to "what clipping was applied", and the
docs page displays them adjacently (lines 62-66 text, 70-79 JSON) with no explanation. A user
diffing a text report against a JSON report — which MultiQC users do — sees a discrepancy and has to
reason about which is authoritative.

*Verdict:* the asymmetry is right, and it is **not** worse than what it replaces (a flat falsehood).
But it must be visible where a reader meets it. See I-3.

**I-3 · Docs step 17 is incomplete in two ways.**

- `library-presets.md:60` — *"Both trimming reports record the preset name **and** the four values
  that were used."* — is itself falsified, not just the 62-66 fence beneath it. Step 17 says only
  "the text-report transcript at `:62-66` … Same fix."
- The explanatory sentence belongs **here**, at 60-79, where the two-value text block and the
  four-value JSON block sit six lines apart. Step 16 puts it in the Overrides section instead, which
  is where the log line lives — a reader who lands on "Reporting and versioning" from the sidebar
  never sees it.

*Recommendation:* rewrite `:60` to something like "Both trimming reports record the preset name and
the values in force; the text report lists Read 2 clipping only for paired runs, while the JSON always
records all four alongside `mode`." That closes I-2's visibility gap and the falsified sentence in one
edit.

**I-4 · V5 is vacuous as written** (proven in §3). Replace `cargo test flag_summary` with the actual
test names. Worth generalising: three other V-rows name a fragment rather than a test
(`cargo test library_plus_` is fine — it will match the four new tests — but `..._omits_r2_...` and
`paired_report_lists_all_four` are abbreviated in the table and one of them, `paired_report_lists_all_four`,
does not match the name step 13 actually gives the test, `paired_report_lists_all_four_values`; that
one happens to still match as a prefix). Given the repo's history, every V-row should name a test that
demonstrably runs.

**I-5 · Step 15 should assert the error text, not merely a non-zero exit.** I confirmed (§4) that the
refusal fires today with the #408 message, and that `main.rs:391` is the only `--rename`+uBAM guard —
so V10 is genuinely non-vacuous as the plan claims. But `assert!(!ok)` alone would silently start
passing for an unrelated reason if any future guard refuses that shape first. Assert on
`--rename` + `ubam` in the message, matching how the `cli.rs` refusal tests are written.

**I-6 · Scope accounting: seven review findings are silently omitted, not decided.**

The plan's "Explicitly out of scope" names three items (CHANGELOG length, docs register, undocumented
aliases). Felix's scope decision was *"substantive findings + test hardening; editorial excluded"*.
Measured against that, these fall through without a disposition:

| id | reviewer | pri | why it is not obviously editorial |
|---|---|---|---|
| **F-5c** | A | Low | *Test hardening* — no `--library` differential on the uBAM or interleaved paths. Step 15 adds a uBAM *refusal* test, which is a different thing. The scope decision includes test hardening. |
| **C-1** | B | **Medium** | Source-comment convention — but B notes one instance is printed to users via `--help`, which makes it user-facing text, not register. Steps 1 and 3 incidentally fix two instances. |
| **L-1** | B | Low | *Reporting behaviour*: a run with explicit `--clip_R1 10` and no preset writes no `Clipping in force:` line at all. B: "worth a decision rather than a silent asymmetry." B2 makes the text report's clipping coverage patchier, so this is now adjacent. |
| **F-6** | A | Low | *Behaviour*: `--library EMSEQ` is refused while `--clip_R1`/`--clip_r1` both work. One `ignore_case = true`. |
| **F-9** | A | Low | *Behaviour*: log line repeats per input file, after adapter detection. Interacts with I-1. |
| F-7 | A | Low | `json_null()` exists, `"library": null` is hand-written. |
| L-2 / L-4 | B | Low | `match … and_then` vs let-chain; a failing `Vec<u8>` compare dumps 39.7 MB. |

I am not proposing any of these — that is Felix's call and O-3 already records his filing bar.
I am saying the plan should carry a one-line disposition for each so the decision is recorded rather
than inferred from absence. F-5c and C-1 are the two where "out" is least self-evident.

**I-7 · `cargo fmt` reformats `report.rs:279`** (diff quoted in §4). Step 2's one-token edit leaves the
tree unformatted and CI's `cargo fmt --all -- --check` fails. One sentence in step 2.

### Optional

**O-1 · O-1's key name is wrong.** There is no `paired` key in the JSON. It is top-level
`json["mode"]`, values `"single-end"` / `"paired-end"` (`report.rs:803-813`). Fix the wording so an
implementer does not go hunting for `parameters.paired`.

**O-2 · B1's description reads as if the new block always fires; under `--clump_only` it never does.**
`--clump_only --hardtrim5 --library emseq` reports `--clump_only does not clip; --library is not
compatible` (`:850`), because the whole `clump_only` block at `:780` precedes both the specialty
mutual-exclusion bails at `:893-903` and the new block. Is that sensible? Marginally: the more
fundamental defect is "two mutually exclusive modes", and the user is told about `--library` instead.
But this is pre-existing #441/`--clump_only` ordering (the same happens today for
`--clump_only --hardtrim5 --clip_R1 5`), it is not introduced here, and fixing it means reordering the
`clump_only` block. Leave it; just add a clause to B1 so the reachability is on record.

**O-3 · "does not clip" is literally true but reads as "does no trimming".** Structurally confirmed:
`specialty.rs` never touches the clip flags, so no preset can reach a clip site in any of the four
modes. But `--clock` and `--implicon` *do* trim — `--clock` extracts an 8 bp UMI + 4 bp fixed sequence
and clips R1 at 13 / R2 at 15. A user told "`--clock` does not clip" may reasonably conclude the mode
performs no trimming. The `--clump_only` family can use that phrasing safely because clump-only really
does nothing to the reads. The #421 family and the #369 `-a2` warnings both chose per-mode accuracy
("performs no adapter trimming", "writes no QC report") for exactly this reason.
`"--clock does not apply end clipping; --library is not compatible"` stays inside the wording family
Felix chose, is exactly true for all four, and costs one word. Felix's call — I am not reopening the
refuse-vs-warn decision, only the accuracy of one phrase.

**O-4 · V13 is close to a tautology and V15 is the expensive one.** Nothing in this change can reach
`--rrbs` output: no `--library` ⇒ `resolve()` early-returns with `preset: None` ⇒ the
`Clipping in force:` line is never written, and `flag_summary`'s signature change is a compile-time
break, not a silent one. V13 will pass; run it if cheap, but the risk lives in V1-V12. V15 ("run the
exact commands in `library-presets.md` and diff against the page") is the one with real value now that
three transcripts change — keep it, and note it needs the `--dont_gzip`-free real invocation the page
shows.

**O-5 · `flags.md:14` does not mention the refusal.** Not falsified by B2 (it says `--library`
"expands into the four clip flags … and nothing else", which stays true). But `flags.md` is the
per-flag reference and now under-documents a hard refusal. One clause, or deliberately out — either
is fine, it just is not currently decided.

---

## 7. Step ordering and executability

**Ordering is coherent.** The only compile-order hazard is step 1 changing `flag_summary`'s signature
before step 2 updates the call sites — and that is a hard compile error, not a silent one, so the
sequence 1 → 2 is right and cannot be got wrong by accident. Steps 3-4 (cli), 5 (report unit tests),
6 (integration), 7 (uBAM), 8 (docs/CHANGELOG) are independent of each other and of 1-2, except that
step 10 must land in the same commit as steps 1-2 or the suite is red between them.

**One sequencing note the plan does not make:** steps 6/9/11/12 change tests whose *current* form is
the thing V8/V11/V12 mutate against. Do the mutation runs (V8-V12) **against the new tests**, not the
old ones; running them against the old tests just reproduces §3's "survives" result. The plan's
V-table says "currently green — B's mutation B", which is right but easy to misread as "verify the old
state".

**Efficiency:** agreed, nothing here is on a hot path. `flag_summary` runs at most a few times per
run; the 3-tuple is stack-allocated and the loop is still four iterations; the `validate()` block
short-circuits on `self.library.is_some()`. No further comment.

---

## 8. Alternatives

**A · Store `paired` inside `ResolvedClips` at construction.** `Cli::effective_clips()` already has
`self.paired`, so `resolve()` could take it and `flag_summary()` could keep its zero-arg signature.
*Pro:* it deletes the entire two-call-site desync hazard the plan calls its load-bearing detail — one
source of truth, set once — and with it C-2 (there would be no call-site argument to get wrong).
*Con:* `report::TrimConfig` **already carries `paired`** (`report.rs:172`), so the field would exist
twice in the same config, and the two could themselves diverge — note `main.rs:2239` hard-codes
`paired: false` for the uBAM single path while `clips.paired` would come from `cli.paired`. It also
makes `resolve()` take a parameter it does not resolve anything with, and touches all five
`library.rs` unit tests. *Verdict:* the plan's threading is the better fit, precisely because
`TrimConfig.paired` exists — but the plan should record that it considered and rejected this, and it
must close C-2, which is the cost of choosing threading.

**B · `flag_summary` returns the parts; callers filter.** Rejected: it duplicates both the filter and
the join at two call sites, which is strictly worse for the desync risk being addressed.

**C · Hoist the preset announcement into `main()` before dispatch** (A's F-9 suggestion, which she
notes "would incidentally give F-1 a natural home"). Settled against — `Cli::validate()` is the right
home because it provably precedes `ensure_output_dir` and every dispatch arm, which is what makes V4
true. Noting it only so the record shows the alternative was on the table.

---

## 9. Assumptions the plan left implicit

Beyond its seven stated ones:

1. **`overrides` order is a stable contract.** Step 6 asserts `lib["overrides"][0..3]` by index. That
   order comes from Rust evaluating struct-literal fields in source order inside `resolve()` — real
   but unstated. Worth one line in the test's doc comment, since a future reorder of the four
   `pick` calls would break the test for a reason that looks unrelated.
2. **`serde_json`'s `Index` returns `Null` for missing keys.** The whole F-5a fix rests on it. Stated
   as a fact in step 9; it is correct (I confirmed by mutation) but it is an assumption about a
   dependency, not about this code.
3. **No test currently asserts the stderr `selected:` line.** The plan's V7 implicitly assumes manual
   verification is adequate here. Made explicit as C-2.
4. **`--library` is unreleased.** Load-bearing for the "not a breaking change" argument *and* for
   C-1's recommendation; stated in §Context but not listed as an assumption, and it is the one that
   would invalidate the delivery plan if the release un-parked.
5. **`cargo fmt` will not reformat anything else.** Verified for `report.rs:279`; not checked for the
   new `cli.rs` block, which should be fine at its indentation but is worth a `cargo fmt` run before
   the first commit.

---

## 10. Action items

### Critical (do before implementing)

1. **C-1** — Correct `CHANGELOG.md:21-22` ("all four values in force"), which B2 falsifies. Recommend
   editing #441's entry in place and folding the refusal clause into it, rather than adding a second
   entry that retracts an unreleased behaviour.
2. **C-2** — Add a stderr assertion so the `main.rs:1121` half of B2 is not silently revertible:
   `assert!(!err.contains("--clip_R2"))` in `explicit_clip_flag_overrides_the_preset_and_is_logged`,
   plus the paired positive in step 13.
3. **C-3** — Replace step 4 item 4's `["--implicon", "8", …]` with the proven
   `["--implicon", "--paired", "--library", "pbat", R1, R2]`; delete assumption 6.

### Important

4. **I-1** — Decide the stderr wording. Recommend rewording `main.rs:1119` to be about application so
   V7 and the gate both stay honest.
5. **I-3** — Rewrite `library-presets.md:60` (falsified prose) and site the text/JSON explanation at
   60-79, not only in Overrides.
6. **I-4** — Make every V-row name a test that demonstrably runs; `cargo test flag_summary` matches
   zero.
7. **I-5** — Step 15 asserts the error text, not just exit status.
8. **I-6** — Add a one-line disposition for F-5c, C-1, L-1, F-6, F-9, F-7, L-2, L-4.
9. **I-7** — Note in step 2 that `report.rs:279` must be reformatted (rustfmt splits it).
10. **I-2** — Add the decisive argument for the asymmetry: gating the JSON would break the F-5b
    agreement invariant the plan is simultaneously introducing.

### Optional

11. **O-1** — `json["mode"]`, not `paired`.
12. **O-2** — One clause in B1 noting `--clump_only` short-circuits before the new block.
13. **O-3** — Consider "does not apply end clipping" for factual accuracy under `--clock`/`--implicon`.
14. **O-4** — V13 is near-tautological; V15 is where the value is now.
15. **O-5** — Decide whether `flags.md:14` gets a refusal clause.
16. **§8-A** — Record why `paired`-in-`ResolvedClips` was rejected.
17. **§7** — Note that V8-V12 must be run against the *new* tests.
18. **§9** — Promote the five implicit assumptions.

---

## 11. What I verified as correct (so absence of a finding means absence of a defect, not of a check)

- All 14 line-number anchors resolve to what the plan says. Eleven are exact; two are off by one or
  two lines at the end of a range; none is wrong about what it points at.
- `flag_summary` has exactly two call sites. No third, none in tests.
- The refusal's placement is genuinely last among fallible checks; all four claimed precedences hold.
- No existing test combines `--library` with a specialty mode — the refusal breaks nothing.
- Three `report::TrimConfig` sites, all carrying `paired` and `library`; the gate reaches all three.
- The four-distinct-values lever works exactly as reasoned, including deterministic `overrides` order.
- T-1 and F-5a are real, proven by surviving mutations at unit *and* integration level.
- T-3's premise is real: the refusal fires today, `main.rs:391` is the only guard, V10 is not vacuous.
- F-5b is a genuine invariant, and holds under `--rrbs` for a reason the plan does not state.
- Specialty modes cannot clip: `specialty.rs` has zero clip references.
- The plan's B2 patch compiles, passes clippy `-D warnings`, and breaks exactly one test out of 751 —
  the one predicted, with the exact predicted output.
- The proposed CHANGELOG entry is one sentence and matches the file's format.
- The proposed source comments comply with CLAUDE.md (one line, current-state fact, no justification);
  step 1's instruction to delete the two justification sentences from `flag_summary`'s doc comment is
  correct and incidentally closes two instances of B's C-1.

**Not checked:** V13's differential build (traced as unreachable rather than run); V15's docs
transcripts beyond the two I reproduced; behaviour on Linux (the ±2 test-count difference the plan
warns about).
