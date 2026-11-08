# Plan Coverage Report

**Mode:** B (adjusted — no `IMPL.md`; `PLAN.md` §5 + §9 serve as the spec)
**Plan:** `plans/11062026_ubam-phred64-preflight/PLAN.md`
**Date:** 2026-11-08
**Tree:** `dev` @ `49152c1` + uncommitted changes (5 files)
**Verdict:** COMPLETE — 7/7 implementation steps DONE, 7/7 validation rows DONE, 0 MISSING, 0 PARTIAL

## Summary

- Total items: 20 (7 × §5 steps, 7 × §9 validation rows, 3 × §9 mutations, 3 × §9 suite gates)
- DONE: 19
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 1 (§9 row 4 — documented in §13 Deviation 1, and the chosen form is one the plan itself offered)
- NOT RE-EXECUTED: 3 mutations (audited by inspection only — see the mutation section)

Both §13 deviations were checked against the code rather than taken at their word; both descriptions
are accurate. The §13 "all three mutations killed their predicted targets" claim was **not**
re-executed (shared tree, no source mutation permitted) but every claimed mutation→test pairing is
sound by inspection, and the two "exactly N failed / M passed" tallies are internally consistent
with the 18-test file.

## Coverage ledger — §5 Implementation outline

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | Widen the pre-flight condition to `fastqc_requested() \|\| OutputFormat::UBam` | §5.1 | DONE | `src/main.rs:388-391`, verbatim the plan's snippet |
| 2 | Two-line comment stating the current rule, no justification | §5.2 | DONE | `src/main.rs:386-387`, verbatim the plan's example; no mention of the old gate, FASTQ behaviour or reasoning |
| 3 | Disambiguate the startup NOTE; update the `--fastqc` sentence; keep the quality-trimming clause | §5.3 | DONE | `src/main.rs:558-570`; all three sub-requirements met |
| 4 | Refresh the stale test comment and rename the test | §5.4 | DONE | renamed at `tests/integration_fastqc_phred64.rs:304`, comment rewritten at `:308-309`; no stale references to the old name remain anywhere outside `plans/` |
| 5 | Add the §9 tests, all refusal tests using `assert_preflight_refused` | §5.5 | DONE | 4 new tests at `tests/integration_fastqc_phred64.rs:548-651`; all 3 refusal tests call `assert_preflight_refused` |
| 6 | CHANGELOG: new entry **and** amend the pre-existing #532 entry | §5.6 | DONE | `CHANGELOG.md:29-33` (new entry) and `:25-26` (the amended #532 sentence). **Both** halves landed |
| 7 | Docs — **two** files | §5.7 | DONE | `docs/src/content/docs/guide/outputs.md:101` and `docs/src/content/docs/guide/quality.md:37,40` |

### Detail on the items worth confirming

**§5.1 / §5.2 — placement and comment.** The condition sits at `src/main.rs:388`, between the
uBAM-input guard (`:369`) and the invariant block, which is before `naming::ensure_output_dir` (now
`:551`, shifted from the plan's `:546` by the 5 added lines). The required ordering holds. The new comment is exactly the two lines the plan drafted:

```rust
// On uBAM output the declared offset feeds the QUAL arithmetic, and bases
// below Q31 floor at 0 rather than shifting.
```

It is appended below the pre-existing 5-line `#532` comment rather than replacing it; §5.2's
"at most two lines" governs the new text, which it satisfies. The pre-existing lines were left
intact, which is correct — they still describe the `--fastqc` half of the gate.

**§5.3 — NOTE disambiguation.** All three of §5.3's requirements are met independently:

1. *Disambiguated.* "Writing FASTQ, the modes that do no quality arithmetic — `--hardtrim5/3` and
   `--clump_only` — pass the scores through untouched; writing uBAM they do not." Fronting the
   qualifier makes it govern both list items unambiguously, and the trailing clause closes the
   reading reviewer A observed misfiring.
2. *Gate sentence updated.* "A mislabelled input is refused before anything is written when the
   output is uBAM, or when `--fastqc` is given."
3. *Quality-trimming clause kept.* "Wherever quality scores are read (quality trimming, and BAM
   QUAL on uBAM output)…" is unchanged. §2.4/§10 required this explicitly — it is the only warning
   covering the FASTQ `-q` band the plan leaves open. **Not shortened.**

**§5.6 — the CHANGELOG amendment, not just the new entry.** The lead flagged this as a thing to
confirm rather than assume. Both landed:

- *New entry* (`CHANGELOG.md:29-33`), under `### Version 2.4.0`, directly after the #532 entry:
  bold lede + one sentence of substance pointing at `docs/guide/quality.md` rather than inlining
  migration guidance. Matches §5.6's resolution of the "one sentence" conflict.
- *Amendment* to the pre-existing #532 entry: `--phred64 without --fastqc is unaffected` →
  `on FASTQ output, --phred64 without --fastqc is unaffected`. This is the exact sentence §5.6
  identified as having become false.

**§5.7 — both docs files.** The lead flagged this as "two files, not one". Both are present:

- `outputs.md:101` — added to the canonical "Rejected at startup, before anything is written (v1
  scope)" list, in the conditional shape §5.7 named `--rename` as the precedent for ("only when the
  input contradicts the declaration … Genuine Phred+64 input is accepted").
- `quality.md` — the `--output-format ubam` paragraph (`:37`) gains "a mislabelled input is refused
  before anything is written", and the `:::caution` (`:40`) widens "With `--fastqc`" to "With
  `--fastqc`, or whenever the output is uBAM". Register is plain and consistent with the
  surrounding prose.

The cross-link `/guide/quality/#phred-encoding` added in `outputs.md` resolves: `## Phred encoding`
is at `quality.md:29`, and the trailing-slash form matches `guide/overview.md:18` and
`guide/flags.md:49`. §13's verification claim here is accurate.

## Coverage ledger — §9 Validation

| # | Verify | Test | Status | Result |
|---|---|---|---|---|
| 1 | Mislabelled input refused before any write | `preflight_fires_on_ubam_output_without_fastqc` | DONE | PASS |
| 2 | Silent hardtrim path closed | `preflight_fires_on_hardtrim_ubam_output` | DONE | PASS |
| 3 | Silent `--clump_only` path closed, **no `--length`, no `-q`** | `preflight_fires_on_clump_only_ubam_output` | DONE | PASS |
| 4 | NOTE no longer claims hardtrim is safe | `note_qualifies_pass_through_by_output_format` | DEVIATED (within plan intent) | PASS |
| 5 | Error names file, character, ASCII value | asserted inside row 1 + `assert_preflight_refused` | DONE | PASS |
| 6 | #362 regression — four existing guards stay green | see below | DONE | 4/4 PASS |
| 7 | FASTQ output unaffected | `preflight_does_not_fire_without_fastqc_on_fastq_output` | DONE | PASS |

**Row 3 — the C-2 trap is genuinely avoided.** The test passes `--phred64 --clump_only
--output-format ubam` and nothing else: no `--length` (which `cli.rs:826` rejects outright, the
unrelated-error trap), and no `-q`. It asserts through `assert_preflight_refused`, which checks the
refusal's *text* (`"--phred64 declares Phred+64 input"`) and that the output dir was never created.
A non-zero exit from the `--length` rejection could not satisfy it. The trap §9 warned about is
closed, and the test carries a comment saying why both flags are absent.

**Row 5 — the message really does carry all three fields.** `src/main.rs:410-420` interpolates
`input.display()`, `lowest.escape_ascii()` and `lowest`. The row-1 test asserts `ASCII 35`, which is
`#` — the character `write_phred33_fixture` (`:243`) writes for every base. The path and the
escaped character are covered by the surrounding format string, which `assert_preflight_refused`
anchors on.

**Row 6 — the four named guards exist at the named locations and pass.** The lead asked for this to
be confirmed rather than assumed:

| Test | Named location | Found at | Result |
|---|---|---|---|
| `phred64_ubam_out_se_stores_true_phred` | `integration_ubam_out.rs:632` | `:632` | PASS |
| `phred64_ubam_out_specialty_stores_true_phred` | `integration_ubam_out.rs:666` | `:666` | PASS |
| `phred64_ubam_out_pe_two_file_stores_true_phred` | `integration_ubam_out.rs:799` | `:799` | PASS |
| `phred64_clump_only_ubam_out_stores_true_phred` | `integration_clump_only_ubam.rs:779` | `:779` | PASS |

All four line numbers are exact. The `quals()` helper §9 cited is at `integration_ubam_out.rs:616`,
also exact. §9's "Baseline: 75 tests green across the three files" checks out: the three files now
hold 44 + 17 + 18 = 79 `#[test]` functions, i.e. 75 before the 4 added here.

**Row 4 — DEVIATED, correctly documented.** §9 row 4 offered two forms: "assert on stderr of test 2
before the refusal, **or** on a FASTQ-output hardtrim run". The implementation took the second.
§13 Deviation 1 records this and gives the reason, and the reason is correct: the NOTE is printed
at `src/main.rs:558-570`, which is *after* `ensure_output_dir` at `:551` and therefore after the
pre-flight at `:388`, so a refused run emits no NOTE at all and the first form is not
implementable. The test asserts the run succeeds (FASTQ output) and that stderr contains
`"writing uBAM they do not"` — i.e. it pins the disambiguation §5.3 added, which is what row 4 is
for. Marked DEVIATED only because a form was chosen; it is within plan intent, not a gap.

## Mutation table — audited by inspection, NOT re-executed

**I did not re-run these.** Two other agents share this tree and the brief forbids mutating source,
so what follows is an inspection audit of whether each claimed mutation→test pairing is *sound*,
not a reproduction of §13's results.

| Mutation | §13 claim | Sound by inspection? |
|---|---|---|
| drop the `\|\| matches!(… UBam)` arm | rows 1–3 fail; 15 pass | **Yes.** The condition reverts to `phred64 && fastqc_requested()`. Rows 1–3 all pass `--phred64` with uBAM output and no `--fastqc`, so the scan is skipped and each run completes; every one of the three asserts `!out.status.success()`, so all three fail. Nothing else in the file changes behaviour. 3 + 15 = 18 = the file's test count. |
| `\|\|` → `&&` | rows 1–3 **plus 5 pre-existing `--fastqc` tests** fail | **Yes.** The condition becomes `phred64 && fastqc_requested() && UBam`. Rows 1–3 fail as above. Additionally every pre-existing refusal test that pairs `--fastqc` with FASTQ output now skips the scan and succeeds where it asserts refusal — so the extra casualties are expected, and §13 is right that this is a guard the plan had not claimed. |
| `lowest < 64` → `lowest <= 64` | `preflight_accepts_ascii_64_which_is_phred64_q0` fails; 17 pass | **Yes, and the §9 correction is confirmed.** That test's fixture is `write_mixed_fixture(…, Some(0), '@')` → 24×`h` then 10×`@`, minimum byte 64; `64 <= 64` bails, and the test asserts success, so it fails. The draft's original target, `preflight_accepts_genuine_phred64_input`, uses `test_files/phred64_test.fastq`, whose quality alphabet I verified is exactly `{B, h}` → minimum 66; `66 <= 64` is false, so that test would survive the mutant. §9's reassignment is arithmetically correct. 1 + 17 = 18. |

Both tallies in §13 are internally consistent with the file's 18 tests, which is weak corroboration
that the runs happened as described — but it is corroboration, not verification.

The current `src/main.rs` holds `lowest < 64` (`:401`), so no mutant survived into the tree, and the
`cargo test --release` run below linked a fresh binary (51s) in which rows 1–3 pass — which they
cannot under either of the first two mutants. §9's "rebuild after restoring" caveat is satisfied.

## Suite gates (§9 "Suite")

| Gate | Command | Result |
|---|---|---|
| Unit + integration tests | `DEVELOPER_DIR=… cargo test --release` | **842 passed, 0 failed, 0 ignored; exit status 0** |
| Formatting | `DEVELOPER_DIR=… cargo fmt --all -- --check` | exit 0, no diff |
| Clippy (cold) | `CARGO_TARGET_DIR=<scratch> … cargo clippy --all-targets --release -- -D warnings` | **exit 0; 45 crates compiled + 148 checked from scratch; 0 warnings, 0 errors** |

**842 confirmed.** Summing the 21 `test result:` lines gives exactly 842 passed / 0 failed across 20
test binaries plus doc-tests, and the explicitly captured exit status is `0`. §13's figure is
accurate. The four new tests and all four row-6 guards appear individually as `... ok`.

**On clippy.** The first warm `cargo clippy` returned exit 0 in 0.49 s — that is a cache hit, not a
check, so it is not evidence of anything. A genuinely cold run was then issued into a scratch
`CARGO_TARGET_DIR`, leaving the shared `target/` untouched: 45 crates compiled and 148 checked from
scratch in 22.76 s, exit 0, zero warning and zero error lines. §13's cold-clippy claim is
**confirmed independently** (§13 reported 29 crates checked; the scratch run checked a superset, so
the numbers differ only in how much cache each started from). §9's own caveat still stands
regardless: local rustc 1.95.0 vs CI's 1.99 means a clean local clippy is not a clean CI clippy.

## §13 "Not done" items — verified as stated, not gaps

Both are explicit, reasoned non-deliverables rather than dropped work, and neither corresponds to
anything §5 or §9 required:

1. **Docs site not built.** Correct as a judgement call; the link target was independently
   re-verified here (`## Phred encoding` at `quality.md:29`, trailing-slash form matching
   `overview.md:18` and `flags.md:49`).
2. **No test pins the uniform-high-quality residual.** The stated reason is sound — such a test
   would encode a limitation as desired behaviour. The limitation itself *is* documented, in
   `quality.md`'s caution ("The check settles only what the bytes prove…"), so it is not
   undocumented.

The §13 "Finding" (the check cannot refuse a uniformly-high-quality Phred+33 file) is a genuine
limitation of the approach rather than an implementation defect, and the two documentation changes
it drove — the CHANGELOG's "an input holding any quality character below ASCII 64 is refused"
instead of "a mislabelled Phred+33 input is refused", and the new caution sentence — are both
present in the diff and both correctly scoped. `outputs.md`'s entry is indeed already phrased
conditionally ("a file that does is refused"), so §13's claim that it needed no change is accurate.

## Verdict

**COMPLETE.** Every step in §5.1–§5.7 landed, including both halves of §5.6 (new entry *and* the
#532 amendment) and both files in §5.7. Every row in §9's validation table has a test that actually
exercises it, and all of them pass, as do the four pre-existing #362 guards in row 6.

All three suite gates were re-run here and all three are clean: 842 tests pass with exit status 0,
`cargo fmt --all -- --check` exits 0, and a genuinely cold clippy (scratch target dir, 45 compiled +
148 checked) exits 0 with no warnings.

One qualification on the evidence, stated plainly rather than folded into the verdict:

- **The three mutations were not re-executed.** All three claimed pairings are sound by inspection
  and both tallies are internally consistent with the file's 18 tests, but §13's mutation results
  rest on the implementer's own run, not on anything reproduced here. The brief ruled out source
  mutation because two other agents share this tree.

Nothing remains for the implementer to address.
