# Plan Review A — `--passthrough` alias check wording (#389)

**Reviewer:** A (independent, fresh context)
**Plan under review:** `plans/08082026_passthrough-alias-wording/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev`
**Date:** 2026-08-08

**Scope note:** the maintainer's decision — keep the case-folded `collision_key`, change wording only — is taken as given and is NOT relitigated here. This review asks whether the plan, *given that decision*, is correct, complete, and truthfully worded.

**Verdict: REVISE before implementing.** The direction is right and most of the plan's factual claims check out against the code, but one instruction (Implementation step 3, the `norm_path()` → `collision_key` substitution) rests on a false premise and, if executed literally, would replace a *true* comment with a *false* one — the exact defect class #389 exists to remove. Details in §1.1.

---

## 1. Logic review

### 1.1 CRITICAL — `norm_path()` is not a stale name, and the proposed substitution introduces a falsehood

The plan asserts twice that `norm_path()` is dead nomenclature:

> Context: "references **`norm_path()`**, a stale pre-D13 name for what is now `collision_key`."
> Implementation step 3: "replace the stale `norm_path()` reference with `collision_key`".

Both are wrong. `norm_path` still exists, is still `pub`, and is a *component* of `collision_key`, not a former name for it:

- `src/io.rs:48-50` — `pub fn norm_path(p: &Path) -> String { p.to_string_lossy().to_ascii_lowercase() }`
- `src/io.rs:84-86` — `pub fn collision_key(p: &Path) -> String { norm_path(&lexical_normalise(p)) }`

So the existing test comment at `src/cli.rs:1874` — "norm_path() is a single `to_ascii_lowercase()` call which is trivially correct" — is **factually accurate as written**. Performing the substitution the plan asks for yields:

> "collision_key() is a single `to_ascii_lowercase()` call which is trivially correct"

which is false: `collision_key` is absolutise (`std::path::absolute`) + `..`-fold against the component stack + lowercase, three operations across `lexical_normalise` and `norm_path`. A wording-truthfulness fix that ships a new untrue comment 900 lines below the one it repaired is a self-inflicted repeat of #389.

**Recommended fix.** Don't rename inside the sentence; re-anchor the sentence. The property the comment is trying to justify — that case-folded equality holds for case-variant spellings — is *already pinned by a real test*: `io::tests::identity_key_is_case_sensitive_but_collision_key_is_not` (`src/io.rs:1181-1194`) asserts `collision_key("Sample_R1.fastq.gz") == collision_key("SAMPLE_R1.fastq.gz")`. A truthful, shorter comment names that anchor instead of re-deriving the argument:

```rust
// Byte-equal spellings are a subset of case-folded equality; the fold itself is
// pinned by io::tests::identity_key_is_case_sensitive_but_collision_key_is_not.
```

This also retires the plan's framing of cross-case coverage as an untested gap (see §4.2).

### 1.2 IMPORTANT — the proposed comment breaks the repo's own comment rule

The drafted comment (Behavior 2) is four lines of reasoning chain:

```rust
// 1.ix — case-folded on purpose (#389), NOT an identity check: a case-variant
// passthrough IS R1/R2 on APFS/NTFS and would be consumed twice; on a
// case-sensitive filesystem the (pathological) two-real-files shape is
// refused loudly rather than risking silent dual-consume elsewhere.
```

The standing instruction for committed code in this project is: default to one line, two maximum with a reason, state the fact rather than the evidence, no reasoning chains — rationale belongs in the commit message or the linked issue. Four lines of derivation is what that rule exists to prevent, and shipping it inside a *wording-cleanup* diff invites the reviewer note the rule was written to pre-empt. The existing comment being five lines is not a defence; this change is the opportunity to fix that too.

Also, "risking silent dual-consume **elsewhere**" is ambiguous as written — "elsewhere" reads as *another code location*, when the intended meaning is *under the alternative key*. Shortening removes the problem.

Suggested two-liner, fact-stated, issue-anchored:

```rust
// 1.ix — case-folded, not an identity check: on APFS/NTFS a case-variant
// passthrough IS R1/R2 and would be consumed twice (#389). len == 2 per 1.ii.
```

### 1.3 IMPORTANT — the new comment silently drops the `1.ii` invariant note

The current comment carries "self.input.len() == 2 here per 1.ii" (`src/cli.rs:934`). That sentence is the *only* thing in the file explaining why the guard at `src/cli.rs:938` — `if self.input.len() == 2` — is not a bug: 1.ii at `src/cli.rs:896-902` already bailed when `len != 2`, so the condition is provably always true. The plan's replacement comment omits it. A future reader then meets an unexplained redundant branch, and the likely next move is either a puzzled review comment or a "simplification" that has to re-derive the invariant.

Keep the clause (it fits in the two-line budget, as above), or take the Optional route in §5.4.

### 1.4 IMPORTANT — the message does not say *which* input matched, though the code knows

The drafted message interpolates only `pt.display()`. The user is told "matches R1 or R2" and left to work out which. The information is free — the code has already evaluated the two comparisons at `src/cli.rs:940-941`, so it can name the matched input and its on-disk spelling.

This matters most in exactly the subcase the reword is about. `--passthrough r1.fq.gz` with inputs `R1.fq.gz R2.fq.gz` currently produces a message naming only `r1.fq.gz`; naming the counterpart (`… matches R1 (R1.fq.gz)`) makes the case-only nature of the match visible at a glance, which is the whole diagnostic point.

There is direct precedent: PR #219 changed the sibling output-collision error specifically to name both colliding paths, and `.github/workflows/ci.yml:635-637` records that history. `io::preflight_output_collisions` (`src/io.rs:120-126`) names both. The passthrough message is the odd one out, and fixing it is inside "wording only" — no key change, no change to the rejection surface.

### 1.5 IMPORTANT — the remediation advice is tuned to the rarer subcase and misleads on the commoner one

Drafted remediation: "If these are genuinely two different files on a case-sensitive filesystem, rename one so the paths differ by more than letter case."

The plan's Behavior 3 claims this "reads truthfully for both subcases", and the Self-Review calls the byte-equal reading "vacuous-but-true". Truthful, yes. *Useful*, no — and for the byte-equal case it is actively misdirecting.

`--passthrough` exists to carry a third stream (the index / cell-barcode read; see the `BS-seq_10K_I1.fastq.gz` fixture and `tests/integration_passthrough.rs`). The overwhelmingly likely real-world trigger for 1.ix is therefore a byte-equal slip — a script or a hand-typed command that passes R1 where I1 was meant. For that user, the only actionable sentence in the message is "rename one so the paths differ by more than letter case", and following it literally *works*: renaming `R1.fq.gz` to `Q1.fq.gz` and re-passing it clears validation and produces a run whose passthrough output is a duplicate of R1. The check's advice has walked the user around the guard rather than to the fix.

The correct byte-equal remediation is "the passthrough stream must be a third file (e.g. the index read), not one of the two reads being trimmed." Cover both subcases, and lead with the third-file point since it is the common one. §5.3 sketches the two-tier form that gets each subcase its own correct sentence at zero behavioural cost.

### 1.6 IMPORTANT — scope is incomplete: two other comments mis-frame this same check

The plan touches one comment and one message. Two further sites describe 1.ix in the terms #389 objects to, and neither is in the plan:

1. **`src/cli.rs:885-889`** — the nine-item envelope preamble: "case-folded collision check (1.ix) uses `crate::io::collision_key` to share the same APFS/NTFS-aware normalisation as the output-collision pre-flight in main.rs (**issue #216 protection**)." #216 is the *output*-collision issue. 1.ix guards *input* dual-consume — a different failure with a different blast radius. Labelling it "#216 protection" and framing it as sharing the output pre-flight's purpose is the precise conflation this plan is fixing forty-five lines below. Leaving it means a reader who starts at the preamble arrives at the corrected comment already holding the wrong model.

2. **`src/io.rs:38-47`** — `norm_path`'s doc comment: "Used by: `Cli::validate()` to catch `--passthrough` **aliasing** R1 or R2 (e.g. `--passthrough r1.fq.gz` while R1 is `R1.fq.gz`)". Two problems. It is stale as a direct-use claim — `Cli::validate` calls `collision_key`, and `norm_path` has no callers anywhere outside `collision_key` and its own tests (verified by repo-wide grep). And it uses the same identity verb ("aliasing") that #389 asks to be dropped from the user-facing string. This is the doc a reader lands on when chasing `collision_key` → `norm_path`, so it belongs in the same diff.

Both are one-line edits. Including them makes the fix complete; excluding them means #389's "either switch the key or say what it actually checks" is only two-thirds done.

### 1.7 Plan claims that check out

Stated explicitly so the implementer does not redo the verification:

- **Line numbers are accurate.** 1.ix comment `src/cli.rs:933-937`; guard and body `938-949`; message `943-947`; the two pinning tests `1870-1888` with the shared comment at `1872-1877`. All as the plan describes.
- **Check ordering is as claimed.** `check_restartable_input` (1.viii) at `src/cli.rs:932` precedes 1.ix, so the passthrough file always exists when 1.ix fires. The plan's edge-case reasoning follows correctly: on a case-sensitive filesystem a *non-existent* case-variant yields "--passthrough file not found" from 1.viii, so 1.ix's false refusal fires only when both files genuinely exist — the pathological shape the decision knowingly accepts.
- **Both platform subcases behave as the plan describes.** On APFS, `--passthrough r1.fq.gz` + input `R1.fq.gz`: 1.viii opens the same inode, 1.ix refuses — guard working. On ext4 with two real files: 1.viii passes, 1.ix refuses — the false refusal. Confirmed by reading the code path, not inferred.
- **`..` / `./` spellings keep resolving** through `lexical_normalise` inside `collision_key`; `io::tests::both_keys_normalise_spelling` (`src/io.rs:1198-1209`) pins it.
- **CHANGELOG has a home for the entry.** Unreleased → `#### Changes` at `CHANGELOG.md:129`. Precedent for a message-only entry exists in that very section (`CHANGELOG.md:145-146`, the `--output-dir` → `--output_dir` message fix), so step 4 is well-founded rather than over-reporting.

---

## 2. Assumptions

### 2.1 A1 ("no test outside `cli.rs:1871-1888` pins the old text") — CONFIRMED

Independently re-verified by repo-wide grep excluding `target/`, `plans/`, `.git/`:

- `"cannot point at R1"` → three hits: `src/cli.rs:944` (the message), `src/cli.rs:1880`, `src/cli.rs:1887` (the two asserts).
- `"aliases an input"` → one hit: `src/cli.rs:945`.
- `tests/integration_passthrough.rs` asserts on `"--passthrough requires --paired"` (line 200) only — 1.i, not 1.ix.
- No `.github/workflows/*.yml` mentions `passthrough` at all.
- Docs and README carry no copy of this message.

One near-miss worth recording so nobody re-flags it: `.github/workflows/ci.yml:638` does grep `"case-insensitive, for APFS/NTFS safety"`, and `ci.yml:635` names the string `"case-insensitive match"` in a comment. That grep targets `io::preflight_output_collisions`' *output*-collision message, reached by a `--paired` four-file run with no `--passthrough`. It cannot be tripped by this change. The plan's Integration claim ("the validation CI matrix does not exercise `--passthrough` message text") therefore holds — but it holds narrowly, and the reason it holds is worth one line in the plan rather than a bare assertion, because the qualifier string is shared vocabulary between the two messages.

### 2.2 A2 (D13 taxonomy unchanged) — CONFIRMED, and consistent with the code

`path_identity_key` (`src/io.rs:78-80`) and `collision_key` (`src/io.rs:84-86`) are both purely lexical; neither touches the filesystem, and the doc comment on `lexical_normalise` (`src/io.rs:52-54`) states the symlink consequence explicitly. No third key and no platform probe is introduced. Consistent with the plan.

### 2.3 The "no behaviour change" claim — true as intended, imprecisely stated

Same key, same comparisons, same position in the chain, same `bail!` → identical accept/reject sets and identical exit codes. The claim is sound in the sense that matters.

But stderr text *is* observable behaviour to anything grepping it (wrapper scripts, pipeline log assertions). The plan states "**No behaviour change**" flatly, in bold, twice. Recommend restating as: *no change to which runs are accepted or rejected; the message text changes.* That phrasing is also what the CHANGELOG entry should say, so a user who greps for the old string learns why it vanished. This is a precision point, not a disagreement — and the project already treats message rewords as normal changelog-worthy changes (§1.7).

### 2.4 Unstated assumption: "the byte-equal subcase is the portable one, therefore the representative one"

The existing test comment and the plan both lean on byte-equal being the testable subcase. That is correct as a *testability* argument (§4.2). The plan then quietly carries it into the *wording* decision, drafting remediation for the case-variant subcase while calling the byte-equal reading "vacuous-but-true". The frequency ordering is the other way round (§1.5): byte-equal is both the testable subcase *and* the likely one. Worth making explicit, because it is the assumption that produces the message's weakest sentence.

### 2.5 Unstated assumption: `pt.display()` is adequate identification

`display()` is lossy for non-UTF-8 paths, and the message shows the path *as spelled*, not as normalised — so a `../x` spelling is echoed back unresolved while the match happened on the absolutised form. Both are pre-existing, both are consistent with the rest of the codebase's diagnostics ("keeps naming the paths as the user spelled them", `CHANGELOG.md:109-110`), and neither needs changing. Noted only so it is a recorded decision rather than an oversight.

---

## 3. Efficiency analysis

Nil, and the plan's one-line dismissal is the right amount of attention. For completeness: 1.ix performs three `collision_key` calls (one for the passthrough, two for the inputs) once per run at argument-validation time, each `std::path::absolute` + a component walk + one `to_ascii_lowercase` allocation over a path-length string. Adding a second interpolation to the message (§1.4) or a second sentence (§1.5) costs nothing measurable — it is on the error path, which formats once and then the process exits.

Nothing in the plan touches an allocation, a loop, or an I/O boundary. No scalability surface. Confirmed there is nothing missed here.

---

## 4. Validation sufficiency

### 4.1 The pin-update risk is adequately covered

Validations 1 and 3 together close the failure mode the plan worried about. A partial pin update — message changed, one assert missed — is caught twice over: the missed assert leaves a `"cannot point at R1"` hit so Validation 3 is non-zero, *and* `cargo test` fails that test. The plan's Self-Review note about adding Validation 3 for exactly this reason is well-judged.

### 4.2 The "cross-case testing is impossible" claim is right about `validate()` and wrong about the property

The retained test comment says true cross-case testing needs a case-insensitive filesystem CI does not guarantee. For the **full `Cli::validate` path** that is correct, and I confirmed the mechanism: on a case-sensitive filesystem a case-variant fixture path does not exist, so 1.viii (`src/cli.rs:932`) fires "--passthrough file not found" and 1.ix is never reached. You cannot drive 1.ix cross-case on Linux CI without creating a real second file.

But the **property 1.ix depends on** — that the key folds case — is already pinned, platform-independently, by `io::tests::identity_key_is_case_sensitive_but_collision_key_is_not` (`src/io.rs:1181-1194`). So there is no coverage gap to apologise for, and no new key-level test to write. The plan should stop describing this as a limitation and instead cite the io test; that is both more accurate and the truthful replacement for the `norm_path` sentence (§1.1).

### 4.3 Gap — nothing asserts the *new* text is present

Validation 3 greps for the *old* strings and expects zero. Nothing confirms the new prefix landed in all three places. Cheap addition: grep the new prefix and expect exactly three hits (one message, two asserts). This catches a message/assert prefix mismatch at grep time rather than at `cargo test` time, and catches the case where someone updates the asserts to a prefix the message does not actually contain in the same wrapped-string form. Note the message is a wrapped Rust string literal with a `\` continuation, so the chosen pinned prefix must not straddle the line break — `"--passthrough matches R1 or R2"` is safe against the drafted wrapping, but any re-wrap during implementation must preserve that.

### 4.4 Gap — Validation 2 is an eyeball, not a check

"Read the rendered message for byte-equal and case-variant inputs" is a human review step with no artefact. It will pass by construction. If §1.4 is adopted, it becomes a real test at no cost: `test_passthrough_rejects_pointing_at_r1` asserts the message names R1, `…_at_r2` asserts R2. That converts the one genuinely under-tested aspect of the new message — that it discriminates which input matched — into a regression pin, and it is the only new test worth writing in this change.

### 4.5 No silent-wrong-result surface exists

Worth stating plainly for the risk register: this change cannot make the program produce wrong output. It cannot alter which runs are refused, cannot alter trimming, and its worst realistic failure is an unhelpful or untrue sentence on stderr. That is why the *wording review* is the substance of the validation here, and it is also why §1.1 and §1.5 are graded as high as they are — the message text is the entire deliverable, so a false or misdirecting sentence is not cosmetic, it is the bug.

---

## 5. Alternatives

All within the maintainer's decision (fold retained, no key change).

### 5.1 Minimal — strip the claim, add nothing

`"--passthrough matches R1 or R2 under case-insensitive comparison: {}"`. Satisfies #389's letter with the smallest possible diff and no new prose to get wrong. Loses all remediation, so a user who hits it byte-equal gets a correct statement and no guidance. Acceptable fallback if the preference is the tiniest diff; weaker than 5.3 for the same review cost.

### 5.2 As drafted

Truthful for both subcases, one remediation sentence aimed at the rarer one. Trade-off analysed in §1.5 — the drafted form's cost is that its only advice is wrong-footed for the likely user.

### 5.3 RECOMMENDED — two-tier message, branch on which subcase matched

The code can distinguish the subcases with a byte comparison it already has the operands for: `pt == &self.input[i]` (or `path_identity_key` equality, if `./` spellings should count as byte-equal) versus keys-equal-but-paths-not. Same key, same rejection, same order — the branch chooses only the sentence:

- **byte-equal:** name the matched read and say the passthrough must be a third file, e.g. the index read.
- **case-only:** state the comparison is case-folded for APFS/NTFS safety, name both spellings, and give the rename remediation.

This resolves §1.4 and §1.5 together, gives each subcase advice that is correct *for that subcase*, and remains entirely inside "wording only" — no behavioural change whatsoever, since both arms `bail!`. It costs one `if` on the error path and makes Validation 2 testable per §4.4. The extra cost over 5.2 is a few lines of diff in a change that is already being reviewed by two reviewers; the benefit is that the message stops giving the common case advice that routes around the guard.

### 5.4 Optional adjunct — delete the redundant `len == 2` guard

`src/cli.rs:938`'s `if self.input.len() == 2` is provably always true given 1.ii (`src/cli.rs:896-902`). Removing it and de-indenting the body is behaviour-preserving and makes §1.3's note unnecessary. Against: it widens a strings-only diff into a control-flow diff, and defensive redundancy next to a data-loss guard is a defensible thing to keep. My recommendation is to **keep the guard and keep the one-clause comment** (§1.2's two-liner does both) — but if the maintainer prefers the guard gone, this is the moment, because the comment justifying it is being rewritten anyway.

### 5.5 Considered and rejected — probe filesystem case sensitivity at runtime

Would let the message be platform-accurate ("these are two files here, so this is refused conservatively" vs "these are one file"). Rejected: it requires filesystem access, and #385 already rejected syscall/symlink semantics for these keys by design (A2). It would also make the message non-deterministic across platforms, breaking the very unit tests this plan is updating. Recorded for completeness; do not pursue.

---

## 6. Action items

### Critical

1. **Fix the `norm_path()` premise before implementing.** `norm_path` is live at `src/io.rs:48`, is a component of `collision_key`, and the existing test comment describing it as "a single `to_ascii_lowercase()` call" is *true*. Correct the plan's Context bullet and Implementation step 3, and replace the sentence by re-anchoring it on `io::tests::identity_key_is_case_sensitive_but_collision_key_is_not` (`src/io.rs:1181`) rather than by substituting the function name — the literal substitution ships a false comment. (§1.1)

### Important

2. **Cut the new comment to two lines**, fact-stated, per the repo's comment rule; move the dual-consume derivation to the commit message. Drop the ambiguous "elsewhere". (§1.2)
3. **Retain the `len == 2 per 1.ii` clause** so the redundant guard at `src/cli.rs:938` stays explained — or take §5.4 and remove the guard. (§1.3, §5.4)
4. **Name which input matched** in the message; the code already knows, and the sibling output-collision message sets the precedent (PR #219, `src/io.rs:120-126`). (§1.4)
5. **Fix the remediation for the byte-equal subcase** — "pass a third file (e.g. the index read)", not "rename one". As drafted, the advice lets the likely user walk around the guard instead of fixing the mistake. Adopting §5.3 resolves this and item 4 together. (§1.5, §5.3)
6. **Extend scope to the two other mis-framing comments**: `src/cli.rs:885-889` (calls 1.ix "issue #216 protection"; #216 is the output-collision issue, 1.ix guards input dual-consume) and `src/io.rs:38-47` (`norm_path` doc claims direct use by `Cli::validate`, and uses the "aliasing" identity verb). One line each; without them #389 is two-thirds done. (§1.6)
7. **Restate "no behaviour change" precisely** as "no change to which runs are accepted or rejected; the message text changes", in both the plan and the CHANGELOG entry. (§2.3)

### Optional

8. **Add a positive grep to Validation** — new prefix present exactly three times (message + two asserts); ensure the pinned prefix does not straddle the string literal's `\` line break. (§4.3)
9. **Turn Validation 2 into assertions** — have `…_at_r1` / `…_at_r2` pin that the message names the correct read. Only worthwhile if item 4 is adopted, and then it is the one test worth adding. (§4.4)
10. **Match the house qualifier idiom** — `(case-insensitive, for APFS/NTFS safety)`, as used at `src/io.rs:110` and `src/io.rs:121` and grepped by `.github/workflows/ci.yml:638`, rather than the drafted "compared case-insensitively, for APFS/NTFS safety". One vocabulary for one concept. (§1.4)
11. **Avoid `: {}.`** — a path immediately followed by a period is ambiguous to copy-paste. House style interpolates paths mid-sentence followed by a comma or a word (`src/io.rs:112`). (§5.2)
12. **Consider making `norm_path` private.** It is `pub` with no callers outside `collision_key` and its own tests; demoting it (or inlining it) would remove the naming confusion that produced item 1 in the first place. Out of this plan's scope — worth a follow-up issue rather than scope creep here. (§1.6)
13. **Draft the CHANGELOG entry in house style** — bolded lead sentence plus issue link, matching the neighbouring entries in Unreleased → `#### Changes` (`CHANGELOG.md:129`). The `--output-dir` entry at `CHANGELOG.md:145-146` is the closest precedent for a message-only change. (§1.7)

### Non-issues (verified, do not re-flag)

- A1 holds: only `src/cli.rs:944-945` and `src/cli.rs:1880`/`1887` carry the text; `tests/integration_passthrough.rs` pins 1.i only; no CI grep touches this message. `ci.yml:638`'s grep targets the output-collision message and cannot be tripped.
- Check ordering, the 1.viii-before-1.ix consequence, and both platform subcases are as the plan describes.
- Efficiency is genuinely nil.
- "Cross-case `validate()` testing needs a case-insensitive filesystem" is correct — but the fold property is already pinned at the key level, so there is no coverage gap (§4.2).
