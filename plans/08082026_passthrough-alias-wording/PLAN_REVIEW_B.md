# Plan Review B — `--passthrough` alias check wording (#389)

**Reviewer:** B (independent, fresh context)
**Plan under review:** `plans/08082026_passthrough-alias-wording/PLAN.md`
**Repo:** `/Users/fkrueger/Github/TrimGalore`, branch `dev` @ `03d793f` (verified as HEAD)
**Date:** 2026-08-08
**Scope note:** The maintainer's decision — keep the case-folded `collision_key`, change wording only — is treated as settled. This review evaluates the plan *given* that decision.

**Verdict (short):** the decision is sound and I verified its premise independently. The plan's factual groundwork is almost entirely accurate, but it has three defects that matter: its automated validation does not actually protect the thing #389 asked for; it misses the one *other* place in `src/` that still carries both the objected-to framing and a now-false attribution; and step 3, followed literally, replaces a stale identifier with a **false** statement.

---

## 1. Logic review

### 1.1 Claims I verified as correct

Every line reference in the plan checks out against the working tree:

| Plan claim | Actual | OK |
|---|---|---|
| check 1.ix at `src/cli.rs:933-949` | comment opens 933, block closes 949 | yes |
| comment to replace at `933-937` | 5 lines, "…catches case-only INPUT aliases…" | yes |
| message at `943-947` | `bail!` spanning 943-947 | yes |
| tests at `cli.rs:1871-1888` | `test_passthrough_rejects_pointing_at_r1` 1871-1881, `_r2` 1883-1888 | yes |
| shared test comment at `1872-1877` | 6 lines, names `norm_path()` | yes |
| base `03d793f` | `git log -1` → `03d793f` | yes |
| 1.ix runs after 1.viii existence check | `check_restartable_input(pt, …)` at 932 | yes |
| `tests/integration_passthrough.rs` does not pin this text | its only stderr pin is `"--passthrough requires --paired"` (1.i) | yes |
| CI does not exercise `--passthrough` message text | `grep -n passthrough .github/workflows/*.yml` → **zero hits** | yes |
| no docs/README pin | no `alias`/`cannot point at` hit for passthrough in `docs/` or `README.md` | yes |

The `case-insensitive match` grep hit at `.github/workflows/ci.yml:635` is a **false alarm for this plan** — it is a historical comment about the *output-collision* message (PR #219 reworded it), and the live CI grep is on `"case-insensitive, for APFS/NTFS safety"` against an output-collision log, not a passthrough run. No CI pin exists on 1.ix. The plan's Integration section is right.

### 1.2 The decision's premise holds — verified independently

The plan asserts that dropping the fold would reintroduce silent dual-consume. I checked whether anything *else* would catch it, because if the output pre-flight already did, the rationale would collapse:

- passthrough output is `<stem>_passthrough.fq(.gz)` (`src/io.rs:321-340`), and it *is* pushed into the pre-flight candidate list (`src/main.rs:763-764`).
- but `io::preflight_output_collisions` (`src/io.rs:96-130`) only compares **output-vs-input** and **output-vs-output**. With `--passthrough R1.fq.gz` (or `r1.fq.gz` on APFS), the planned outputs are `R1_val_1`, `R1_val_2`, `R1_passthrough` — none of which alias an input or each other.

So nothing downstream catches input-vs-input aliasing. **1.ix is the sole guard**, and keeping the fold is the correct call. Good — the plan's rationale is not just asserted, it survives checking.

### 1.3 C3 — step 3, followed literally, writes a *false* comment

This is the sharpest defect. Plan line 15 says the test comment "references **`norm_path()`**, a stale pre-D13 name for what is now `collision_key`."

**That premise is wrong.** `norm_path` is not a stale name — it still exists, still `pub`, at `src/io.rs:48-50`, and is still exactly one `to_ascii_lowercase()`:

```rust
pub fn norm_path(p: &Path) -> String {
    p.to_string_lossy().to_ascii_lowercase()
}
```

`collision_key` *calls* it: `norm_path(&lexical_normalise(p))` (`io.rs:85`).

The real defect in the test comment is different, and the difference is load-bearing. The comment currently argues:

> `norm_path()` is a single `to_ascii_lowercase()` call which is trivially correct

That sentence is **true of `norm_path`** but the check no longer calls `norm_path` — it calls `collision_key`. Step 3 instructs "replace the stale `norm_path()` reference with `collision_key`". Do that mechanically and the committed comment reads "`collision_key()` is a single `to_ascii_lowercase()` call which is trivially correct" — which is **false**: `collision_key` also absolutises via `std::path::absolute` and folds `..` against the component stack (`lexical_normalise`, `io.rs:55-73`). Step 3 would replace a stale-but-true statement with a current-but-false one, in a change whose entire purpose is removing a false statement from this exact block.

Fix: rewrite the clause, don't swap the identifier. The byte-equal-is-a-subset argument still holds and can be stated without claiming what the function does — e.g. "byte-equal paths yield equal keys under any normalisation, so this covers the folded path too; the fold itself is pinned by `io::tests::test_norm_path_case_folds`." That last pointer is worth having, because the fold genuinely *is* tested there (`io.rs:767-784`) — which also answers the comment's "CI can't test cross-case" worry properly.

### 1.4 C2 — the goal is not met: `src/io.rs:38-47` survives untouched

The plan's Goal is "make the check's words match its semantics", and its A1 grep concluded the old framing lives in exactly four lines of `cli.rs`. It does not. `norm_path`'s own doc comment says:

```rust
/// Case-folded (ASCII lowercase) string view of a path for collision detection
/// on case-insensitive filesystems (APFS/NTFS). Used by:
///   * `Cli::validate()` to catch `--passthrough` aliasing R1 or R2 (e.g.
///     `--passthrough r1.fq.gz` while R1 is `R1.fq.gz`).
```

Two problems in one bullet:

1. It carries **the same "aliasing" framing #389 objects to** — and here it is on a `pub fn`'s rustdoc, i.e. more durable than an inline comment.
2. It is **factually stale as an attribution**: `Cli::validate` does not call `norm_path`. `grep -rn norm_path --include='*.rs'` shows exactly three non-test call sites — its definition (48), `collision_key` (85), and the test (768-783). Post-#385, `validate` reaches the fold only *through* `collision_key`, which the second bullet already describes.

Leaving this in means the fix is cosmetic: the objected-to sentence still ships, one file over. A secondary, milder instance sits at `io.rs:769-771` ("the output-collision pre-flight + `--passthrough` validation both use") — that one is transitively true and can stay, though it reads oddly once bullet 1 above is corrected.

Root cause: **A1's grep was phrase-narrow.** It searched two exact strings from `cli.rs`; the same idea spelled "aliasing R1 or R2" was invisible to it. See I5.

### 1.5 The "no behaviour change" claim

Holds in the sense that matters: no key, no comparison, no ordering, no control flow changes; the rejection surface is bit-identical. Two honest caveats worth putting in the plan rather than leaving implicit:

- **stderr text is user-visible behaviour.** Anyone's wrapper script grepping `"cannot point at R1 or R2"` breaks. That is the intended effect, not a bug, but it is why the CHANGELOG entry is not optional. The plan already has the CHANGELOG step, so this is a wording nit on the claim, not a gap.
- If I2's "name both paths" refinement is adopted, the diff stops being string-only (see §5).

### 1.6 Edge cases

The plan says "none new" and that is right. One pre-existing ordering artifact becomes slightly more visible and is worth a sentence: R1/R2 existence is checked at `cli.rs:987-989`, **after** the passthrough block, so on a case-sensitive filesystem with a genuinely-missing `R1.fq` and `--passthrough r1.fq`, the user sees the new "rename one so the paths differ" advice rather than "Input file not found". Pre-existing, low harm, not worth reordering — but "no new edge cases" is more precisely "no new edge cases; one pre-existing ordering artifact that the new remediation sentence makes marginally more confusing".

---

## 2. Assumptions

### Stated

- **A1 — "No test outside `cli.rs:1871-1888` pins the old text."** True for *test* pins (I re-ran the grep over `src/`, `tests/`, `docs/`, `.github/`, `README.md`, `CHANGELOG.md`: only `cli.rs:944-945`, `1880`, `1887`). But A1 was used to conclude the *comment* scope too, and there it fails — see C2. Narrow the claim to "no test pins" and add a separate scope grep for the framing.
- **A2 — "D13 taxonomy stays as-is."** Valid and consistent with the recorded decision. Nothing in the plan touches either key.

### Implicit assumptions the plan does not surface

- **The message is only ever read by humans.** Nothing machine-parses TrimGalore stderr in-repo (verified), but the assumption should be stated, because it is what licenses a free-form rewrite.
- **`plans/` is out of scope for the "zero stale hits" grep.** The plan does not say this, and it must — see I4.
- **The rejection is desirable in the common case.** The plan reasons hard about the rare case-sensitive-FS shape and barely about the overwhelmingly likely one (user fat-fingered R1 as the passthrough file). That asymmetry leaks into the proposed message — see I3.
- **A 4-line comment is acceptable here.** The project's own convention says otherwise — see I1.

### Ambiguity

- "CHANGELOG.md — Unreleased/changed" does not name an existing heading. `Unreleased` currently has **four** subsections: `#### Bug fixes` (line 6), `#### Changes` (129), `#### Fixes` (272), `#### Infrastructure (contributor-facing)` (458). With both "Bug fixes" and "Fixes" available, an implementer can easily file a no-behaviour-change edit as a bug fix. Name `#### Changes` explicitly.

---

## 3. Efficiency analysis

Nothing to analyse — no allocation, loop, or I/O change; `collision_key` is called the same two times on the same paths. The plan's one-line dismissal is proportionate.

Two non-performance mechanical hazards in the same territory, since this is where "nothing to measure" tempts an unreviewed diff:

- **Backslash-continuation whitespace.** The existing message uses `\` line continuations, which eat the next line's leading whitespace, so the trailing space must sit *before* the backslash. A three-line prose message is three chances to produce `case-insensitivelyfor` — and with the pin as currently specified (prefix only), **no test would catch it**. This is a second, independent reason to adopt C1.
- **`cargo fmt` and long string literals.** `fmt` will not rewrap a string literal, so the author owns the line breaks; keep them under 100 columns to match the surrounding block or clippy's line-length-adjacent lints in `-D warnings` mode become the discovery mechanism.

---

## 4. Validation sufficiency

The plan's four validations catch a partial test-pin update, formatting, and clippy. They do **not** protect the deliverable. Two gaps:

### C1 (Critical) — the proposed pin does not protect the #389 fix

Validation 1 relies on the two unit tests, and step 3 pins them to `contains("--passthrough matches R1 or R2")`. That prefix says nothing about honesty. A future "tidy up the error messages" commit could ship:

```
--passthrough matches R1 or R2: sample_R1.fq.gz aliases an input file
```

…and both tests stay green. The exact defect #389 was filed about would be reintroduced with a green suite. For a change whose entire deliverable *is* the wording, the wording has to be the contract.

Pin the honest part, positively and negatively:

```rust
assert!(err.contains("--passthrough"), "got: {err}");
assert!(
    err.contains("case-insensitiv"),          // the comparison is stated, not assumed
    "1.ix must say the comparison is case-insensitive (#389); got: {err}"
);
assert!(
    !err.contains("aliases an input"),        // no unconditional identity claim
    "1.ix must not assert identity (#389); got: {err}"
);
```

The negative assertion is the one that actually holds the line, and the `#389` in the failure message tells the next person why the test exists. Add it to both tests, or factor the three asserts into a small helper the two tests share.

### Gap 2 (Important) — no validation covers the scope question

Validation 3 greps for the two *old* strings. Nothing greps for the *framing*, which is why C2 slipped through. Add:

```bash
grep -rn "aliasing R1 or R2\|aliases an input\|cannot point at R1" \
  src/ tests/ docs/ .github/ README.md CHANGELOG.md
```

### I4 (Important) — Validation 3's expected result is wrong as written

"grep for `"cannot point at R1"` and `"aliases an input"` post-change → zero hits" is unachievable and, worse, actively misleading. Unscoped, those two phrases hit **this plan itself** (4 lines), `plans/06062026_passthrough-mode/PLAN.md:114`, and about eight more files under `plans/08062026_se-output-collision-preflight/` and `plans/08072026_paired-report-preflight/`. An implementer chasing "zero hits" would either report a false failure or start editing historical plan artifacts — which the global workflow rules put off-limits. Scope the grep to `src/ tests/ docs/ .github/ README.md CHANGELOG.md` and state the expectation as "zero hits outside `plans/`".

### What is adequately covered

- The fold itself is already pinned by `io::tests::test_norm_path_case_folds` (`io.rs:767-784`), including the `R1.fq.gz`/`r1.fq.gz` equality. The plan does not mention this, but it means the "cross-case can't be tested on CI" caveat in the test comment is narrower than it sounds — the *key* is tested, only the end-to-end filesystem behaviour is not. Worth citing in the rewritten comment (see C3).
- Validation 2 is a human read, which is the right instrument for prose. Keep it, but make it a checklist of the two subcases rather than a vibe check: (a) byte-equal — does the message still tell the user what to do? (b) case-variant on ext4 — is the rename advice reachable and correct?

---

## 5. Alternatives (within the settled decision)

### 5.1 Name both paths, not just the passthrough one — I2

The proposed message interpolates only `pt.display()`, then says:

> If **these** are genuinely two different files on a case-sensitive filesystem, rename one…

"these" and "one" have no antecedent — only one path was printed. The user must guess whether R1 or R2 matched, which is precisely the information needed to act. The repo already set the opposite precedent for the sibling message: PR #219 reworked the output-collision error specifically to name **both** colliding paths, and `ci.yml:634-637` records that as a deliberate improvement.

Minimal restructure that yields the matched path:

```rust
let pt_key = crate::io::collision_key(pt);
if let Some(matched) = self
    .input
    .iter()
    .find(|p| crate::io::collision_key(p) == pt_key)
{
    anyhow::bail!(
        "--passthrough must not be one of the R1/R2 inputs; {} and {} compare equal \
         case-insensitively (for APFS/NTFS safety). If they are genuinely two different \
         files on a case-sensitive filesystem, rename one so the paths differ by more \
         than letter case.",
        pt.display(),
        matched.display()
    );
}
```

Rejection surface is unchanged (`1.ii` guarantees `input.len() == 2`, so iterating equals testing `[0]`/`[1]`), and the redundant `if self.input.len() == 2` wrapper falls away with it.

**Trade-off, stated plainly:** this converts a string-only diff into a small logic diff, which raises the review bar on a change Felix has already flagged as single-reviewer-weight. If keeping the diff string-only is worth more than naming the match, then fix the dangling pronoun instead of the plumbing:

> …rename either the passthrough file or the matching input so the paths differ by more than letter case.

Both are acceptable; I'd take the first, because "which input did I collide with" is the user's actual next question.

### 5.2 Keep the imperative for the common case — I3

"`--passthrough` **cannot point at** R1 or R2" was doing real work: it told the 99% user (who mistakenly passed R1 as the passthrough file) what to change. The proposed "`--passthrough` **matches** R1 or R2" is a statement of fact followed by advice that only applies to the pathological case. The message gets more honest and less useful at the same time.

"must not be one of the R1/R2 inputs" (as in 5.1) keeps the normative content, keeps the common-case remediation implicit-but-clear, and still makes no identity claim — the identity question is deferred to the conditional second sentence, which is exactly the structure #389 asked for.

### 5.3 Shorten the comment to project convention — I1

The proposed replacement is four lines carrying a reasoning chain ("…rather than risking silent dual-consume elsewhere" — and "elsewhere" is vague). The repo's own convention is explicit: default to one line, two maximum, state the fact not the evidence, and put the derivation in the commit message. The full dual-consume argument belongs in the commit body and in #389, both of which will outlive the comment.

```rust
// 1.ix — deliberately case-folded, not an identity check (#389): on APFS/NTFS a
// case-variant passthrough is R1/R2 and would be read twice.
```

Two lines, states the fact, keeps the issue pointer. Everything else the plan wants to say is already written down in the plan and the issue.

### 5.4 Trim the duplicate rationale at `cli.rs:885-889` — O1

The `--passthrough` envelope preamble already says:

> case-folded collision check (1.ix) uses `crate::io::collision_key` to share the same APFS/NTFS-aware normalisation as the output-collision pre-flight in `main.rs` (issue #216 protection).

Not false, but "issue #216 protection" is an *output* concern attached to an *input* check — the same conflation #389 is about — and after this change the file carries two rationales for 1.ix, nine lines apart, framed differently. Drop the 1.ix clause from the preamble (the new inline comment covers it) or reduce it to a pointer.

### 5.5 Considered and rejected

- **A third key (`input_alias_key`) or platform-conditional check.** Out of bounds per A2 and already rejected in #385; adds a taxonomy entry to buy a pathological case.
- **Downgrade 1.ix to a warning on case-only matches.** Superficially attractive (case-sensitive FS proceeds, APFS user gets told) but it cannot distinguish the two without a syscall, and a warning on APFS means shipping the silent dual-consume with a note. Correctly not on the table.

---

## 6. Action items

### Critical

- **C1 — pin the honesty, not just the prefix.** Add a positive assertion on the case-insensitivity qualifier and a **negative** assertion `!err.contains("aliases an input")` to both `test_passthrough_rejects_pointing_at_r1`/`_r2`, with `#389` in the failure message. Without this the plan's deliverable is unprotected and a future message tidy-up can reintroduce the exact defect with a green suite. (§4)
- **C2 — fix `src/io.rs:38-47`.** The `norm_path` rustdoc's first "Used by" bullet both repeats the objected-to "aliasing R1 or R2" framing and is factually wrong (`Cli::validate` calls `collision_key`, not `norm_path`). Drop or re-attribute that bullet; without it the plan's stated Goal is not met. Add the file to the Implementation outline. (§1.4)
- **C3 — rewrite the test comment's clause; do not just swap the identifier.** `norm_path` is *not* a stale name (it lives at `io.rs:48`, still one `to_ascii_lowercase()`), so plan line 15's premise is wrong, and following step 3 literally produces the **false** statement "`collision_key()` is a single `to_ascii_lowercase()` call" — `collision_key` also absolutises and folds `..`. Restate the byte-equal-subset argument without describing the function's internals, and point at `io::tests::test_norm_path_case_folds` for the fold. (§1.3)

### Important

- **I1 — cut the new comment to ≤2 lines** per the project's comment convention; move the dual-consume derivation to the commit message and #389. Concrete rewrite in §5.3.
- **I2 — name both paths in the message** so "these"/"one" have an antecedent and the user learns which input matched; follows the PR #219 precedent for the sibling output-collision error. Restructure in §5.1, with a string-only fallback if the diff must stay string-only.
- **I3 — keep the normative clause** ("must not be one of the R1/R2 inputs") so the common-case user still learns what to change; the current draft offers remediation only for the rare case. (§5.2)
- **I4 — fix Validation 3's expectation.** Scope the grep to `src/ tests/ docs/ .github/ README.md CHANGELOG.md`; unscoped it hits this plan plus about nine files under `plans/`, so "zero hits" is unachievable and invites editing historical artifacts. (§4)
- **I5 — widen the scope grep** beyond two exact phrases to catch the framing (`aliasing R1 or R2`), and correct A1 to claim only "no *test* pins outside `cli.rs`". This is the root cause of C2. (§2)

### Optional

- **O1 — trim the duplicate 1.ix rationale** at `cli.rs:885-889`, or its "(issue #216 protection)" tag, so the file carries one framing rather than two. (§5.4)
- **O2 — name the CHANGELOG subsection explicitly** (`#### Changes`, line 129); `Unreleased` has both `Bug fixes` and `Fixes`, and a no-behaviour-change edit should not land under either.
- **O3 — watch the `\` continuation whitespace** in the multi-line format string; a missing space before a backslash silently joins two words and, with the prefix-only pin, no test catches it. Largely mitigated by C1.
- **O4 — consider a distinct qualifier phrasing.** "(compared case-insensitively, for APFS/NTFS safety)" is near-identical to the CI-pinned output-collision string `"case-insensitive, for APFS/NTFS safety"` (`ci.yml:637`). No conflict today (different invocations, different logs), but a future grep-tightening could cross-match.
- **O5 — soften the "no new edge cases" claim** to acknowledge the pre-existing ordering artifact: R1/R2 existence is validated after 1.ix (`cli.rs:987-989`), so on a case-sensitive filesystem with a missing R1 the rename advice can appear in place of "Input file not found". (§1.6)
- **O6 — record in the plan that the fold is already unit-tested** (`io.rs:767-784`), which narrows the test comment's "CI can't test cross-case" caveat to end-to-end filesystem behaviour only.

### On review weight

The plan asks whether single-reviewer weight suffices. On the evidence: the diff is small, but its two-string surface concealed three Critical items — an unprotected deliverable, a missed file that defeats the stated goal, and an instruction that writes a false comment if followed literally. Small diff, not small review. If I2's restructure is adopted, the diff also stops being string-only.
