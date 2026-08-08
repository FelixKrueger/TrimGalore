# Plan Review B — demote `io::norm_path` to private (#399)

**Plan:** `plans/08082026_normpath-visibility/PLAN.md`
**Base verified:** `dev` @ `ce13761` (matches the plan's stated base; `git status --porcelain src/ CHANGELOG.md Cargo.toml` clean).
**Reviewer:** B (independent; no shared state with Reviewer A).

**Verdict up front:** the change is safe and correct, and I verified it compiles rather than
taking the plan's word for it. Two of the plan's supporting *reasons* are wrong on the facts
(the validation net is not `cargo build`, and a dead-code warning is impossible here), and one
substantive omission is worth fixing inside this change: the case-folding trade-off note is
about `collision_key`'s behaviour, not `norm_path`'s, and privatization deletes it from the
rendered public docs. Details below.

---

## 1. Logic review

### 1.1 Is `norm_path` genuinely callerless outside `src/io.rs`? — Yes, verified four ways

The plan's grep claim holds. I checked the places a naive `io::norm_path` grep misses:

| Surface | Result |
|---|---|
| `git grep norm_path -- src tests benches examples build.rs Cargo.toml docs '*.md'` | only `src/io.rs` (def `:45`, sole call `:82`, test `:765-780`) + one `//` comment at `src/cli.rs:1903` |
| Untracked files too (`grep -rn`, excluding `plans/`, `target/`, `optimus_prime/`, `target_repro_b/`) | same set — nothing new |
| `src/lib.rs` | 17 bare `pub mod` lines, **no `pub use` re-exports**, no prelude. So `io::norm_path` *is* public surface today (`trim_galore::io::norm_path` resolves), but nothing consumes it. |
| `benches/`, `examples/` | do not exist |

Two details that matter and that a `io::norm_path`-shaped grep would have missed:

- `src/specialty.rs:12` and `src/clump_only.rs:52` both do `use crate::io as naming;`, so a caller
  there would read `naming::norm_path`. Historical plan text does describe `naming::norm_path`
  call sites, so this was a live possibility, not a hypothetical. Neither file contains the
  identifier — grepping the bare symbol (as I did) covers the alias; grepping the qualified path
  would not have.
- `src/cli.rs:1903` is a `//` comment (not `///`) naming the **test**
  `io::tests::test_norm_path_case_folds`, not the function. The test keeps its name, so this is
  unaffected. The plan says this too and is right.

### 1.2 Does it compile? — Yes, empirically

I did not rely on reasoning. I extracted `git archive HEAD` into a scratch copy (the shared tree
was not touched), applied `pub fn norm_path` → `fn norm_path` there, and ran
`cargo check --all-targets --locked`:

```
Finished `dev` profile [unoptimized + debuginfo] target(s)   # exit 0, 0 warnings
```

All 11 integration tests were genuinely compiled, not skipped — `.rmeta` artifacts exist for
`integration_adapter2` … `integration_ubam` in the scratch target dir. (My first log looked
suspiciously short at 40 lines; that was my own `| tail -40`, not a truncated build.)

**A2 is proven, not assumed:** `src/io.rs:565-567` is `#[cfg(test)] mod tests { use super::*; … }`,
and the green `--all-targets` run compiles that module, so the private fn is reachable from the
tests exactly as the plan says.

### 1.3 Rustdoc intra-doc links — not a hazard here, and not gated by CI either

Two independent reasons this cannot bite:

1. **No intra-doc links to `norm_path` exist anywhere.** The only cross-file mention is the `//`
   comment above. Plain backticks are not links — rustdoc resolves only bracketed forms
   (`[`norm_path`]`), and there are none. The `private_intra_doc_links` lint direction is also
   inverted from the worry: `norm_path`'s own doc references `collision_key`, i.e. a private item
   pointing at a public one, which is always fine.
2. **There is no rustdoc job in CI.** `ci.yml:160` "Docs build (Astro)" is the documentation
   *website* (`working-directory: docs`, `npm ci && npm run build`), and `docs.yml` deploys that
   same Astro site. No `cargo doc`, no `RUSTDOCFLAGS`, no `--document-private-items` anywhere in
   `.github/workflows/`. So no rustdoc warning can gate this PR.

Nothing beyond a compile error can break. The one *non-breaking* consequence is real and is my
Important finding — see §1.5.

### 1.4 `fn` (private) vs inlining into `collision_key` — `fn` is right, but for a different reason than the plan gives

I agree with the resolution. I do not agree with the argument.

**The plan's argument is the weaker one available.** It says inlining "dissolves the named seam
`test_norm_path_case_folds` pins". But the fold is *already* pinned at the public boundary, twice:

- `src/io.rs:1177-1191` `identity_key_is_case_sensitive_but_collision_key_is_not()` — asserts
  `collision_key` folds while `path_identity_key` does not, and carries a REGRESSION note about
  the #216 CI guard.
- `src/io.rs:1194-1206` `both_keys_normalise_spelling()` — pins spelling normalisation on both keys.

So losing `test_norm_path_case_folds` would not leave the fold unpinned, and "the seam the test
pins" overstates what that test uniquely protects.

**The argument that does survive scrutiny** (worth putting in the plan or the commit message
instead): `collision_key` **absolutises before folding**, so re-expressing the fold test in terms
of `collision_key` makes its relative-path assertions cwd-dependent. Today's test asserts
`norm_path(Path::new("foo.fq.gz")) == "foo.fq.gz"`; through `collision_key` that becomes
`"<lowercased cwd>/foo.fq.gz"`. Two of the test's four assertions use relative paths and would
not translate cleanly. A named pure-fold function is what lets the fold be asserted as a property
without dragging the working directory into the test. That is a real reason to keep the seam, and
it is cheaper than inlining regardless.

### 1.5 The trade-off note is mis-homed, and privatization makes that consequential — **Important**

Current doc layout in `src/io.rs`:

- `norm_path` (`:38-44`) carries the substantive paragraph: *"Pragmatic trade-off: on opt-in
  case-sensitive APFS volumes this may false-positive, but the penalty is a loud early error
  rather than silent data loss."*
- `path_identity_key` (`:72-77`) explains its case-**sensitive** rationale well.
- `collision_key` (`:79-83`) is two lines and carries **none** of the trade-off.

The note is mis-homed *already*, before this change. `norm_path` is a pure
`to_string_lossy().to_ascii_lowercase()`; "may false-positive" is not a property of lowercasing a
string. It is a property of *deciding output collisions on a folded key* — which happens in
`collision_key` and `preflight_output_collisions`, not in `norm_path`.

Privatization turns that latent mis-homing into a concrete loss:

- **rustdoc omits private items by default**, so after this change the only statement of the
  trade-off vanishes from the generated docs entirely. Every remaining *public* path-key doc is
  silent on it.
- An in-source reader who arrives at `collision_key:81` must scroll *up* past `path_identity_key`
  and `lexical_normalise` to find the reasoning behind the key they are reading.

**Recommendation:** move the trade-off paragraph onto `collision_key` as part of this change, and
reduce `norm_path`'s doc to one line. Two hunks instead of one, still zero behaviour change. This
converts the change from "hide the confusing thing" into "put each fact where its reader is",
which is a materially better answer to #399's own framing. Note this also retires the plan's
stated reason for keeping the named fn ("the doc comment's trade-off note lives well on the named
fn") — replace it with the cwd argument in §1.4.

### 1.6 Does #399 resolve the re-audit's "`norm_path`'s doc omits `path_identity_key`" finding, or hide it?

Split the finding in two; the answer differs per half.

- **The misuse half — genuinely resolved.** The reason to cross-reference `path_identity_key` from
  `norm_path` was to stop a caller reaching for the folded helper when they wanted the
  case-preserving key. That is exactly the mistake #389 had to correct, and after privatization no
  code outside `io.rs` *can* make it. The guarantee moves from "a doc a reader might skim" to
  "the compiler rejects it" (`error[E0603]`, demonstrated in §3.1). That is strictly stronger than
  the cross-reference would have been.
- **The informational half — only resolved if the note moves.** A reader of the public API still
  needs to learn that two keys exist with deliberately different case semantics. Most of that is
  already served: `collision_key` says "Case-**folded**", `path_identity_key` says
  "Case-**sensitive**, because …", the two are adjacent, and a test names the contrast. What is
  *not* served after privatization is the trade-off sentence.

**So:** treating #399 as the fix for the re-audit's finding is right about the misuse surface and
right about the doc gap **only if §1.5 is adopted**. If the trade-off note stays on a private fn,
#399 hides that finding rather than resolving it. That is the crux of this review.

### 1.7 The "(private)" doc framing — mild objection (Optional)

Plan step 1 proposes the first doc line "gains '(private)' framing". I would not. The `fn`
keyword three lines below already states the visibility, and docs that narrate their own
visibility go stale the next time visibility changes. `CLAUDE.md`'s comment convention ("state the
fact, not the evidence", default to one line) argues for *shortening* this doc block, not adding a
parenthetical to it. Combined with §1.5 the natural result is a single line:
`/// Case-folded (ASCII lowercase) view of a path; the folding half of `collision_key`.`

---

## 2. Assumptions

| Plan assumption | Status |
|---|---|
| **A1** — zero callers outside `io.rs`, grep-verified on `src/`, `tests/` | **Confirmed**, and widened: also `benches/`/`examples/` (absent), `build.rs`, `Cargo.toml`, `docs/`, all untracked non-`plans/` files, and the `use crate::io as naming` aliases in `specialty.rs`/`clump_only.rs`. |
| **A2** — `mod tests` reaches private items | **Confirmed** by compile, not just by language rule (`io.rs:565` is `#[cfg(test)] mod tests` with `use super::*`). |

Unstated assumptions the plan relies on, all of which I checked and all of which hold:

- **No rustdoc gate in CI.** True (§1.3). Had `cargo doc -D warnings` existed, the plan's
  validation table would have had a gap.
- **No doc-tests reference the symbol.** True, and structurally impossible to matter here: the
  crate's only doc code fence is a ```` ```text ```` block at `src/specialty.rs:235-240`, so there
  are no compiled doc examples at all. Worth knowing *why* this matters: `--all-targets` excludes
  doc-tests, so a doc example using `norm_path` would slip past clippy and surface only under
  `cargo test`. Non-issue on this base.
- **The coverage CI job has no threshold.** True — `ci.yml:271-310` runs `cargo llvm-cov` and
  uploads LCOV with no fail-under gate, so a line-attribution shift cannot fail CI.
- **The crate has no `#![deny(...)]`/`#![warn(missing_docs)]` lint attributes** that could react to
  a visibility change. Confirmed: `src/lib.rs` and `src/main.rs` carry no inner attributes at all.

The plan's line-number citations are all accurate on this base (`:45-47` def, `:82` caller,
`:765-780` test).

---

## 3. Validation sufficiency

The proposed command set is sufficient. The plan's *explanation* of why it is sufficient is wrong
in two places, and one of them could cause a future reader to weaken the check.

### 3.1 `cargo build` is NOT the safety net — clippy `--all-targets` is (**Important**)

Validation row 1 reads "`cargo build` + clippy `-D warnings` → clean (an external caller would be
a compile error)". The first half is false. `cargo build`/`cargo check` compile lib + bin only,
**not** `tests/`. I ran the positive control rather than reasoning about it — added a new
`tests/zz_probe.rs` calling `trim_galore::io::norm_path` against the private variant:

```
A) cargo check --locked                → EXIT 0   ("Finished dev profile")   ← MISSES it
B) cargo check --all-targets --locked  → EXIT 101
   error[E0603]: function `norm_path` is private
     --> tests/zz_probe.rs:6:33
```

Removing the probe returned `--all-targets` to green (exit 0), so the harness demonstrably can
both pass and fail — the §1.2 green result is trustworthy.

Practically the plan is fine, because it *does* run `cargo clippy --all-targets` and `cargo test`,
either of which catches this. The fix is to the wording, so nobody later trims the run to
`cargo build` believing that is the guard. Since `tests/` is precisely where a hidden consumer of
a `pub` library item would live, this is the one row where the stated reasoning inverts the risk.

### 3.2 The dead-code claim is backwards (**Important**, wording only)

Plan step 2: "a dead-code warning would prove a missed caller assumption — none expected,
`collision_key` uses it". `dead_code` fires when **nothing** uses an item — the opposite of a
missed caller — and it cannot fire here at all, because `collision_key` calls it. The clean run
confirms: **0 warnings**. The lint that speaks in the failure mode is `E0603`, a hard error, and
only on targets that actually get compiled. Recommend deleting the parenthetical rather than
rewording it.

### 3.3 Rows 2 and 3

Both fine. `cargo test test_norm_path_case_folds` is unaffected by visibility (child module), and
full `cargo test` is the belt-and-braces catch for the doc-test blind spot in §3.1 that does not
apply here anyway. `cargo fmt --all -- --check` is listed and is worth keeping — the edit shortens
a line and could in principle disturb rustfmt's doc-comment wrapping if §1.5's doc rewrite is
adopted.

### 3.4 Nothing material is missing

I looked for a gap and did not find one. The whole risk surface of a `pub` → private demotion is
"does every compilation unit that could reference it still compile", and `clippy --all-targets` +
`cargo test` cover every unit this crate has. No runtime behaviour exists to test.

---

## 4. Efficiency

Nil, as the plan says, and correctly dismissed. `norm_path` stays a separate function called once
per path; privatization marginally widens the optimizer's freedom (a non-exported function is
easier to inline away) but the call is already within one crate and one translation unit, so the
effect is immaterial and not worth mentioning in the commit message.

---

## 5. CHANGELOG: skip is right, the stated reason is not

The plan's conclusion (no entry) is correct. Its justification — "consistent with the repo treating
only behaviour/message changes as entries" — is not accurate, and the real precedent is stronger:

- **A slot for non-user-visible changes exists.** `CHANGELOG.md:481` — the Unreleased section
  already has `#### Infrastructure (contributor-facing)`, with two entries (the new `docs-build`
  CI job; the `justfile` docs recipes that ran `cd Docs`). v2.1.0-beta.2 had
  `#### Infrastructure (contributor-facing, no runtime effect)` (CI gates, Dependabot config). So
  "only behaviour/message changes" is false as stated.
- **But this change is below that section's granularity floor.** Those entries are things a
  contributor trips over — a new PR gate, a recipe broken on Linux, Dependabot routing. A
  visibility keyword on an internal helper is not.
- **The decisive precedent:** `grep -n "collision_key\|norm_path\|path_identity_key" CHANGELOG.md`
  returns **zero hits** across 1300+ lines that document the entire #216 / #383 / #384 / #388 /
  #389 / #391 collision-key arc — all of which touched these very functions. This changelog has
  never named these helpers, so it certainly should not start for a `pub` removal.
- Supporting: `7ea2741` (docs-only, #387) merged with no CHANGELOG entry, so not every PR gets one.

**Recommendation:** keep the decision, fix the reason to something like *"below the Infrastructure
section's granularity floor; the CHANGELOG has never named these helpers, across the whole
collision-key arc."* An accurate precedent is worth more than a convenient one here, because the
inaccurate version could be cited later to skip an entry that *does* belong.

---

## 6. Alternatives

| Alternative | Assessment |
|---|---|
| **`fn` (plan's choice)** | Correct. Sole caller is in the same module, so module-private is the tightest visibility that works. |
| **`pub(crate) fn`** | Strictly worse *here*. It would keep the symbol reachable from `specialty.rs`/`clump_only.rs` (which already alias `crate::io as naming`) — i.e. it preserves most of the confusion surface #399 exists to remove, while buying nothing, since no cross-module caller exists or is wanted. Reject explicitly if it comes up in review. |
| **Inline into `collision_key`** | The issue's alternative. Viable but worse — see §1.4; the cwd-dependence argument, not the plan's seam argument, is what kills it. Also loses a cheap property test for no gain. |
| **`#[inline]` + private** | Pointless noise; LLVM handles a one-liner in-crate. |
| **Private + move the trade-off doc to `collision_key`** | **My recommendation** (§1.5). Same behaviour, same risk, one extra hunk, and it is the version that actually resolves the re-audit finding rather than hiding it. |
| **Do nothing / close #399** | Not recommended, but worth naming the honest cost: the change buys clarity, not correctness. Given #389's planning error traced directly to the misleading `pub`, and given the compiler now enforces the misuse guard (§1.6), the clarity is worth one keyword. |

---

## 7. Action items

### Critical
None. The change compiles, behaviour is unchanged, and the plan's proposed commands do catch the
only failure mode.

### Important

1. **Move the case-folding trade-off paragraph from `norm_path` to `collision_key`** (`src/io.rs`
   `:42-44` → `:79-83`), and shorten `norm_path`'s doc to one line. Without this, privatization
   deletes the trade-off from the rendered public docs and #399 hides the re-audit's doc finding
   instead of resolving it (§1.5, §1.6). Small, zero-behaviour, and it is the difference between a
   good change and a merely harmless one.
2. **Fix validation row 1: the net is `clippy --all-targets` / `cargo test`, not `cargo build`.**
   Empirically `cargo check` (same target set as `cargo build`) passes with an external caller in
   `tests/` present; only `--all-targets` raises `E0603` (§3.1). Keep both commands — just stop
   crediting the wrong one, so the step is not trimmed later.
3. **Delete the dead-code parenthetical in step 2.** `dead_code` signals *no* caller, cannot fire
   while `collision_key` calls the fn, and the clean run emits 0 warnings (§3.2).
4. **Correct the CHANGELOG rationale** (keep the decision). The repo *does* have an
   `Infrastructure (contributor-facing)` slot in the current Unreleased section; the real reason to
   skip is granularity plus the zero-hit history for these three symbol names (§5).

### Optional

5. **Replace the plan's "seam" justification with the cwd-dependence argument** (§1.4). The current
   one is refuted by `identity_key_is_case_sensitive_but_collision_key_is_not` and
   `both_keys_normalise_spelling`, which already pin the fold publicly; the surviving argument is
   that `collision_key` absolutises, so a fold test written through it becomes cwd-dependent.
6. **Drop the "(private)" doc framing** (§1.7). The `fn` keyword states it; self-narrating
   visibility goes stale and cuts against `CLAUDE.md`'s one-line comment convention.
7. **Say `pub(crate)` was considered and rejected** in the commit message — it is the obvious
   reviewer question and the answer (it keeps the `naming::` alias surface open for nothing) is
   worth one clause.
8. **Optionally note in the plan that no rustdoc gate exists** (§1.3). It is the assumption that
   would have made this change riskier than it looks, and recording that it was checked saves the
   next person the same audit.

---

## Appendix — how the claims here were verified

- `git grep` + `grep -rn` over `src/ tests/ build.rs Cargo.toml docs/` and all untracked
  non-`plans/` files; separate check of the `use crate::io as naming` aliases.
- `src/lib.rs` read in full (17 `pub mod`, no re-exports).
- `.github/workflows/{ci,docs,release}.yml` grepped for `cargo doc` / `rustdoc` / `RUSTDOCFLAGS` /
  `document-private` — no hits; `ci.yml:160` job confirmed to be Astro.
- **Isolated compile** (`git archive HEAD` → scratchpad, `pub fn` → `fn`, fresh
  `CARGO_TARGET_DIR`): `cargo check --all-targets --locked` exit 0, 0 warnings; `.rmeta` present
  for all 11 integration tests.
- **Positive control**: probe caller in `tests/` → `cargo check` exit 0, `cargo check
  --all-targets` exit 101 `E0603`; probe removed → exit 0 again.
- `CHANGELOG.md` heading map + `grep` for the three helper names (zero hits); last 25 commits
  classified by whether they touched `CHANGELOG.md`.
- The shared working tree was **not** modified — every experiment ran in
  `…/scratchpad/npcheck1`.
