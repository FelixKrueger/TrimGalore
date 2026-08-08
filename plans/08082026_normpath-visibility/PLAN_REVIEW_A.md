# Plan Review A — demote `io::norm_path` to private (#399)

**Reviewer:** A (independent; Reviewer B reviewed the same plan separately)
**Plan:** `plans/08082026_normpath-visibility/PLAN.md`
**Base verified:** `dev` @ `ce13761` (matches the plan's stated base)
**Verdict:** Approve with two doc-content changes folded in (Important). No Critical findings.

---

## What I verified (rather than trusted)

Every factual claim in the plan was checked. The change was also **compiled for real** on an
isolated copy of the crate in the scratchpad — the shared working tree was not touched.

| Claim | Method | Result |
|---|---|---|
| Zero callers outside `src/io.rs` | Repo-wide grep incl. all file types, `tests/`, `examples/`, `build.rs`, `Cargo.toml`, `docs/`, `*.md` | **Confirmed.** Only `src/io.rs:45` (def), `:82` (`collision_key`), `:765-780` (test). Plus one *comment* at `src/cli.rs:1903` citing the test by name — unaffected. |
| Grep could actually have found `.md`/docs hits | Positive control: same grep shape for `collision_key` | **Confirmed** — returns `SESSION_HANDOFF.md` and many `plans/*.md`, so markdown was genuinely in scope. The `norm_path` negative is trustworthy. |
| No re-export widens the surface | Read `src/lib.rs` | **Confirmed.** 17 bare `pub mod` lines, **no `pub use`**. So `norm_path` is reachable today only as `trim_galore::io::norm_path`, and nothing uses that path. |
| Making it private compiles | `cargo clippy --all-targets -- -D warnings` on a scratchpad copy with `pub fn` → `fn` | **Clean, exit 0** (honest exit code — first run was pipe-masked by `tail`, re-run without the pipe). Covers lib + bin + all 11 integration tests + the example. |
| The compile-error net is real | Positive control: added an external `trim_galore::io::norm_path` caller in `tests/` | **Fires:** `error[E0603]: function `norm_path` is private`, exit 101. Removed afterwards. |
| A2 — `mod tests` reaches a private item | `cargo test --lib test_norm_path_case_folds` on the private-fn copy | **Passes** (`test io::tests::test_norm_path_case_folds ... ok`; 409 other lib tests filtered out). |
| No rustdoc intra-doc-link hazard | Grep for bracketed `[norm_path]` forms; audited CI | **None exist.** The file uses plain backticks throughout, never intra-doc links. Also: **CI has no rustdoc job at all** — `docs-build` is the *Astro* docs site (`.github/workflows/ci.yml:160`), and no workflow runs `cargo doc` or sets `RUSTDOCFLAGS`. So the hazard is both absent and ungated. |

**Bottom line on §1 of the brief: yes, it compiles, and nothing beyond a compile error can break.**
The rustdoc private-intra-doc-link hazard is real in principle but does not apply — see Important-2
for the one way the implementer could *introduce* it.

---

## Logic review

The plan is internally consistent and the mechanism is sound. Three corrections:

**L1 — Validation 1's stated rationale is self-contradicting (wording, not substance).**
The plan writes: *"`cargo build` (a dead-code warning would prove a missed caller assumption —
none expected, `collision_key` uses it)"*. A `dead_code` warning **cannot** fire here, for exactly
the reason given in its own parenthesis: `collision_key` at `io.rs:82` is a non-test caller. The
check is right; the reason is wrong. The actual safety net is **`E0603` (privacy violation)** for an
external caller, which I confirmed fires. Fix the sentence so the next reader does not inherit a
false model of what the validation proves.

**L2 — The plan's CHANGELOG justification is factually wrong about repo convention.**
The plan skips the CHANGELOG as *"consistent with the repo treating only behaviour/message changes
as entries"*. That is not this repo's convention. The **Unreleased** section (CHANGELOG.md lines
4–500) already carries `#### Infrastructure (contributor-facing)` at line 481 — CI-job and
`justfile` changes with no runtime effect — and history has `#### Tests`, `#### Documentation`, and
even a literal `#### Infrastructure (contributor-facing, no runtime effect)` (line 1242). So the
repo *does* record non-user-visible changes under a dedicated heading.

Whether *this* change earns an entry is a separate judgement, and I think **skipping is still
defensible** — the existing Infrastructure entries are things a contributor trips over (a CI job,
a broken recipe), and a one-keyword visibility demotion is below that bar. But see Optional-1: the
crates.io angle is the one argument that could tip it.

**L3 — The change resolves the confusion for API readers, but only *hides* it for the audience
that actually made the #389 error.** This is the substantive finding, and it answers the brief's
§3 directly. Current doc blocks:

- `norm_path` (io.rs:38-44) — carries the case-folding **trade-off note** and says it is *"what
  every collision check uses"*. Mentions `collision_key`. **Does not mention `path_identity_key`.**
- `path_identity_key` (io.rs:72-74) — *"Case-**sensitive**, because … `X` and `x` are two files"*.
- `collision_key` (io.rs:79-80) — *"Case-**folded**, because outputs differing only in case alias
  each other on APFS/NTFS"*. **Carries no trade-off note.**

The re-audit's concern was that `norm_path`'s block never mentions the case-preserving sibling.
Demoting it to private genuinely resolves that **for anyone reading the public API** (docs.rs or
`cargo doc`): they now land on `path_identity_key` and `collision_key`, which sit six lines apart
and explicitly contrast case-sensitive vs case-folded with the reason for each. That is a good
adjacent pairing and a real improvement.

It does **not** resolve it for the in-file reader — which is precisely who made the #389 mistake.
`norm_path` remains at line 45, remains the *first* fold-related item in the file, remains
documented, and still says nothing about `path_identity_key` while claiming to be what "every
collision check uses" — an invitation to generalise to "every path comparison", which is the #389
error itself. **The keyword is not the fix for that; the doc content is.** Treating #399 as closing
the re-audit's doc concern is therefore only half-right. See Important-1 and Important-2.

---

## Assumptions

- **A1 (zero external callers) — verified true**, by exhaustive grep *and* by a real compile.
  Stronger than the plan claims: it holds across `src/`, `tests/`, `examples/`, `build.rs` and all
  markdown/docs.
- **A2 (`mod tests` reaches private items) — verified true** by actually running the test.
- **A3 (unstated) — "the published crate's library surface is not a contract."** The plan's Context
  bullet 3 asserts *"The crate is a binary product; `lib.rs` exists for unit-test reach, not as a
  public API contract, so no external-consumer concern."* The **stance** is reasonable; the
  **premise is asserted without acknowledging that the crate is published to crates.io**:
  `release.yml:421-441` has a "Publish to crates.io" job, `Cargo.toml` has no `publish = false`,
  no `[lib]` restriction, and full registry metadata (keywords, categories, readme). So the
  implicit lib target ships, docs.rs renders `trim_galore::io::*`, and removing a `pub fn` from a
  2.x crate is a technically semver-breaking API removal.

  I am **not** arguing against the change — nothing sane depends on TrimGalore as a library, and
  this is exactly the confusion the issue wants gone. I am arguing the plan should *say so
  accurately* rather than deny that a published surface exists. See Optional-2 for the durable fix.
- **A4 (unstated) — no behavioural coupling to the fold's name.** True: `--passthrough` validation
  (`cli.rs:934-952`) and the pre-flight (`io.rs:99,104`) both reach the fold through
  `collision_key`; neither names `norm_path`. Confirmed by grep.

---

## Efficiency

Nil, correctly. No call sites move, no allocation or complexity change; `fn` vs `pub fn` is a
visibility annotation with no codegen consequence (and the function was never a cross-crate
inlining candidate, since nothing outside the crate called it). Nothing to add.

---

## Validation sufficiency

**Sufficient for the highest-risk failure mode**, and I proved that rather than assuming it.

- The single risk is a hidden external caller. `cargo clippy --all-targets` compiles lib, bin,
  every `tests/*.rs`, and the example — the complete set of compiled callers. An external caller
  produces `E0603`, which I demonstrated fires. There is no silent-wrong-result mode available to a
  visibility change, so there is nothing for a behavioural test to catch that the compiler misses.
- Validations 2 and 3 (`cargo test`) are cheap and appropriate; I ran validation 2 on the private
  variant and it passes.

Two gaps, both minor:

- **V-gap 1 — nothing exercises rustdoc.** Since CI has no `cargo doc` job, a
  `private_intra_doc_links` warning would never be caught. Irrelevant as the plan stands (no such
  links exist), but it becomes relevant the moment the implementer adds a cross-reference per
  Important-2 — so add the constraint: **use plain backticks, matching this file's existing style,
  not `[…]` intra-doc links.** Optionally run `cargo doc --no-deps` once, locally.
- **V-gap 2 — local `--all-targets` and CI's target set differ.** This working tree has an
  **untracked** `examples/fastqc_only.rs`, so a local `--all-targets` run compiles a file CI does
  not have. It does not reference `norm_path`, so it changes nothing here — but a green local run
  is not automatically a statement about the tracked tree.

---

## Alternatives

| Option | Assessment |
|---|---|
| **A. `fn` (plan's choice)** | **Recommended.** True minimum, zero behavioural surface, keeps the named seam. |
| **B. Inline into `collision_key`** (the issue's alternative) | Correctly rejected, and the plan's reason understates its own case. `collision_key` absolutises, so `collision_key(Path::new("foo.fq.gz"))` returns a CWD-dependent absolute path — the existing exact-string assertions (`== "foo.fq.gz"`, `== "/some/dir/sample_r1.fq.gz"`) could not survive the rewrite, only weaker fold-equality pairs. Worth noting the *behavioural* pairing is already pinned independently by `identity_key_is_case_sensitive_but_collision_key_is_not` (io.rs:1178), so inlining would lose the fold *primitive's* unit test, not the fold's contract. Still: more churn, no benefit. |
| **C. `pub(crate) fn`** | Not mentioned by the plan; worth one line. It removes the item from the published/docs.rs surface *exactly as `fn` does*, while permitting a future second in-crate caller without re-widening. Adds nothing today (no other module wants the raw fold), so `fn` is right — but this is the natural fallback if a second in-crate caller ever appears, and it satisfies #399 equally. |
| **D. Keep `pub`, fix only the docs** | Rejects the issue's premise, so no. But note the asymmetry it exposes: per L3, **the doc edit carries most of the durable value and the keyword carries the API-surface value.** Shipping the keyword alone delivers the lesser half. |

I agree with the plan's choice of A. My disagreement is narrow: the plan says the trade-off note
*"lives well on the named fn"*, and I think that is wrong on the merits — see Important-1.

---

## Action items

### Critical
None.

### Important

**Important-1 — Move the case-folding trade-off note from `norm_path` to `collision_key`.**
Not merely because visibility changes who reads it, but because **`collision_key` is where it
belongs on substance**: `norm_path` is a pure string fold that decides nothing, so it cannot
false-positive anything. The false-positive risk and the "loud early error rather than silent data
loss" justification are properties of *the collision decision*. Two independent reasons to move it:

1. Post-demotion it vanishes from the published API docs (docs.rs renders `collision_key`, which
   currently states only *what* and *why folded*, never the *cost*).
2. Even for in-file readers the note is 35+ lines above `collision_key`, with `lexical_normalise`
   in between, so someone landing at line 81 does not see it today either.

Leave `norm_path` with a one-line "what it does". This turns the change from behaviour-neutral into
a net documentation improvement, at no extra risk.

**Important-2 — Add the missing `path_identity_key` cross-reference; do not rely on the keyword to
resolve the re-audit's concern.** Per L3, privacy hides `norm_path` from API readers but leaves the
in-file reader with the same under-specified block that enabled #389. The plan's outline currently
adds only *"(private)"* framing, which addresses none of this. `norm_path`'s block should note that
the fold is the **output-collision** key and that input identity deliberately uses the
case-**preserving** `path_identity_key` — the two answer different questions. Keep it to one line
per the repo's comment convention, and **use plain backticks, not `[…]` intra-doc links** (V-gap 1:
no CI job would catch a broken/private link).

Also worth the implementer's glance while in the area: the *test's* comment at io.rs:766-768 says
`norm_path` is *"the case-folded normalisation that the output-collision pre-flight + `--passthrough`
validation both use."* True only indirectly — neither calls `norm_path`, both go via `collision_key`
— and it omits that the `--passthrough` check *also* uses the case-sensitive key for a different
question (`cli.rs:934-952` uses both). Same class of imprecision #399 exists to remove. Optional to
fix, three lines from an edit the plan is already making.

**Important-3 — Correct the two wrong justifications (L1, L2) in the plan text before
implementing.** Both are one-sentence fixes, and both currently teach a false model: that
`dead_code` is the safety net (it cannot fire), and that this repo does not changelog
non-user-visible work (it has a dedicated section for it). Neither changes the diff.

### Optional

**Optional-1 — Reconsider the CHANGELOG once, on the correct grounds.** The real argument for an
entry is not "internal change" but "**removes an item from the published crate's API**" (A3). If
Felix wants the crates.io surface tracked, a one-liner under the existing
`#### Infrastructure (contributor-facing)` heading is the right home. My own read: skip it — below
the bar set by the neighbouring entries — but skip it for that reason, not the plan's.

**Optional-2 — The durable fix for this whole class: declare the library surface unstable.** All 17
modules in `lib.rs` are `pub` with no `pub use` curation, so every internal helper in the crate is
published API and every future demotion re-raises A3. A single `//!` note in `src/lib.rs` — the
library surface is an implementation detail of the binary and carries no semver guarantee — would
pre-authorise this change and all its successors. Out of scope for #399; worth its own issue.

**Optional-3 — Note the untracked-`examples/` caveat** (V-gap 2) if the implementer reports a green
local `--all-targets` run as evidence about the tracked tree.

---

## Self-assessment of this review

The load-bearing negative ("no callers outside `src/io.rs`") was confirmed three ways: exhaustive
grep, a positive control proving the grep could see markdown, and a real compile with a positive
control proving the compiler rejects an external caller. The plan's mechanism is correct and safe;
my findings are about the plan's *reasoning* (two wrong justifications) and about an opportunity it
declines (the doc content, which is where the #389 confusion actually lives). Nothing here should
delay implementation.
