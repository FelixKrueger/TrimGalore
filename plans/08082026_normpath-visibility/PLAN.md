# Plan: demote `io::norm_path` to private, and re-home the trade-off note (#399)

**Issue:** [#399](https://github.com/FelixKrueger/TrimGalore/issues/399) — filed from #389's review round.
**Base:** `dev` @ `ce13761`.

## Revision history

- **r2 (2026-08-08):** incorporates dual plan review (PLAN_REVIEW_A.md, PLAN_REVIEW_B.md — both *approve, no Criticals*). Both reviewers independently compiled the demotion in isolated copies of the crate with positive controls, and both raised the same four corrections: r1's validation credited the wrong command, its dead-code parenthetical was backwards, its CHANGELOG rationale misdescribed repo convention, and its reason for keeping the named function was the weaker of the two available. Scope gains **one doc hunk** (Important in both reviews): move the case-folding trade-off note to `collision_key`. Behaviour still unchanged.
- r1 (2026-08-08): initial plan — visibility keyword + a "(private)" doc tweak.

## Goal

`io::norm_path` is `pub` with no callers outside `src/io.rs`. The public name is what enabled #389's planning error (it was mistaken for a stale alias of `collision_key`). Demote it to module-private, and — because privatization removes it from rendered docs — move the substantive case-folding trade-off note onto `collision_key`, where it belongs on substance and where public readers will land.

## Context

- `io.rs:45-47` — the function. Sole production caller `io.rs:82` (`collision_key`); test at `io.rs:765-780`. `cli.rs:1903` is a `//` comment naming the *test*, unaffected.
- **Callerless-outside-`io.rs` verified four ways across both reviews**, including the surfaces a naive grep misses: `src/lib.rs` has 17 bare `pub mod` and **no `pub use`**; `benches/`/`examples/` don't exist in the tracked tree; `build.rs`, `Cargo.toml`, `docs/`, and untracked non-`plans/` files are clean. **B's sharpest catch:** `specialty.rs:12` and `clump_only.rs:52` both `use crate::io as naming`, so a caller there would read `naming::norm_path` — and historical plan text *does* describe `naming::norm_path` call sites, making this a live possibility rather than a hypothetical. Grepping the bare symbol covers the alias; grepping the qualified path would not have.
- **Both reviewers compiled it for real** (isolated `git archive` copies; this tree untouched): `cargo check/clippy --all-targets` clean, 0 warnings, all 11 integration tests genuinely compiled (`.rmeta` present).
- **The safety net is `E0603`, not `dead_code`, and only on `--all-targets`.** B ran the positive control: with an external caller planted in `tests/`, `cargo check` exits **0** (it compiles lib + bin only) while `cargo check --all-targets` exits **101** with `error[E0603]: function norm_path is private`; removing the probe returned it to green, so the harness demonstrably both passes and fails.
- **No rustdoc gate exists in CI** (both verified): `ci.yml:160` "Docs build (Astro)" is the documentation *website*; no `cargo doc`, `RUSTDOCFLAGS`, or `--document-private-items` in any workflow. So the private-intra-doc-link hazard is both absent (no bracketed `[…]` references to the symbol anywhere — the file uses plain backticks throughout) and ungated.
- **Doc-comment layout today:** `norm_path` (`:38-44`) carries the substantive *"Pragmatic trade-off: on opt-in case-sensitive APFS volumes this may false-positive, but the penalty is a loud early error rather than silent data loss"*; `path_identity_key` (`:72-77`) explains its case-sensitive rationale; `collision_key` (`:79-83`) is two lines and carries **none** of the trade-off.
- **The crate publishes a library surface.** `release.yml:421-441` has a crates.io publish job, `Cargo.toml` has no `publish = false` and no `[lib]` restriction, so docs.rs renders `trim_galore::io::*` and removing a `pub fn` from a 2.x crate is technically a semver-breaking API removal (A §A3). r1 asserted the surface "is not a public API contract"; that stance is reasonable — nothing sane depends on TrimGalore as a library, and this is precisely the confusion #399 targets — but the plan should state it accurately rather than deny the surface exists. Durable fix filed separately (see Follow-ups).

## Behavior

No behavioural change: same fold, same caller, same tests, same rendered output. What changes: one visibility keyword, and which doc block carries the trade-off paragraph.

## Implementation outline

1. `io.rs:45` — `pub fn norm_path` → `fn norm_path`.
2. **Move the trade-off paragraph** from `norm_path` (`:42-44`) onto `collision_key` (`:79-83`). Rationale (both reviewers, independently): `norm_path` is a pure `to_string_lossy().to_ascii_lowercase()` — "may false-positive" is not a property of lowercasing a string, it is a property of *deciding output collisions on a folded key*, which happens in `collision_key`. Privatization would otherwise delete the only statement of the trade-off from the rendered public docs, and an in-source reader arriving at `collision_key` must currently scroll *up* past `lexical_normalise` and `path_identity_key` to find it.
3. **Reduce `norm_path`'s doc to one line** — e.g. `/// Case-folded (ASCII lowercase) view of a path; the folding half of `collision_key`.` **Do not** add "(private)" framing (r1's step 1, dropped per B §1.7): the `fn` keyword three lines below states the visibility, and docs that narrate their own visibility go stale the next time visibility changes.
4. While in the block, add the one-line note that input identity deliberately uses the case-**preserving** `path_identity_key` — the two keys answer different questions (A Important-2). **Plain backticks, not `[…]` intra-doc links** (no CI job would catch a broken or private link).
5. Optional drive-by, three lines from an edit already being made: the *test's* comment at `io.rs:766-768` says `norm_path` is what "the output-collision pre-flight + `--passthrough` validation both use" — true only indirectly (both go via `collision_key`), and it omits that the `--passthrough` check also uses the case-sensitive key for a different question (`cli.rs:934-952` uses both). Same class of imprecision #399 exists to remove.
6. **Verification:** `cargo clippy --all-targets --release -- -D warnings`, `cargo test test_norm_path_case_folds`, full `cargo test`, `cargo fmt --all -- --check` (step 3 shortens a doc line and could disturb rustfmt's wrapping).
7. **No CHANGELOG entry** — decision unchanged from r1, rationale corrected. The repo *does* have a slot for non-user-visible work (`CHANGELOG.md:481`, `#### Infrastructure (contributor-facing)`, with CI-job and `justfile` entries), so r1's "the repo only records behaviour/message changes" was false. The real reasons: this sits below that section's granularity floor (its entries are things a contributor trips over), and **`grep` for `collision_key|norm_path|path_identity_key` across 1300+ CHANGELOG lines returns zero hits** — the changelog has never named these helpers across the entire #216/#383/#384/#388/#389/#391 arc, so it should not start for a visibility keyword. Supporting: `7ea2741` (docs-only, #387) merged with no entry.

## Efficiency / Integration

Nil / none. Visibility is an annotation with no codegen consequence; the sole call is already intra-crate.

## Assumptions

- **A1:** zero callers outside `io.rs` — verified by exhaustive grep (including the `naming::` aliases), by real compiles, and by a positive control proving the greps could see markdown.
- **A2:** `mod tests` reaches private items — verified by *running* `test_norm_path_case_folds` on a private-fn copy, not by citing the language rule.
- **A3:** the published library surface is uncurated and treated as an implementation detail of the binary; this change is technically a semver-breaking API removal that no plausible consumer notices. Stated rather than denied (A §A3).
- **A4:** no behavioural coupling to the fold's *name* — `--passthrough` validation and the pre-flight both reach the fold via `collision_key`; neither names `norm_path`.
- **A5:** no doc-tests reference the symbol. Structurally impossible to matter here (the crate's only doc fence is a ```` ```text ```` block at `specialty.rs:235-240`), but worth recording *why* it matters: `--all-targets` excludes doc-tests, so a doc example using `norm_path` would slip past clippy and surface only under `cargo test` (B §2).

## Validation

| # | What | How | Expected |
|---|------|-----|----------|
| 1 | No missed caller | `cargo clippy --all-targets --release -- -D warnings` — **this** is the net, not `cargo build`: `check`/`build` compile lib + bin only and miss an external caller in `tests/` entirely (B §3.1, proven by positive control) | clean; an external caller would be `E0603` |
| 2 | Fold still pinned | `cargo test test_norm_path_case_folds` | passes (child module reaches the private fn) |
| 3 | No collateral | full `cargo test` + `cargo fmt --all -- --check` | green |
| 4 | Trade-off note landed where readers arrive | read `collision_key`'s block after the move | states both *what* and the false-positive cost; `norm_path` reduced to one line |

**Caveat for the implementer** (A V-gap 2): this working tree has an **untracked** `examples/fastqc_only.rs`, so a local `--all-targets` run compiles a file CI does not have. It doesn't reference `norm_path`, but a green local run is not automatically a statement about the tracked tree.

## Questions or ambiguities

- **[Resolved r2 — reviewer contradiction, settled without escalation]** A offered `pub(crate) fn` as an equally-valid variant ("removes it from the published surface exactly as `fn` does, permits a future second in-crate caller"); B called it *strictly worse here* because it keeps the symbol reachable from `specialty.rs`/`clump_only.rs` via their existing `crate::io as naming` aliases — preserving most of the confusion surface #399 exists to remove, for no benefit, since no cross-module caller exists or is wanted. **B's argument wins, and it is consistent with A's own reasoning** that the in-file/in-crate reader is the audience who made the #389 error. Plan keeps plain `fn`; the commit message should note `pub(crate)` was considered and rejected, since it is the obvious reviewer question.
- **[Open, non-critical]** Exact doc copy is reviewable prose.

## Follow-ups (out of scope, filed separately)

- All 17 modules in `lib.rs` are `pub` with no `pub use` curation, so every internal helper is published API and every future demotion re-raises A3. A single `//!` note declaring the library surface an implementation detail with no semver guarantee would pre-authorise this change and its successors (A Optional-2).

## Implementation notes (2026-08-08, branch `fix/399-normpath-private`)

- `315557d` — demotion + doc re-homing as planned (all 7 outline steps, incl. the optional step-5 test-comment fix).
- `1a2a384` — review batch. Both code reviewers APPROVED with no Criticals and independently raised the same Medium: the relocated trade-off said "opt-in case-sensitive APFS volumes", naming the rarest instance of its own condition, when Linux ext4/xfs are case-sensitive by default and every CI run exercises the rejection there (`io.rs:1170-1174` already said "which on Linux is false"). Widened. Also took both reviewers' L1/L4: the test comment claimed callers never reach the fold directly, five lines above five direct calls — now scoped to production callers, and it recovers B's "`collision_key` absolutises first" fact, which is what explains why the fold is tested directly.
- **Reviewer contradiction resolved:** A suggested making the `path_identity_key` reference a real intra-doc link (its target is public, so no private-link risk); B leaned to plain backticks for local consistency (`io.rs` has zero bracketed links). Kept backticks — B's consistency argument, and A explicitly called its own point the weaker of the two.
- **Not taken:** A's M2 and both reviewers' rustdoc observations → filed as #405 (12 pre-existing `cargo doc` errors, incl. `alignment.rs:11` — a public doc linking to a private item, the realized form of the hazard this plan reasoned about hypothetically, and completely ungated since no workflow runs `cargo doc`).
- **Coverage: COMPLETE** (19/19, zero gaps).
- **Verification note:** I ran the `E0603` positive control on the branch itself, not only in the reviewers' isolated copies — `cargo check` passes with an external caller planted in `tests/` while `--all-targets` exits 101, confirming the plan's corrected claim about which command is the guard.

## Self-Review (r2)

- **What r1 got wrong, owned:** three of its justifications were false or weak — `dead_code` cannot fire (its own parenthetical said why), `cargo build` is not the guard (it never compiles `tests/`, which is exactly where a hidden consumer of a `pub` item would live), and the repo does changelog non-user-visible work under a dedicated heading. r1's reason for keeping the named function ("the trade-off note lives well on the named fn") is now *inverted* by r2: the note moves precisely because it does **not** belong there. The surviving reason to keep the named seam is B §1.4's: `collision_key` absolutises, so re-expressing the fold test through it would make its relative-path assertions cwd-dependent.
- **What the reviews added beyond corrections:** the doc re-homing, which converts a behaviour-neutral hide into "put each fact where its reader is" — and which is the difference between #399 *resolving* the re-audit's doc finding and merely concealing it. The misuse half of that finding is genuinely resolved and strengthened (from "a doc a reader might skim" to `E0603`); the informational half is resolved only if the note moves.
- **Remaining risks:** none structural. A visibility change has no silent-wrong-result mode; the single failure mode is a compile error the validation demonstrably catches.
