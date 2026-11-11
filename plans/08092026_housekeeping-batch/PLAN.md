# Plan: housekeeping batch — #401, #402, #405, #407 in one PR

**Issues:** [#401](https://github.com/FelixKrueger/TrimGalore/issues/401), [#402](https://github.com/FelixKrueger/TrimGalore/issues/402), [#405](https://github.com/FelixKrueger/TrimGalore/issues/405), [#407](https://github.com/FelixKrueger/TrimGalore/issues/407)
**Batching decision (Felix, 2026-08-09):** one PR. All four are mechanical, none change behaviour, and none interact — four separate plan→review→implement pipelines for ~30 lines is where the process stops paying for itself.
**Base:** `dev` @ `90dda11`.

## Goal

Clear four accumulated documentation/test-hygiene items in a single reviewed round, taking the open queue from 7 to 3. No runtime behaviour changes anywhere in this batch. One item (#405) turns out **not** to be purely internal — see Part C.

## Part A — #401: the #385 loose ends

Two items, both verified against the tree.

**A1. `CHANGELOG.md:206-211` contradicts what shipped.** The entry ends *"All four now use the same key as the pre-flight."* D13 then split the key in two, and #389 refined the `--passthrough` half again. Today: the pre-flight uses the case-folded `collision_key`; the three input-identity checks use the case-preserving `path_identity_key`; `--passthrough` uses both. Nothing in the CHANGELOG discloses the split (grep for `path_identity_key` / "identity key" → zero hits). The *bug* the entry describes is genuinely fixed (both keys share `lexical_normalise`), so the entry is **stale, not false about the outcome** — correct the final sentence only, and do not rewrite the history around it.

**A2. The D8→D13 reversal has no unit guard.** `CODE_REVIEW_D13.md` finding 3 asked for a four-line test that was never added:

```rust
let cli = Cli::parse_from(["trim_galore", "--paired",
    "Sample_R1.fastq.gz", "Sample_R2.fastq.gz", "SAMPLE_R1.fastq.gz", "SAMPLE_R2.fastq.gz"]);
assert!(cli.validate().is_ok(), "case-only variants are distinct inputs");
```

D8 over-applied case-folding to the input-identity checks, which rejected two genuinely distinct files on a case-sensitive filesystem; D13 reversed it. The `io.rs` tests pin the *keys*; what broke was the wiring in `cli.rs`, and the only thing covering that today is the #216 CI validation job (`.github/workflows/ci.yml:625-630`). No files on disk, so it runs identically on macOS and Linux. Goes in `cli.rs`'s test module near the other `validate()` tests.

## Part B — #402: declare the library surface an implementation detail

`src/lib.rs` is 17 bare `pub mod` lines with no `pub use` curation and no module doc at all (verified: the file opens directly with `pub mod adapter;`). The crate publishes to crates.io, so docs.rs renders `trim_galore::io::*` and every future visibility tightening — #399's `norm_path` demotion was the first — is technically a semver-breaking API removal.

Add a `//!` header stating that the library target exists so integration tests and the binary can reach internals, that its surface is an implementation detail, and that it carries no semver guarantee. Two or three lines. This pre-authorises the whole class rather than re-litigating it per demotion.

## Part C — #405: the 12 rustdoc errors (NOT purely internal)

`RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` fails on `dev` with 12 diagnostics. Exact locations, verified:

| Site | Diagnostic | Renders in `--help`? |
|---|---|---|
| `src/alignment.rs:11` | public doc for `alignment` links to **private** item `myers_proves_no_match` | no (module doc) |
| `src/cli.rs:412,416,420,424,428,432` | unresolved link to `Deprecated` (×6) | **no** — all six carry `hide = true`; verified `--help` contains zero "Deprecated" |
| `src/clump_only.rs:874` | unresolved link to `i` (×2) — from `R1[i], R2[i]` in a doc comment | no (fn doc) |
| `src/report.rs:33` | unresolved link to `match_len` — from `adapter_length_counts[adapter_idx][match_len]` | no (field doc) |
| `src/cli.rs:381,386` | unclosed HTML tag `N` — from `.<N>bp_5prime.fq(.gz)` | **YES** |

**The two `<N>` sites are the reason this part needs a decision rather than a sweep.** They are clap doc comments, so they double as `--help` text, and `--help` currently renders `.<N>bp_5prime.fq(.gz)` literally. Any rustdoc fix changes user-facing output:

- **Backticks** — `` `.<N>bp_5prime.fq(.gz)` `` — silences rustdoc and renders with visible backticks in `--help`. **Consistent with house style**: help text throughout `cli.rs` already uses backticks freely (`--memory`'s help renders `` `1G` ``, `` `512M` ``). **Recommended.**
- Drop the angle brackets (`.Nbp_5prime.fq(.gz)`) — neutral in both renderings, but slightly less clear that `N` is a placeholder.
- Escape (`\<N\>`) — silences rustdoc, but leaks backslashes into `--help`. Rejected.

The other ten are genuinely invisible to users: backtick the bare identifiers (`` `i` ``, `` `match_len` ``, `` `[Deprecated]` `` → `` `Deprecated` `` or plain text), and for `alignment.rs:11` either backtick `myers_proves_no_match` (dropping the link) or note that the item is private. Backticking is the minimal, consistent choice — it matches what `io.rs` does throughout.

**Deliberately NOT in scope:** adding a `cargo doc` CI gate. #405 raises it as a second, separable part, and adding a new CI job is a bigger decision than this batch should carry — without it, though, the next broken link accumulates silently, which is exactly how these twelve arrived. Recommend filing that as its own issue if Felix wants it; noted under Questions.

## Part D — #407: `--preserve-tags` help contradicts the code

`src/cli.rs:238` says *"Ignored for FASTQ input."* But `main.rs:286-296` makes it a **hard error** when all inputs are FASTQ *and* `--output-format ubam` is set ("--preserve-tags has no effect with all-FASTQ inputs; either remove the flag or convert at least one input to uBAM via 'samtools import'"). Without `--output-format ubam` it really is just a warning — so the text is half right, which is what makes it misleading: "ignored" gives no reason to expect an abort.

One clause: `Ignored for FASTQ input (a hard error when combined with --output-format ubam, since there are no source tags to propagate).`

## Behavior

No runtime behaviour changes in any part. What changes: one CHANGELOG sentence, one new unit test, a `lib.rs` module doc, 12 doc-comment fixes (two of which alter `--help` text by adding backticks), and one help-text clause.

## Implementation outline

1. **A1** — correct `CHANGELOG.md:206-211`'s final sentence to describe the two keys and why they differ (case-folding is right for output paths, wrong for input identity). Do not touch the surrounding historical entries.
2. **A2** — add the four-line test to `cli.rs`'s test module, named for what it guards (e.g. `test_case_only_input_variants_are_distinct`), with a comment naming D8/D13 and #383 so the next reader knows why it exists.
3. **B** — add the `//!` header to `src/lib.rs`.
4. **C** — the ten invisible sites: backtick the bare identifiers. The two `<N>` sites: backticks per the recommendation above. Then re-run `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` and confirm it exits 0 — the count must go to zero, not merely down.
5. **D** — the one clause at `cli.rs:238`.
6. **CHANGELOG for the batch itself** — one bullet under `#### Changes` for the `--preserve-tags` help correction (user-facing) and the `<N>` backticks (user-facing). A1 is itself a CHANGELOG edit; B and the other ten C items are internal and need no entry (per #399's established reasoning: below the `#### Infrastructure` granularity floor). Avoid the `#### Bug fixes` / `#### Fixes` duplication hazard (#401 itself documents it).
7. `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test`, and the rustdoc run from step 4.

## Efficiency / Integration

Nil throughout — docs, one test, one changelog sentence. No call sites, no runtime paths, no output bytes.

## Assumptions

- **A1:** the six `Deprecated` sites are invisible to users — verified via `hide = true` on all six and a `--help` grep returning zero.
- **A2:** the two `<N>` sites *are* visible — verified by rendering `--help`.
- **A3:** nothing pins the `--preserve-tags` help text or the `<N>` strings. **To be verified at implementation** (grep `tests/`, `.github/`); `integration_no_args_help.rs` is known to pin only `"Usage: trim_galore"` and `"--adapter"`, but the `<N>` strings have not been grepped yet and this batch must not repeat #400's deferred-grep mistake.
- **A4:** the new A2 test passes on the current tree — i.e. it documents behaviour that is already correct post-D13, rather than exposing a live bug. **Must be confirmed before committing**; if it fails, that is a genuine regression and this batch stops and reports rather than adjusting the test.
- **A5:** `cargo doc` exiting 0 after part C means all 12 are fixed, not that some became warnings — the run uses `-D warnings`, so this holds by construction.

## Validation

| # | What | How | Expected |
|---|------|-----|----------|
| 1 | Rustdoc clean | `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` | exit 0, zero diagnostics (was 12) |
| 2 | A2 test guards the right thing | run it, then **invert the guard locally** (re-key one identity check on `collision_key`) and confirm it fails | passes as shipped; fails under the D8 mutation |
| 3 | `--help` still correct after the `<N>` change | rebuild, render `--help`, read both hardtrim entries | placeholder still legible; no stray escapes |
| 4 | Nothing pinned was broken | grep `tests/ .github/` for the changed help substrings (A3) before editing; full `cargo test` after | zero hits pre-edit; suite green after |
| 5 | A1 is accurate | read the corrected CHANGELOG sentence against `io.rs`'s two key functions and `cli.rs`'s four checks | describes the split truthfully |
| 6 | No collateral | fmt + clippy `-D warnings` + full `cargo test` | green |

## Questions or ambiguities

- **[Open — Felix's call] The `<N>` fix is user-visible.** Backticks are recommended (house style), but it does change `--help` output for two flags. Say if you'd rather drop the angle brackets or leave those two rustdoc errors alone.
- **[Open] A `cargo doc` CI gate** is deliberately out of scope (see Part C). Worth its own issue if wanted — without it these accumulate silently.
- **[Open] Batch atomicity.** If any part turns out non-trivial under review, I'd rather split it out than grow the batch. Flagging so that's the expected response, not a surprise.

## Self-Review

- **What the investigation changed:** #405 looked like a pure backtick sweep and isn't — two of its twelve sites render in `--help`, so the batch carries a user-visible change and needed a recommendation rather than a default. That is the one place this plan could have quietly shipped a `--help` regression.
- **Traps checked:** A3 does **not** defer the pinning grep as a hand-wave — it is called out as unverified and assigned to implementation, because deferring exactly this grep is what #400's review flagged. A4 makes the new test's pass/fail meaning explicit up front, so a failure is treated as a regression discovery rather than a test to tweak. Validation 2 is an expected-fail control on the one new test, per the session's repeated finding that a check which cannot fail proves nothing.
- **Edge cases:** none material — no branching logic is added anywhere in the batch.
- **Remaining risk:** only the `<N>` register question, which is Felix's to settle.
