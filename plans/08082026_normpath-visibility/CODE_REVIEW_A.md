# Code Review A — #399 `norm_path` visibility + doc re-homing

**Reviewer:** A (independent; Reviewer B reviews the same diff separately)
**Branch:** `fix/399-normpath-private` @ `315557d` — base `dev` @ `ac08a97`
**Diff:** `git diff ac08a97..315557d` — one file, `src/io.rs`, 3 hunks, +10/−12
**Plan:** `plans/08082026_normpath-visibility/PLAN.md` (r2)
**Verdict:** **APPROVE.** No Critical, no High. One Medium worth fixing before merge (a doc statement that is narrower than the truth), one Medium that belongs in a separate issue, five Lows.

---

## Summary

The change does what it claims and nothing else. `io::norm_path` is now module-private, the substantive case-folding trade-off has moved onto `collision_key` where the decision is actually made, and `collision_key` gained an accurate pointer to the case-preserving `path_identity_key`. Behaviour is unchanged; I verified that by compiling and running, not by reading.

Two things are worth the caller's attention.

**First, the one real defect.** The relocated trade-off paragraph says the false positive happens "on opt-in case-sensitive APFS volumes". That scoping is too narrow: Linux ext4/xfs are case-sensitive *by default*, so ordinary Linux — the platform most Trim Galore runs on, and the platform CI runs on — is the primary home of this false positive, not a macOS opt-in edge case. The repo already knows this: the regression-test doc at `io.rs:1170-1174` uses the general framing and ends "which on Linux is false". The text was moved verbatim, so this is not a regression introduced here, but the move promotes it to the primary public location for the fact, and #399 exists to remove exactly this class of imprecision. One-line generalisation, trivial fix.

**Second, a caveat about how I verified.** `cargo clippy --all-targets --release -- -D warnings` in the shared working tree returned exit 0 in **0.32 s** — a pure cache hit, not a compile. I did not accept that as evidence. I re-ran the check from scratch in an isolated `git archive` copy of `315557d` (tracked files only, so no untracked `examples/fastqc_only.rs`, matching what CI sees) with a private target dir, and additionally ran the positive control myself. Results below. The negative is trustworthy; it would not have been on the cached run alone.

The plan's central question — "does privatizing delete the only public statement of the trade-off?" — is verifiably answered *no*: I built the rustdoc and confirmed `norm_path` has vanished from the rendered surface entirely while `collision_key`'s page now carries the paragraph.

---

## Verification log

Everything here I ran or read myself.

| # | Check | Result |
|---|-------|--------|
| 1 | `cargo fmt --all -- --check` | **clean** |
| 2 | `cargo clippy --all-targets --release -- -D warnings` (shared tree) | exit 0 — but `Finished in 0.32s`, i.e. **cache hit, not a compile**. Not relied upon. |
| 3 | **Cold-cache** `cargo clippy --all-targets -- -D warnings`, isolated `git archive 315557d` copy, private target dir, tracked files only | **exit 0, clean** — a genuine compile of exactly what CI compiles |
| 4 | **Positive control:** external caller `trim_galore::io::norm_path` planted in `tests/` of the isolated copy | **`error[E0603]: function norm_path is private`, exit 101** — the harness demonstrably fails |
| 5 | Probe removed → re-check | **exit 0** — back to green, so #4 was caused by the probe |
| 6 | `cargo test` (full) **on the pinned isolated copy**, with `fn norm_path` (private) confirmed at `io.rs:39` in that tree | **561 passed, 0 failed** across 1 lib + 11 integration binaries + doc-tests, incl. `io::tests::test_norm_path_case_folds ... ok` |
| 7 | `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` | fails on **12 pre-existing** diagnostics, **none from `io.rs`** — this diff adds zero |
| 8 | Rendered docs: is `norm_path` gone? | **yes** — no `fn.norm_path.html`, zero grep hits across the entire generated doc tree |
| 9 | Rendered docs: did the trade-off survive? | **yes** — `io/fn.collision_key.html` renders both the "different question" pointer and the trade-off paragraph |
| 10 | CI rustdoc gate (`cargo doc` / `RUSTDOCFLAGS` / `rustdoc` / `--document-private-items` in `.github/workflows/`) | **absent** (grep exit 1) — confirms the plan |
| 11 | Every production caller of `collision_key` read for the "loud early error" claim | `io.rs:104-110` bails; `cli.rs:949-962` bails. **No caller where a false positive could silently lose data.** |
| 12 | Repo-wide `git grep norm_path -- ':!plans/'` | only `src/io.rs` + one `//` comment naming the *test* |

Two notes on how this evidence was obtained, both of which changed what I was willing to claim.

**The shared working tree moved off the reviewed commit mid-review.** By the time I finished, `/Users/fkrueger/Github/TrimGalore` was on `fix/400-clump-help-text` @ `0676fdf`, and `315557d` is **not** an ancestor of it. `git diff 315557d HEAD -- src/io.rs` is +12/−10 — the exact inverse of #399's −10/+12 — i.e. that branch carries the *pre*-#399 `pub fn norm_path`. Consequences worth knowing:

- Every load-bearing check above (#3, #4, #5, #6, #7, #8, #9) ran against a `git archive 315557d` copy, so it is commit-pinned and unaffected. I re-ran the full test suite on that pinned copy specifically to remove any doubt, after confirming `fn norm_path` is private in it.
- `fix/400-clump-help-text` was cut from `dev` (`ac08a97`) and **does not contain #399**. That is fine for parallel work, but if both land the caller should merge them independently rather than assuming one carries the other. #400 also adds a CHANGELOG entry (+15 lines) where #399 deliberately adds none.
- Any check another agent runs in the shared tree right now is not a statement about #399.

**One of my own runs was a false negative I had to discard.** My first `cargo test` was piped through `tail -40`, which truncated the log to a single binary's summary and hid the fold test entirely — it looked like a 23-test suite. The 561-test figure comes only from runs with the complete log captured.

---

## Issues by area

### 1. Logic / correctness — sound

- **Demotion compiles, nothing reachable is broken.** Verified by #3/#4/#5 above, not by inference. `--all-targets` is genuinely the net: the positive control proves it catches an external caller that plain `cargo check` would miss.
- **The fold is still pinned.** `io::tests::test_norm_path_case_folds` runs and passes, so the child module does reach the private fn (plan assumption A2 verified by execution). The assertions are unchanged; `collision_key` was correctly *not* substituted for `norm_path` in the test, because `collision_key` absolutises and would make the relative-path assertions cwd-dependent.
- **The `path_identity_key` cross-reference is accurate.** It is genuinely case-preserving (`lexical_normalise` then `to_string_lossy`, no folding anywhere), and it is genuinely what input identity uses — `cli.rs:537`, `549-550`, `635`, `638`, `934`, `941`. The sentence is true as written.
- **The relocated trade-off is now attached to the right function.** This is the substantive win. `norm_path` is `to_string_lossy().to_ascii_lowercase()` — "may false-positive" is not a property of lowercasing a string. It is a property of deciding output collisions on a folded key, which is what `collision_key` does. I checked the claim's second half against every production caller (#11): both comparison sites bail loudly, so "the penalty is a loud early error rather than silent data loss" holds universally for this key. There is no caller that feeds `collision_key` into a dedup or map where a false positive would silently discard work.

#### M1 — the trade-off's scope is too narrow *(Medium; trivial fix)*

`src/io.rs:77-79`:

```rust
/// Pragmatic trade-off: on opt-in case-sensitive APFS volumes this may
/// false-positive, but the penalty is a loud early error rather than
/// silent data loss.
```

"Opt-in case-sensitive APFS volumes" names only the macOS case. But:

- Linux ext4/xfs/btrfs are case-sensitive **by default**. On Linux, `Sample_R1_trimmed.fq` and `sample_r1_trimmed.fq` are two genuinely distinct files, and this key folds them together, so the pre-flight refuses a run that would have been safe. That is the false positive, and it lives on the dominant platform — not on an opt-in volume.
- The repo already documents it that way. `io.rs:1170-1174` (the `identity_key_is_case_sensitive_but_collision_key_is_not` doc): the #216 CI guard "feeds four genuinely distinct files (`Sample_R1` / `SAMPLE_R1`) on a case-sensitive filesystem and asserts the *output* pre-flight refuses them … which on Linux is false." So CI exercises this trade-off on Linux on every run.

The rest of the file's "APFS/NTFS" framing is correct where it appears, because there it explains why folding is *necessary* (those filesystems are case-insensitive). The false-positive clause is the mirror image and needs the mirror-image scope. Suggested replacement, same length:

```rust
/// Pragmatic trade-off: on a case-sensitive filesystem (Linux ext4, opt-in APFS)
/// two outputs differing only in case are distinct, and this rejects them — a loud
/// early error in preference to silent data loss.
```

Not a regression from this diff (moved verbatim), but this diff is the moment it becomes the canonical public statement of the fact.

#### L1 — "these two paths" for a one-path function *(Low; pre-existing, no action)*

`io.rs:73` opens "Would these two paths be the same *output* file?" while the signature takes one `&Path` and returns a key; the two-path comparison is the caller's. Reads fine as a rhetorical framing and predates this diff.

### 2. Errors / rustdoc — clean here, but the crate is not

- **No bracketed intra-doc link exists anywhere in `io.rs`.** The only two `[…]` constructs are full markdown URLs to issue #381 (`io.rs:532`, `io.rs:583`). The plan's instruction to use plain backticks rather than intra-doc links was followed, so **the demotion cannot have created a private-link error** — and I confirmed that positively by building the docs: of 12 rustdoc diagnostics, **zero come from `io.rs`**.
- **CI would not catch it either way** — no rustdoc job exists (#10), independently confirming the plan.

#### M2 — the crate already has 12 rustdoc errors, including this exact hazard *(Medium; out of scope — file separately)*

`RUSTDOCFLAGS="-D warnings" cargo doc --no-deps` fails today on `dev`, independent of this change:

| Location | Diagnostic |
|---|---|
| `src/alignment.rs:11` | **public documentation for `alignment` links to private item `myers_proves_no_match`** |
| `src/cli.rs:401,405,409,413,417,421` | unresolved link to `Deprecated` (×6) |
| `src/clump_only.rs:874` | unresolved link to `i` (×2) |
| `src/report.rs:33` | unresolved link to `match_len` |
| `src/cli.rs:370,375` | unclosed HTML tag `N` (from `.<N>bp_5prime.fq`) |

The `alignment.rs:11` entry is the realized instance of precisely the hazard the plan reasoned about hypothetically — a public doc pointing at a private item. Since `release.yml` publishes to crates.io and `Cargo.toml` has no `publish = false`, docs.rs renders all of this. It is not this PR's job to fix, but it is the natural companion to the lib.rs-surface follow-up the plan already filed, and it means a future `[norm_path]` link would join the list silently. Worth an issue.

#### L2 — the `path_identity_key` reference is not clickable *(Low; recommend leaving as-is)*

Rendered output shows `path_identity_key` as plain code, not a link. An intra-doc link would be safe here (the target is `pub`, so no private-link risk) and rustdoc would catch a typo in it. But the file uses backticks uniformly, and consistency is the better argument. Keep the backticks — noting only that the plan's stated reason ("no CI job would catch a broken link") is the weaker justification of the two available.

### 3. Structure / style — compliant

- **The 2-line test comment is within convention.** CLAUDE.md sets one line as default, two as maximum. It states the fact, carries no measurements and no reasoning chain. Compliant.
- **`///` blocks are a different register from `//` comments in this file, and correctly treated as such.** `is_gzipped` carries 21 doc lines, `ensure_output_dir` 17. Against that, `collision_key` at 7 lines (was 2) and `norm_path` at 1 (was 7) are both in keeping. The one-line rule governs code comments, and the diff respects it there.
- **Widths fine.** New lines run 83–88 chars against a pre-existing file maximum of 90; `cargo fmt --check` is clean.

#### L3 — moved paragraph kept its old wrap *(Low; trivial fix)*

`io.rs:77-79` wraps at ~70 chars while lines 73-75 immediately above wrap at 83-88, giving the block a ragged right edge halfway down. Reflowing the paragraph to match (2 lines instead of 3) would tidy it. Folds naturally into the M1 rewrite.

#### L4 — "never directly" sits above five direct calls *(Low; trivial fix)*

`io.rs:765-766` says callers "reach it through that key, never directly", and the next five lines call `norm_path` directly — the test itself is the exception. True of production callers, so `production callers` would remove the wrinkle.

#### L5 — mild garden-path in the test comment *(Low; optional)*

"The fold `collision_key` is built on (…)" elides the relative pronoun, so with the backticks drawing the eye it can momentarily parse as `collision_key` being the subject of "is built on". "The fold that `collision_key` is built on" costs one word.

#### Did anything become less clear by moving?

One thing, mildly. The old test comment named *where the fold matters* ("the output-collision pre-flight + `--passthrough` validation both use"); the new one names only the immediate consumer. A reader of the test now learns the plumbing but not the stakes. Two things offset it: the old claim was true only indirectly (both reach the fold via `collision_key`), which is the imprecision #399 exists to remove; and the stakes now live on `collision_key`, which is the better home. **Net clarity: improved.** No other regression — the in-source reader who previously found the trade-off at the fold now finds it 35 lines later at the decision point, which is where they need it.

### 4. Dangling references — none

- **No `pub` item's doc mentions `norm_path`.** Checked every doc block in `io.rs` and grepped `src/` — nothing points a public reader at the now-private symbol for the trade-off or for anything else. This was the specific risk of moving the paragraph out, and it is clean.
- **`cli.rs:1903` is live, not dangling.** It reads "fold pinned by `io::tests::test_norm_path_case_folds`" — it names the *test*, which still exists under exactly that name (it appears in the passing test list), and it is a `//` comment, not a doc comment. Correctly left alone.
- **Nothing outside `src/` refers to the symbol.** `git grep norm_path -- ':!plans/'` returns `src/` only: no `docs/`, no `build.rs`, no `Cargo.toml`, no README. The untracked `examples/fastqc_only.rs` does not reference it either (and my isolated verification excluded it, so the green result speaks about the tracked tree).
- **Rendered docs corroborate.** Zero occurrences of `norm_path` anywhere in the generated doc tree, and `fn.norm_path.html` is gone from the `io` module listing, which now exposes 19 public functions.

### 5. Efficiency — nil, as expected

Visibility is a compile-time annotation with no codegen consequence; the sole call site is unchanged and was already a direct intra-module call. If anything the change is marginally *helpful* — the symbol is no longer externally reachable, so it is a cleaner inlining candidate — but the function is one `to_ascii_lowercase` and the effect is immaterial. Nothing to report.

---

## Plan conformance

All seven implementation steps are satisfied.

| Step | Status |
|---|---|
| 1. `pub fn` → `fn` | done |
| 2. Move trade-off paragraph to `collision_key` | done, verbatim |
| 3. Reduce `norm_path` doc to one line, no "(private)" framing | done — the plan's suggested wording used verbatim |
| 4. Note that input identity uses the case-preserving key; plain backticks | done, no `[…]` links introduced |
| 5. Optional drive-by on the test comment | done, with a deviation worth naming (below) |
| 6. Verification commands | all run; clippy re-run cold because the shared-tree result was cached |
| 7. No CHANGELOG entry | honoured — diff touches `src/io.rs` only |

**Deviation on step 5, benign.** The plan asked the test comment to stop implying `norm_path` is used directly *and* noted it omits that `--passthrough` uses both keys for different questions. The implementation dropped the `--passthrough` mention rather than qualifying it, and the "two keys, two questions" fact landed in `collision_key`'s doc under step 4 instead. That is the better placement — a test comment is not where a reader looks for the key-semantics distinction — so I would not change it back. Noted only because a coverage audit comparing plan text to diff will spot the missing `--passthrough` clause in the test comment and should know it moved rather than vanished.

**Commit message.** It documents the `pub(crate)`-considered-and-rejected reasoning the plan required, and its factual claims check out independently: "561 tests green" matches my run exactly, and I reproduced the E0603 positive control from scratch. Detail sits in the commit message rather than the source, per CLAUDE.md.

**Semver (Low, already covered).** Removing a `pub fn` from `trim-galore` 2.3.0 is technically a breaking API change; `lib.rs` still exposes 17 bare `pub mod` with no `pub use` curation and no `//!` disclaimer. Plan A3 states this accurately rather than denying it, and the durable fix is already filed as a follow-up. No action for this PR.

---

## Recommendations by priority

**Critical** — none.

**High** — none.

**Medium**

1. **M1 — generalise the trade-off's scope** (`io.rs:77-79`). Replace "on opt-in case-sensitive APFS volumes" with a phrasing that includes Linux ext4, which is where this false positive actually lives and where CI exercises it. Rewrite supplied above; folds L3's reflow in at the same time. *Trivial fix.* This is the only item I would ask for before merge.
2. **M2 — file an issue for the 12 pre-existing rustdoc errors**, `src/alignment.rs:11` first (a public doc linking to a private item — the realized form of the hazard this plan reasoned about). Natural companion to the lib.rs-surface follow-up. *Out of scope for this PR.*

**Low**

3. **L4** — "callers" → "production callers" in the test comment (`io.rs:765`). *Trivial fix.*
4. **L3** — reflow the moved paragraph to match its neighbours' ~85-char wrap. *Trivial fix; subsumed by M1.*
5. **L5** — "The fold that `collision_key` is built on" to remove the garden-path. *Trivial fix; optional.*
6. **L2** — leave `path_identity_key` as backticks (consistency with the whole file beats clickability). *No action.*
7. **L1** — `collision_key`'s "these two paths" opener; pre-existing. *No action.*

**Process notes for the caller.**

- **The shared tree is no longer on this commit.** It is on `fix/400-clump-help-text` @ `0676fdf`, whose `src/io.rs` is the pre-#399 version. Any gate that runs in the shared tree from here on is not testing #399. Re-checks must pin the commit (`git archive 315557d` into a scratch dir, private `CARGO_TARGET_DIR`) — which is how every load-bearing result in this review was produced.
- **`fix/400` does not contain #399.** Merge the two independently; neither branch carries the other.
- **A sub-second clippy "success" in the shared tree is a cache hit, not evidence.** The trustworthy form is a cold target dir on a tracked-files-only copy, which also excludes the untracked `examples/fastqc_only.rs` that CI does not have.
- **Don't pipe verification output through `tail`.** It cost me one false negative here: a truncated `cargo test` log looked like a 23-test suite and hid the very test the plan asks to confirm.
