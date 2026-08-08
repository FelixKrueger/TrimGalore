# Code Review A — #389 --passthrough alias-check wording (commit 1f43363)

**Reviewer:** A (independent)
**Base:** dev @ 7ea2741 · **Head:** 1f43363 (`fix/389-passthrough-wording`)
**Files:** `src/cli.rs`, `src/io.rs`, `CHANGELOG.md` (diff stat matches plan scope; no stray changes)
**Verification run:** `cargo test passthrough` → 26/26 green (incl. the two updated pins and the new case-variant test); `cargo fmt --all -- --check` clean; `cargo clippy --all-targets -- -D warnings` clean (dev profile; impl notes report release clean); plan Validation greps 3/4 re-run and matched the documented deviations exactly.

## Summary

The implementation matches the approved plan (r2) including its documented deviations. The two-tier branch provably preserves the accept/reject surface, each message is truthful for its subcase, the `./`-spelling lands in the identity branch, and the new portable case-variant test is sound on both ext4 and APFS (verified by reading the validate() check order, not just the comment). The backslash continuations are correct at byte level and confirmed by rendering the real message through the built binary. Findings are three Lows, all polish: a matched-input selection corner on a pathological (and independently rejected) input, two 3-line test comments vs the ≤2-line house convention, and one word in the `norm_path` rustdoc. **Verdict: APPROVE — ship as-is; Lows are optional follow-ups.**

## Issues by area

### 1. Logic

**L1 — Accept/reject surface preserved (verified, OK).**
Old code bailed on `collision_key(pt) == collision_key(input[0]) || … input[1]`; new code bails iff `find` on the same predicate over the same two inputs returns `Some`, and both arms of the inner branch `bail!`. Key fact making the branch safe: `collision_key(p)` is exactly `path_identity_key(p).to_ascii_lowercase()` (both are string views of the same `lexical_normalise(p)` output — `io.rs:75-83`), so identity-equal ⇒ collision-equal and the discriminator can never widen or narrow the found set. Position in the chain is unchanged (after 1.viii, before the input-existence loop at `cli.rs:994-997`). The retained redundant `len == 2` guard is kept and explained per plan.

**L2 — Message truthfulness per subcase (verified, OK).**
Because collision-equal ∧ identity-differ ⟺ "differ only by ASCII letter case" (from the L1 identity above), the case-variant message's two claims — case-insensitive comparison, and "rename one so the paths differ by more than letter case" — are precise, not approximate. The identity branch fires exactly on lexical identity (absolutised, `.`/`..` folded, case preserved), which is the same file on every filesystem, so "{pt} is input {matched}" is truthful everywhere. `./R1.fq` vs `R1.fq` takes the identity branch (`lexical_normalise` folds `CurDir`, `io.rs:65`) and renders "./R1.fq is input R1.fq" — truthful. Symlink spellings never reach 1.ix (different names ⇒ different collision keys) — unchanged behaviour, correctly out of scope.

**L3 — LOW: matched-input selection ignores a byte-identical second input (trivial fix).**
`.find()` takes the *first* collision-key match and the identity discriminator is evaluated only against it. If the passthrough is byte-identical to input 2 but a case-variant of input 1, the case-variant branch fires and names input 1. Reproduced against the built binary:

```
$ trim_galore --paired --passthrough $D/r1.fq  $D/R1.fq  $D/r1.fq
Error: --passthrough matches input $D/R1.fq case-insensitively (for APFS/NTFS safety): $D/r1.fq. … if they are genuinely two files, rename one …
exit: 1
```

Here pt IS input 2, so the identity message ("must be a third file … is input $D/r1.fq") is the truthful one. Reachable because `validate_paired_input`'s within-pair check (`cli.rs:537`) is identity-based (case-sensitive), so a case-variant R1/r1 pair passes it. Impact bounded on three sides: the run is still rejected (exit 1, surface unchanged); the message's own "if they are genuinely two files" hedge partially covers it; and per the #388 CHANGELOG entry, a case-variant R1/R2 pair is refused by the report-collision pre-flight on every filesystem anyway — so this only mislabels the subcase of a doubly-invalid command line. **Trivial fix** (found set unchanged since identity-equal ⇒ collision-equal, so the surface is provably preserved):

```rust
let pt_id = crate::io::path_identity_key(pt);
if let Some(matched) = self
    .input
    .iter()
    .find(|p| crate::io::path_identity_key(p) == pt_id)
    .or_else(|| self.input.iter().find(|p| crate::io::collision_key(p) == pt_key))
```

(the existing identity/case branch below then works unchanged).

**L4 — New portable test is sound on BOTH ext4 and APFS (verified, OK).**
Confirmed by reading the actual check order, not the test comment: (1) `validate_paired_input` at `cli.rs:621` runs first but its within-pair duplicate check is `path_identity_key`-based, so `<tmp>/R1.fq` vs the R2 fixture passes; (2) 1.i–1.vii don't apply; (3) 1.viii (`check_restartable_input`, `cli.rs:487-501`) is a pure `fs::metadata` stat of the *passthrough* file only, which the test writes; (4) 1.ix bails before the input-existence loop (`cli.rs:994-997`) ever sees the unwritten `<tmp>/R1.fq`. The comparison is lexical end-to-end, so whether the OS treats `R1.fq` as an existing alias of the written `r1.fq` (APFS) or as a nonexistent path (ext4) is irrelevant. Test observed green in this review's run (on APFS; impl notes same).

**L5 — Informational (no action): test temp-dir hygiene.**
If `validate()` unexpectedly returned `Ok`, `unwrap_err()` panics before `remove_dir_all`, leaking `<temp>/tg_pt_case_<pid>`. Pid-keyed naming is safe for the in-process parallel test runner (one test uses the prefix; threads share the pid). Acceptable for a unit test; noted only.

### 2. Efficiency

Nil, as the plan states. Error path only: the `find` performs the same ≤2 `collision_key` evaluations the old OR did; the found case adds two `path_identity_key` evaluations and one branch before the process exits. No hot-path change; no allocation change worth naming.

### 3. Errors

**E1 — Backslash continuations correct (verified, OK).**
All four continuation lines (`cli.rs:942, 949, 950, 951`) carry the space *before* the `\`; checked at byte level (`awk` over the raw lines) and end-to-end by rendering the message through the built binary — "APFS/NTFS safety", "the same file", "genuinely two files" all correctly spaced. Neither pinned prefix ("--passthrough must be a third file", "case-insensitively") straddles a continuation, per plan Validation 4's rider.

**E2 — Stale-pin sweep clean (verified, OK).**
Plan Validation 3 grep over `src/ tests/ docs/ .github/ README.md CHANGELOG.md` returns exactly the three deliberate hits the implementation notes document: the two negative assertions (`cli.rs:1890, 1937`) and the CHANGELOG sentence describing the *old* text. `tests/integration_passthrough.rs` has no 1.ix pin (confirmed empty grep). Validation 4 counts match the documented deviation ("must be a third file" ×2 via the shared helper; "case-insensitively (for APFS/NTFS" ×1, bare adverb pinned in the new test).

**E3 — No unhandled conditions introduced.** Both arms `bail!`; no new fallible operations on the error path (`display()` is infallible; `path_identity_key`/`collision_key` fall back to the raw path on `std::path::absolute` failure, pre-existing behaviour).

### 4. Structure

**S1 — LOW: two test comments exceed the ≤2-line house convention (trivial fix).**
`cli.rs:1897-1899` (3 lines; the plan itself specified two) and `cli.rs:1914-1916` (3 lines). Both compress cleanly: drop "the case-variant branch has its own lexical test below" from the first (the adjacent test is self-evident), and the "(#389)" tail plus "input existence is validated later" clause from the second can merge into two lines. The rewritten preamble (`cli.rs:885-887`, 3 lines) is a pre-existing block *shortened* from five — fine.

**S2 — LOW (nit, optional): `norm_path` rustdoc says collision_key "is what every collision check hashes" (`io.rs:39-40`).**
1.ix compares via `find`, it doesn't hash. "uses" or "keys by" would be exact. One-word trivial fix; harmless as-is.

**S3 — OK: naming and placement.** `pt_key` (was `pt_norm`) now matches the `*_key` family; `matched` naming follows the PR #219 both-paths precedent; helper name `assert_passthrough_identity_rejection` is accurately scoped to the identity branch and centralizes the three pins without over-abstracting (the case-variant test correctly does not reuse it — different prefix). CHANGELOG entry sits under Unreleased → `#### Changes` (not "Bug fixes"), above #383, modeled on the message-only entry style, and its behaviour claim ("Which runs are accepted or rejected is unchanged") is accurate per L1.

**S4 — OK: companion sites.** `io.rs` rustdoc no longer claims direct use by `Cli::validate()` and drops the "aliasing" identity verb; the retained "pragmatic trade-off" paragraph stays accurate — and its case-sensitive-APFS false-positive reader is now exactly the audience the new "if they are genuinely two files" clause serves. The preamble's em-dash is comment-only (the CI em-dash grep targets `--version` output; the new user-facing strings contain none — no interference).

## Recommendations

| # | Priority | Item | Fix |
|---|----------|------|-----|
| 1 | Low | L3 — prefer an identity-equal input when selecting `matched`, so a passthrough byte-identical to input 2 but case-matching input 1 gets the identity message naming the right file | **Trivial fix** — `or_else` chain shown in L3; surface provably unchanged |
| 2 | Low | S1 — trim the two 3-line test comments (`cli.rs:1897-1899`, `1914-1916`) to ≤2 lines per house convention | **Trivial fix** — deletions only |
| 3 | Low | S2 — `io.rs:40` "hashes" → "uses" (1.ix compares, it doesn't hash) | **Trivial fix** — one word, optional |

No Critical, High, or Medium findings.

## Verdict

**APPROVE.** The change does exactly what the plan promised — same accept/reject surface (provable from `collision_key = lowercase(identity_key)`), truthful subcase-specific messages, sound portable test, clean sweeps — and the three Low findings are polish that can land now as trivial fixes or ride a follow-up without risk.
