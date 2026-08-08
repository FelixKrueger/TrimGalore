# Plan Coverage Report

**Mode:** B
**Plan(s):** `plans/08082026_hardtrim-conflict/PLAN.md` (v1, 2026-08-08)
**Date:** 2026-08-08
**Verdict:** COMPLETE

## Summary

- Total items: 12
- DONE: 12
- PARTIAL: 0
- MISSING: 0
- DEVIATED: 0

Audited against the uncommitted diff on `fix/386-hardtrim-conflict` (`CHANGELOG.md`,
`src/cli.rs`; 31 insertions, 1 deletion — no other file touched).

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `conflicts_with = "hardtrim3"` on `--hardtrim5` | §2 Step 1 | DONE | `src/cli.rs:371`, exactly the planned line; one side only, as specified |
| 2 | Conflict test: both flags → parse error naming both | §2 Step 2 | DONE | `test_hardtrim5_and_hardtrim3_together_rejected_at_parse`, `src/cli.rs:1228`; asserts both names via `try_parse_from().unwrap_err()` — parse level, as planned |
| 3 | Regression guard: each flag alone parses **and** validates | §2 Step 2 | DONE | `test_each_hardtrim_alone_still_accepted`, `src/cli.rs:1239`; loops `--hardtrim5 20` and `--hardtrim3 15`, calls `validate()` |
| 4 | CHANGELOG entry under `#### Changes` | §2 Step 3 | DONE | `CHANGELOG.md:122-127`, inside `### Unreleased` → `#### Changes` (heading at :112); states accepted-and-ignored → refused, links #386 |
| 5 | Edge: either flag order → same clap error | §3 | DONE | Test covers one order (as the plan scopes it); both orders verified manually — exit 2 each, message names both, subject/other swap symmetrically |
| 6 | Edge: each flag alone unchanged | §3 | DONE | Manual run wrote `BS-seq_10K_R1.20bp_5prime.fq.gz` and `…15bp_3prime.fq.gz`, exit 0 both |
| 7 | Edge: hardtrim + `--paired`/uBAM unchanged | §3 | DONE | No dispatch code touched (diff is CHANGELOG + cli.rs only); pre-existing `hardtrim5_plus_ubam_runs_clean`, `ubam_out_hardtrim5_writes_bam`, `hardtrim_still_accepts_a_mixed_pair`, 8 collision tests all pass |
| 8 | Edge: deliberate Perl v0.6.11 departure documented | §3 | DONE | Rationale stated in the CHANGELOG entry ("Perl v0.6.x behaved the same way") |
| 9 | `cargo fmt --all -- --check` | §4 | DONE | Clean |
| 10 | `cargo clippy --all-targets --release -- -D warnings` | §4 | DONE | Clean |
| 11 | Full `cargo test` | §4 | DONE | **542 tests, 0 failed** (408 lib + 134 integration); matches §6's claimed 542 |
| 12 | Negative control (§6 claim: remove attribute → conflict test fails, alone-tests pass) | §4 / §6 | DONE | Verified by construction, not re-executed (repo left unmodified): the clap attribute is the *sole* mutual-exclusion mechanism — `Cli::validate()` has no hardtrim5-vs-hardtrim3 check (`grep hardtrim src/cli.rs` shows only range/`--clumpify`/`--clump_only`/`--passthrough` guards), so removal makes the parse succeed and `unwrap_err()` panic, while items 3's asserts never touch the attribute |

## Test verification

| Test name | File | Status |
|-----------|------|--------|
| `test_hardtrim5_and_hardtrim3_together_rejected_at_parse` | `src/cli.rs:1228` | PASS |
| `test_each_hardtrim_alone_still_accepted` | `src/cli.rs:1239` | PASS |
| Full suite (542) | crate + `tests/` | PASS (0 failed) |

## Manual reproduction (§4)

| Invocation | Expected | Observed |
|---|---|---|
| `--hardtrim5 20 --hardtrim3 15 R1` | non-zero, both names | exit 2, `error: the argument '--hardtrim5 <HARDTRIM5>' cannot be used with '--hardtrim3 <HARDTRIM3>'` |
| `--hardtrim3 15 --hardtrim5 20 R1` (reverse) | same error | exit 2, names swapped, otherwise identical |
| `--hardtrim5 20 R1` | still trims | exit 0, `BS-seq_10K_R1.20bp_5prime.fq.gz` |
| `--hardtrim3 15 R1` | still trims | exit 0, `BS-seq_10K_R1.15bp_3prime.fq.gz` |

`.fq.gz` rather than the plan's `.fq` is the gzipped-fixture form of the documented
`.<N>bp_5prime.fq(.gz)` pattern, not a deviation.

## Verdict

**Verdict:** COMPLETE

All three §2 steps implemented as written, all four §3 edge-case rows trace to a test or a
stated rationale, and all §4 validation items pass. §6's "no deviations" and its 542-test
figure both check out against the working tree; the recorded skip of plan review was a
user call and is not a coverage gap.
