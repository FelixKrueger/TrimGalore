# Plan — Reject `--hardtrim5` + `--hardtrim3` together (#386)

**Issue:** [#386](https://github.com/FelixKrueger/TrimGalore/issues/386)
**Branch point:** `dev` @ `c4599f5` (or on top of #392's branch if not yet merged)
**Revision:** v1 (2026-08-08)

## 1. Goal

`trim_galore --hardtrim5 20 --hardtrim3 15 x.fastq` exits 0 having written only the 5′
output — the 3′ request is silently dropped, because `main()` dispatches on `--hardtrim5`
first and returns (`main.rs`, the `if let Some(n) = cli.hardtrim5` block; #385's A3
documented the early return). Turn the silent drop into a usage error.

## 2. Implementation

**Step 1 — `src/cli.rs:371`:** add `conflicts_with = "hardtrim3"` to `--hardtrim5`'s
clap attribute. One side suffices; clap reports the conflict symmetrically whichever
order the flags appear in. Precedent: `conflicts_with` is used 7× in `cli.rs` already.

**Step 2 — tests, `src/cli.rs` mod tests:**
- both flags → parse error mentioning both flag names (clap `try_parse_from` errors
  before `validate()`, so test at parse level);
- each flag alone still parses and validates (regression guards — `--hardtrim5 20` and
  `--hardtrim3 15`).

**Step 3 — CHANGELOG**, `#### Changes`: previously accepted-and-ignored, now refused.

## 3. Edge cases

- `--hardtrim5 20 --hardtrim3 15` in either order → same clap error. (Step 2 tests one
  order; clap's conflict handling is symmetric and not ours to re-test exhaustively.)
- Each flag alone, and each with `--paired`/uBAM output → unchanged; no dispatch code is
  touched, so no behaviour beyond the parse-time rejection changes.
- Perl v0.6.11 accepted both and also ran only one (verified in #386's filing) — this is
  a deliberate departure identical in kind to the duplicate-input rejection (#383 A6):
  refusing an invocation whose silent behaviour was a trap.

## 4. Validation

The two Step-2 tests; `cargo fmt --check`; `clippy -D warnings`; full `cargo test`;
manual: the §1 reproduction must exit non-zero with both flag names in the message, and
`--hardtrim5 20 x.fastq` alone must still produce `x.20bp_5prime.fq`. Negative control:
remove the attribute → the conflict test fails, the alone-tests still pass.

## 5. Self-review

One attribute plus tests; no interaction with the #383/#388 pre-flights (rejection is at
parse time, before validate and dispatch). The only judgement call is clap-level vs
`Cli::validate()`-level rejection: clap chosen because the conflict is static (no
cross-field logic needed), the message is free, and 7 existing conflicts set the
precedent. Risk: none identified.

## 6. Implementation notes

Implemented on `fix/386-hardtrim-conflict` off `dev` @ `c4599f5`. All three steps as
specified; no deviations. Plan review was skipped on the user's fast-track call (the
implement trigger was given after reading the plan; option presented explicitly).

Gates: **542 tests** (540 + 2 new), fmt clean, clippy clean. Negative control: attribute
removed → conflict test FAILED, alone-tests passed; reverted → green. Manual: both flags
→ exit 2, message names both; `--hardtrim5` alone still writes `x.20bp_5prime.fq`.

### Post-code-review round

Single reviewer (scaled for diff size, noted in its report) + coverage audit (**COMPLETE
12/12**). Reviewer: 0 Critical/High, 2 Medium, 3 Low. Its verification exceeded the plan:
both flag orders, `=` form, conflict-vs-range precedence, `--help` rendering, and a sweep
proving parse-time is the binary's entire Cli construction surface.

Applied:
- **M1 — moved the rejection from clap to `Cli::validate()`**, reversing the plan's §5
  placement call. The reviewer's decisive point: the *other* hardtrim exclusions
  (`--clump_only`, `--passthrough`, `--clumpify`) are all validate-level with bespoke
  text, so clap's generic message was the family-inconsistent choice — and the remedy
  (sequential invocations feeding output to input) is genuinely non-obvious. A validate
  check is symmetric by construction, so the reviewer's stated trade-off (losing clap's
  free symmetry) does not in fact apply; both orders tested.
- **M2 — the docs `## Compatibility` section** (which the site auto-deploys from `dev`)
  now names the one hard incompatibility the mode has.

Iteration log: the first M1 patch aborted on a stale (pre-rustfmt) test anchor, writing
nothing — and its negative control then no-op'd against the unpatched tree and "passed"
by disabling a check that did not exist. Caught because the control was expected to fail.
Re-applied with disk-read anchors and an asserted control target; the real control fails
the conflict test and the revert restores green.

Final: **542 tests**, fmt clean, clippy clean; both orders exit 2 with the remedy in the
message; each flag alone unaffected. The docs changelog copy lagging #383/#384 was noted
by the reviewer as pre-existing, not this change's regression.
