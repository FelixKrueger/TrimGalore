# Progress: Output-collision pre-flight for the four uncovered dispatch paths (#383)

**Last updated:** 2026-08-08

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` v2; §13 carries deviations D1–D13 |
| Plan Review | ✅ Complete | Two independent Criticals, both reproduced before adoption |
| Impl Plan | ❌ Excluded | No separate `IMPL.md`; `PLAN.md` §5 is the 13-step outline |
| Implementation | ✅ Complete | 530 tests, fmt + clippy clean, no new doc warnings |
| Code Review | ✅ Complete | 4 passes, 3 reports. 0 Critical, 1 unanimous High — **applied** |
| Coverage | ✅ Complete (see caveat) | 3 audits, all COMPLETE. `COVERAGE.md` authoritative for the pre-D7 tree |
| Coverage re-audit | 🚧 Implementing | `REAUDIT.md` (2026-08-08, D7–D13 vs `dev` @ `ce13761`) — **INCOMPLETE, 3 items**: 1 MISSING, 2 PARTIAL. No behavioural gap |

## What shipped

Four dispatch paths gained the #216 output-collision pre-flight — SE trim (FASTQ + uBAM) and
`--hardtrim5/3` (FASTQ + uBAM) — and all collision checks now route through one helper,
`io::preflight_output_collisions`, at 11 call sites. Four defects in the check itself were fixed
along the way, three of them found by review rather than by the plan:

1. the key was a raw path string, so `./x` and `x` named one file under two keys (defect 1);
2. outputs were never compared against the run's own inputs (defect 2);
3. `..` still defeated the key, so #383 reproduced verbatim through its own fix (D7);
4. only *primary* outputs were compared against inputs, so a report or a `--demux` per-barcode
   file could still overwrite a named input while every primary stayed distinct (D9 — the
   unanimous review finding).

`collision_key` is now the crate's single definition of "the same file", used by the pre-flight
and by all four `cli.rs` identity checks (D8).

**Corrected 2026-08-08:** the sentence above describes D8, which **D13 partly reversed** — the
key is split. `path_identity_key` (case-preserved) serves the three input-identity checks;
`collision_key` (case-folded) serves the pre-flight, and the `--passthrough` check uses both
after #389. Both keys share `lexical_normalise`, so D8's behavioural win survives intact.
`CHANGELOG.md:206-211` carries the same stale claim and still needs the fix (`REAUDIT.md` Gap 1).

## Coverage-audit index

| File | Items | Verdict |
|---|---|---|
| `COVERAGE.md` — **authoritative** | 65 | COMPLETE (1 DEVIATED, documented) |
| `COVERAGE_round2.md` | 25 | COMPLETE |
| `COVERAGE_round3.md` | 67 | COMPLETE (4 DEVIATED, documented or covered) |
| `COVERAGE_round1.md` | — | INCOMPLETE, 2 items — both closed same day |
| `REAUDIT.md` — **authoritative for the final tree** | 20 | INCOMPLETE, 3 items (1 MISSING, 2 PARTIAL) |

**Caveat:** all three COMPLETE audits ran against the tree *before* the post-code-review round
(they report 521 tests; the tree now has 530). That round is documented as D7–D12 in `PLAN.md`
§13 with its own gates and reproductions, but **no independent coverage audit has run against the
final state.** The one deviation `COVERAGE_round3.md` raised — hardtrim's `-o` test shape — is
exactly what D12's new no-`-o` test closes.

**Caveat resolved 2026-08-08** by `REAUDIT.md`, which audits D7–D13 plus the four
`CODE_REVIEW_D13.md` findings against `dev` @ `ce13761`. Verdict INCOMPLETE — 3 items, none
behavioural: the CHANGELOG's "one definition of the same file" entry was never revised after
D13 split the key; the D8→D13 regression still has no unit-level `Cli::validate` guard (only
the CI validation job covers it); and `norm_path`'s doc block still omits `path_identity_key`
(tracked in #399). Every requirement that changes what the binary does is implemented and
tested. The three historical `COVERAGE*.md` files were left untouched.

## Code-review index

| File | Author | Independence |
|---|---|---|
| `CODE_REVIEW_A.md` | Reviewer A | independent; its H3 pre-dates the round-1 gap fixes |
| `CODE_REVIEW_B.md` | blind Reviewer B, plus Appendix B from a second B pass | blind — read the plan only |
| `CODE_REVIEW_A_supplementary.md` | warm agent (audited the plan first) | anchored — discount agreement with the plan |

A fourth pass — the one that found `..` and the `--clock`/`--implicon` hint gap — was overwritten
before it could be snapshotted, because several failed agents resumed hours later and rewrote the
report paths they had been given. Its findings survive in `PLAN.md` §13 (D7, D10) and two were
reproduced against the binary before adoption; its report file is lost.

## Outstanding

- Symlink aliasing — the sole remaining residual of A8, disclosed in the CHANGELOG with its cost.
- `--fastqc_args "-o DIR"` — A2's named exception; costs a regenerable QC artifact, not reads.
- ~~`PLAN.md` §11 Open 3 (`--hardtrim5 N --hardtrim3 M` silently ignores the 3′ trim) and Open 4
  (a `modes/hardtrim.md` sentence about CWD output) — both want their own issue.~~ **Both closed
  2026-08-08:** Open 3 → #386 / PR #393 (`03d793f`); Open 4 → #387 / PR #394 (`7ea2741`), which
  covered all three CWD modes rather than hardtrim alone.
- ~~A coverage audit against the final tree.~~ **Done** — `REAUDIT.md`, 2026-08-08.

Raised by the re-audit (see `REAUDIT.md` Gaps, in priority order):

- **Gap 2, MISSING, untracked.** The D8→D13 regression has no unit-level guard. A four-line
  `Cli::parse_from([… Sample_R1, Sample_R2, SAMPLE_R1, SAMPLE_R2]).validate().is_ok()` test in
  `cli.rs` would close it. Currently only `.github/workflows/ci.yml:625-630` catches the
  overreach, so a repeat would pass a full local `cargo test`.
- **Gap 1, PARTIAL.** `CHANGELOG.md:206-211` still claims all four identity checks "use the same
  key as the pre-flight" — untrue after D13, and the split is disclosed nowhere in the CHANGELOG.
- **Gap 3, PARTIAL, tracked in #399.** `norm_path`'s doc block names `collision_key` but not
  `path_identity_key`, so the split is invisible from the function a reader hits first.

Also from the re-audit: `CODE_REVIEW_D13.md` findings 1 and 2 — the two with reachable harm —
were both closed by later PRs (#388 / PR #392 implemented the recommended paired-report fix;
#389 / PR #395 went further than the finding asked on the `--passthrough` key semantics).

## History

- 2026-08-08: Coverage re-audit → 🚧 Implementing (`REAUDIT.md`: INCOMPLETE, 3 items — 1 MISSING,
  2 PARTIAL, none behavioural). Closes the "no audit against the final state" caveat. §11 Opens 3
  and 4 recorded closed (#386/#393, #387/#394); `CODE_REVIEW_D13.md` findings 1 and 2 recorded
  closed (#388/#392, #389/#395). Targeted suites green on `dev` @ `ce13761`: 45 integration
  collision tests, 54 `io::tests`
- 2026-08-07: Post-review round applied (D7–D12: `..` folded, one identity definition, secondary
  outputs guarded, CWD hint generalised, V5 fixture made falsifiable); 530 tests
- 2026-08-07: Coverage → ✅ COMPLETE ×3; Code Review → ✅ Complete (4 passes, 0 Critical)
- 2026-08-07: Round-1 coverage gaps closed; D6 recorded
- 2026-08-06: Implementation → ✅ Complete (13 steps; D1–D5, 4 negative controls)
- 2026-08-06: Plan → v2; Plan Review → ✅ Complete (two independent Criticals)
- 2026-08-06: Plan → ✅ Complete
