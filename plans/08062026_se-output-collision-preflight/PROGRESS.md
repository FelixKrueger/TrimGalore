# Progress: Output-collision pre-flight for the four uncovered dispatch paths (#383)

**Last updated:** 2026-08-07

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` v2; §13 carries deviations D1–D12 |
| Plan Review | ✅ Complete | Two independent Criticals, both reproduced before adoption |
| Impl Plan | ❌ Excluded | No separate `IMPL.md`; `PLAN.md` §5 is the 13-step outline |
| Implementation | ✅ Complete | 530 tests, fmt + clippy clean, no new doc warnings |
| Code Review | ✅ Complete | 4 passes, 3 reports. 0 Critical, 1 unanimous High — **applied** |
| Coverage | ✅ Complete (see caveat) | 3 audits, all COMPLETE. `COVERAGE.md` authoritative |

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

## Coverage-audit index

| File | Items | Verdict |
|---|---|---|
| `COVERAGE.md` — **authoritative** | 65 | COMPLETE (1 DEVIATED, documented) |
| `COVERAGE_round2.md` | 25 | COMPLETE |
| `COVERAGE_round3.md` | 67 | COMPLETE (4 DEVIATED, documented or covered) |
| `COVERAGE_round1.md` | — | INCOMPLETE, 2 items — both closed same day |

**Caveat:** all three COMPLETE audits ran against the tree *before* the post-code-review round
(they report 521 tests; the tree now has 530). That round is documented as D7–D12 in `PLAN.md`
§13 with its own gates and reproductions, but **no independent coverage audit has run against the
final state.** The one deviation `COVERAGE_round3.md` raised — hardtrim's `-o` test shape — is
exactly what D12's new no-`-o` test closes.

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
- `PLAN.md` §11 Open 3 (`--hardtrim5 N --hardtrim3 M` silently ignores the 3′ trim) and Open 4
  (a `modes/hardtrim.md` sentence about CWD output) — both want their own issue.
- A coverage audit against the final tree.

## History

- 2026-08-07: Post-review round applied (D7–D12: `..` folded, one identity definition, secondary
  outputs guarded, CWD hint generalised, V5 fixture made falsifiable); 530 tests
- 2026-08-07: Coverage → ✅ COMPLETE ×3; Code Review → ✅ Complete (4 passes, 0 Critical)
- 2026-08-07: Round-1 coverage gaps closed; D6 recorded
- 2026-08-06: Implementation → ✅ Complete (13 steps; D1–D5, 4 negative controls)
- 2026-08-06: Plan → v2; Plan Review → ✅ Complete (two independent Criticals)
- 2026-08-06: Plan → ✅ Complete
