# Progress: poly-A tails survive behind a poly-G tail (`--poly_a`)

**Last updated:** 2026-10-21

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md` rev 2 (rev 1 at `PLAN.md.v1`) |
| Plan Review | ✅ Complete | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` — both: ship B1, drop B2 |
| Implementation | ✅ Complete | `fix/polya-behind-polyg`, 2 commits, **not pushed** |
| Code Review | ✅ Complete | `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md` — code sound, findings in docs + tests |
| Coverage | ✅ Complete | `COVERAGE.md` — 30/34 DONE, 2 CHANGELOG gaps, both now closed |
| Review fixes | ✅ Applied | all Critical/High/Medium findings, below |

## Final state

| Commit | Scope |
|---|---|
| `a347b4f` | the reorder + 7 `trim_read` cases, `3 ≤ n ≤ 10`, `tests/integration_poly.rs` (5 cases) |
| `d7d2155` | `modes/poly.md`, `migration.md`, `--poly_a` help text, CHANGELOG |

History rebuilt twice, nothing having been pushed: once because the original docs commit asserted
the Perl adapter was `A{10}`, and once to re-scope after the reporter's data refuted the diagnosis.

## Review findings applied

| Finding | Source | Fix |
|---|---|---|
| **HIGH** — `migration.md` inverted the Perl adapter: `A{20}`, not `A{10}` | both code reviewers | rewritten; `PLAN.md` §Context and §Decision corrected, the "3.6×" argument withdrawn |
| Perl `--polyA` *did* apply `--length` (`:1243`) | Review A | bullet corrected |
| `s/\s+/_/g` rewrote *every* read ID, not only trimmed ones | Review A | bullet corrected |
| `poly.md` 3'seq row prescribed the workaround the fix removes | Review B | row now says leave auto-detection on; `--2colour` demoted to fallback |
| Paired R2 test keyed a HashMap on duplicate read IDs — 13,000 → 10,012 silently | Review B | positional comparison, all 13,000 asserted |
| Its anti-vacuity guard was satisfied by records where only one end moved | Review B | counts records where **both** ends moved, checked against the input |
| `--2colour` test had no record count — passed on an empty file | Review A | count asserted |
| 10-line module header, 3 multi-line test comments, one 3-line source comment | both | cut to the one-line rule; the false "only an end-to-end run pins that order" removed |
| No test for `--poly_g` without `--poly_a` | Review B | `poly_g_without_poly_a_still_trims_the_g_run` |
| @marenmsk not credited | Coverage | credited |
| Report-attribution shift undocumented | Coverage | one line in the entry |
| `poly_a` has 4 construction sites, not 1; cited `TempDir` pattern does not exist | Coverage | `PLAN.md` corrected |

Found while applying them: `cargo clippy` on a cold cache rejected an unused `RecordSource` import
that three warm-cache runs had passed. Both code reviewers flagged their clippy as cache-replay —
they were right to.

## Verification after the fixes

| Check | Result |
|---|---|
| Full suite | ✅ 824 tests (516 lib + 308 integration). `PROGRESS.md` previously said 219 — wrong |
| `clippy --all-targets -D warnings`, cold | ✅ clean after the unused-import fix |
| `fmt --check` | ✅ clean |
| Default-flag byte-identity vs `dev` | ✅ SE, PE, `--clock` md5-identical (re-run after the fixes) |
| Docs build | ✅ 32 pages, `#stage-order` anchor present |
| Reported case | ✅ `A40+G15` → 60, `A60+G20` → 40, `A20+G20` → 60 |

## Open

1. **Diagnosis still unconfirmed on the reporter's data** — awaiting her `Reads with poly-A/T
   trimmed` count. The fix is correct regardless; whether it resolves *her* run is not established.
2. Felix pushes the branch and opens the PR.
3. Perl validation matrix, reproducibility and `--version` provenance run in CI only.

## Follow-ups parked in PLAN.md

Terminal-homopolymer skip (both plan reviewers proposed it independently; 0.08–0.82% collateral vs
the withdrawn B2's 25%), fixed-point iteration of the two poly stages, the 5'-prefix case on R2
(`G15 + T40` clips 0), and #245 item B.

---

## Outcome: the diagnosis was refuted (2026-10-21)

@marenmsk's trimming reports settle it against the plan's own Assumption 1, which said a high
poly-A count would invalidate the diagnosis:

| Run | Reads with poly-A/T trimmed |
|---|---|
| `--quality 20 --phred33 --polyA` | 59,257 (25.9%); 1,556,459 bp |
| `--2colour 20 --polyA` | 56,726 (24.8%); 1,461,121 bp |

The poly-G mechanism predicted `--2colour` would raise that count sharply by exposing hidden tails.
It fell, and a quarter of her reads were already being trimmed. **Poly-G ordering is not her
problem.**

Her actual case is the 3'-anchored scan meeting an A-run with sequence behind it. Measured on
`illumina10K_with_polyA.fastq.gz`: the scan leaves 0 reads ending in ≥5 A but 37 containing `A{20}`
somewhere, where `-a "A{20}"` leaves 0 of the latter and 33 of the former. An adapter match removes
everything 3' of the match; the scan cannot.

**Resolved as working-as-designed, documented, not code-changed.** The reasoning, in the order it
carries weight:

1. Aggressive trimming must be opt-in, because a 3'seq user can choose it deliberately while a user
   whose coverage is depleted has no way to notice.
2. The depletion is real and measured — genomic runs of 10+ A's. On `illumina_10K`, `-a "A{10}"`
   removes 213,308 bp and drops 1,358 reads below `--length 20`, against `A{20}`'s 28,064 bp and 138
   reads. Hence the docs prefer `A{20}`.
3. It is invisible to the person paying it: missing coverage, not a QC flag. Worst in bisulfite and
   EM-seq data, where G→A conversion manufactures A-tracts — ~11.6% of 150 bp windows of a converted
   human genome carry a 12-A run (@lyijin and @alice-mcatamney, #180).

Issue history found while answering: #71 (2019) and #74 (2020) record the Thermo Fisher Collibri
provenance of `--polyA`, and in #71 Felix already recommended `-a {A}10` over the flag. #180 (2023)
is where the genomic-depletion concern was first raised and later measured.

Reply posted: <https://github.com/FelixKrueger/TrimGalore/issues/504#issuecomment-6315097983>.
Issue left OPEN for Felix to close.

## The reorder stands on its own merits

`a347b4f` fixes a real bug that `dev` still has — `A20+G20` demonstrably survives there. It is now
scoped as "found while investigating #504", not as its fix, in both the commit message and the
CHANGELOG, which carries two entries: one for the reorder, one for the documentation.

## Open

1. Felix pushes `fix/polya-behind-polyg` and opens the PR; `wiki/init` is also still unpushed.
2. Felix closes #504 when the docs are live.
3. Perl validation matrix, reproducibility and `--version` provenance run in CI only.
