# PROGRESS — `--passthrough` mode (Multiome/10X cell-barcode carrier)

**Plan:** `PLAN.md` (v2)
**Reviews:** `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`
**Started:** 2026-06-06
**Owner:** Felix Krueger

---

## Pipeline status

| Stage                      | Status      | Notes |
|----------------------------|-------------|-------|
| 1. Plan v1 written         | ✅ done     | 2026-06-06 |
| 2. Manual review (v1)      | ✅ done     | Approved with the gitignore recommendation; agreed to advance to agent review |
| 3. Agent review (dual)     | ✅ done     | Dual independent `plan-reviewer` agents, 2026-06-06; reports A + B in this dir |
| 4. Plan v2 revision        | ✅ done     | Folded in 6 consensus + 5 critical-unique + 2 contradiction resolutions; see PLAN.md Revision history |
| 5. Manual re-review (v2)   | ⏳ pending  | User reviews v2; may request further changes before implementation trigger |
| 6. Implementation          | ⏳ pending  | Awaits explicit trigger (`implement` / `/code-implementation`) |
| 7. Verification (dual)     | ⏳ pending  | Dual `code-reviewer` + `plan-manager` post-implementation |

---

## Implementation step tracker

Per the plan v2 "Implementation outline":

- [ ] Step 1 — CLI field + 9-item validation block + case-folded R1/R2 collision check (`src/cli.rs`)
- [ ] Step 2 — Output naming helper + collision pre-flight extension (`src/io.rs`, `src/main.rs`)
- [ ] Step 3 — `PartialEq, Eq` derive + stats fields (`src/report.rs::PairValidationStats`)
- [ ] Step 4 — `read_id_prefix()` helper with `/[123]` stripping (`src/trimmer.rs` or `src/fastq.rs`)
- [ ] Step 5 — Serial path (`src/trimmer.rs::run_paired_end`), 8-arm three-way EOF match
- [ ] Step 6 — Parallel path reader thread (`src/parallel.rs::read_pairs_round_robin`, `PairedWork`)
- [ ] **Step 6a — Reader-thread error propagation hardening** (B-Crit-1; channel payload `Result<PairedBatchResult>`)
- [ ] Step 7 — Parallel path worker + output (`process_paired_batch`, `process_pairs`, `PairedBatchResult`), with explicit gzip/plain fork mirroring + atomic 3-write flush
- [ ] Step 8 — Wire `main.rs` (open readers/writers, sanity-check, dispatch)
- [ ] Step 9 — Report (text + JSON passthrough block)
- [ ] Step 10 — FastQC integration
- [ ] Step 11 — Tests + docs (8 validation tests, fixture, `--help`, CLAUDE.md)

---

## Validation checklist (from plan v2 §Validation)

- [ ] §1 — End-to-end smoke (serial): 50-record fixture; row-by-row lockstep + counts.
- [ ] §2 — **Serial/parallel parity** (load-bearing): **≥10,000 record fixture** (multi-batch), row-by-row R1/R2/passthrough ID-lockstep + decoded record-identity + `PairValidationStats` equality.
- [ ] §3 — Truncated passthrough → loud named-file error.
- [ ] §4 — Shuffled passthrough IDs → hard error at first mismatch.
- [ ] §5 — Each CLI rejection (`--retain_unpaired`, `--clumpify`, `--clock`, `--implicon`, `--hardtrim5/3`, `--demux`, multi-pair, non-paired, file-not-found, **case-folded R1/R2 collision**) covered by a unit test.
- [ ] §6 — **`--dont_gzip` plain-output coverage** *(new in v2)*.
- [ ] §7 — **`--cores 1` dispatcher integration** *(new in v2)*.
- [ ] §8 — **Reader-thread error propagation** *(new in v2; covers B-Crit-1)*.

---

## v2 revision summary

Plan v2 (2026-06-06) folds in dual plan-reviewer agent feedback:

**Consensus items (both reviewers flagged independently):**
- AB1: Parity test fixture ≥10K records (multi-batch).
- AB2: `read_id_prefix()` strips trailing `/[123]` in v1.
- AB3: `--dont_gzip` plain-output test (Validation §6).
- AB4: `--cores 1` dispatcher integration test (Validation §7).
- AB5: Row-by-row lockstep assertion in parity test (Validation §2).
- AB6: Worker plumbing for `process_pairs` spelled out (immutable third slice; gzip/plain fork mirrored).

**Reviewer A critical:**
- A-Crit-1: `#[derive(PartialEq, Eq)]` on `PairValidationStats` — Validation §2 won't compile without it.
- A-Crit-2: Passthrough-longer-than-R1/R2 EOF arm in `read_pairs_round_robin`.
- Empty-input edge case corrected (sanity_check rejects).
- Sync-check cost corrected to ~90 ms / 1M reads.

**Reviewer B critical:**
- B-Crit-1: Reader-thread error propagation hardening (Step 6a; channel-payload change).
- B-Crit-2: Assumption #9 softened — "record-identical" not "byte-identical".
- Output gzip driven from R1 (documented).
- 8-arm three-way EOF match enumerated.
- `passthrough_records_checked` accounting site decided (worker stamps from batch length).

**Contradictions resolved:**
- Stats placement: extend `PairValidationStats` (B's framing — lifecycle matches).
- Severity of case-fold R1/R2 collision: real logic gap, required check (B's framing — consistent with #216).

---

## Open decisions deferred to implementer

(From plan v2 §Questions or ambiguities — Still open, non-blocking)

- FastQC for cell-barcode reads will look poor — documented; `--fastqc-skip-passthrough` opt-out deferred until users complain.
- Auto-cleanup of partial outputs on reader error — deferred; v1 leaves partial files + clear error.

---

## Out of scope (deferred to a follow-up release)

- Multi-pair input with passthrough.
- Clumpy reordering with passthrough.
- Passthrough output suffix customization.
- Specialty-mode combinations.
- Per-pair passthrough flags.
- `--retain_unpaired` rescue semantics (clean v2 path exists — emit on any surviving mate — deferred until requested).
- Auto-cleanup on reader error.

---

## Notes / caveats

- `plans/` is now tracked in Git as of this session (root `.gitignore` updated; `plans/.gitignore` ignores `archived/` only).
- No CI `validation` matrix entry planned (Perl 0.6.11 has no counterpart).
- Reader-thread channel payload change (`PairedBatchResult` → `Result<PairedBatchResult, anyhow::Error>`) in Step 6a is the most subtle existing-code change in the plan; existing R1/R2-only tests should still pass but worth running the full `cargo test` between Step 6 and 7.
- This remains purely additive — no behavioural change in the absence of `--passthrough`.
