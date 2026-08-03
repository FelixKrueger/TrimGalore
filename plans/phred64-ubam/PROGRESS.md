# PROGRESS — `--phred64` uBAM fixes (issue #358)

**Plan:** [PLAN.md](PLAN.md) — **v2**, revised after dual plan review
**Issue:** [#358](https://github.com/FelixKrueger/TrimGalore/issues/358)
**Base:** `dev` @ `2e3da11`
**Branch (proposed):** `fix/phred64-ubam`

## Pipeline status

| Step | Status | Notes |
|---|---|---|
| Investigation | ✅ done | Both bugs reproduced; scope widened 3 → 7 production write sites |
| Plan v1 | ✅ done | — |
| Manual review | ✅ done | Felix approved proceeding to agent review |
| Agent review (dual `plan-reviewer`) | ✅ done | **Both APPROVE WITH REVISIONS.** A: 1 critical / 4 important. B: 3 critical / 5 important |
| Plan v2 (revisions folded) | ✅ done | All 3 critical + all 8 important + worthwhile optional |
| Implementation | ⛔ **blocked** | Requires the explicit trigger ("implement" / `/code-implementation`) |
| Verify (dual `code-reviewer` + `plan-manager`) | ⛔ blocked | — |

## What dual review changed — two material errors in v1

**1. Implementation ordering was backwards (both reviewers, independently).**
v1 said the Bug 2 guard was "independent of steps 1–7." It is a **correctness prerequisite**. `--clump_only` never interprets quality and `BamReader` normalises to Phred+33, so today's hardcoded `saturating_sub(33)` is *correct* on BAM-input paths even with `--phred64`. Re-verified:

```
--clump_only --phred64 --output-format ubam ubam_test.bam
  input QUAL md5  = 31878eb2…
  output QUAL md5 = 31878eb2…   identical -> correct today
```

Applying the writer fix without the guard subtracts 64 from Phred+33 bytes → saturates → **all-Q0 BAM, every read retained, zero exit status, no warning.** Worse than either original bug: Bug 2 fails loudly, this fails silently. v2 inverts the order and records that one-PR is a correctness constraint, not a shape preference.

**2. The fixture-encoding claim was false (both reviewers).**
v1 asserted "every fixture in `test_files/` is Phred+33." `truncated.fq.gz` is **Phred+64** — min byte 66 / max 99 → Q2..Q35 at offset 64 (textbook Illumina 1.5) vs Q33..Q66 at offset 33; `HWUSI-EAS611…#0/1` GA/GAII read names; trailing `B` runs are the Illumina 1.5 read-segment QC indicator. v1 also miscounted the ambiguous set (named six files, only three have min 73). Heuristic withdrawn; the conclusion survives on the grep check alone (zero `phred64` hits in `tests/` or CI).

## Other revisions folded in

| From | Change |
|---|---|
| A + B | Comment being rewritten has a **second** falsehood — "FastqReader rejects sub-33 input" is untrue; `saturating_sub` is load-bearing, not defensive |
| A | Validation promoted from manual shell steps to **integration tests** (7 of 13) — the bug class is faulty *wiring*, which a unit test can't pin |
| A | Docs retargeted to `guide/quality.md`; `outputs.md`'s list is headed "Rejected at CLI-validate time", which the guard deliberately is not |
| A | Folded in `bam.rs:920` `+` → `saturating_add` — real `u8` overflow on mixed-`0xFF` QUAL, two lines from a touched site |
| B | Field renamed **`input_phred_offset`** — `bam.rs:60` already defines `PHRED_OFFSET = 33` with *different* semantics (read-side, fixed) |
| B | `debug_assert_eq!` is compiled out under `cargo test --release` (what CI runs) — no longer cited as coverage |
| B | §3.1's "the defect is read-side" justification is **false** for hardtrim / clump-only, where the flag is inert and invocations work today. Replaced with the input-format argument; contradiction with O-2 resolved |
| B | Guard preserves `clump_only`'s documented lossless-qual invariant (`clump_only.rs:721-726`) — now recorded |
| B | FastQC-on-uBAM is an unlisted beneficiary (pre-fix plots centre on Q71) |
| B | Archived bad output **cannot** be repaired by re-running — needs −31, not −64. Regenerate from original FASTQ |
| A + B | Counts/refs corrected: **14** call sites (not 13), **409** tests (not 407), `cli.rs:970-972`, `main.rs:237-260`; `specialty.rs` needs no clippy allow; `#### Fixes` already exists |
| A + B | **Assumption 8 resolved** — all four sites take one added parameter, zero intermediate threading |

## Behaviour regression to disclose

Two invocations work correctly today and become hard errors. Needs a CHANGELOG *behaviour-change* line, not just a fix line:

- `--clump_only --phred64 <bam> [--output-format ubam]`
- `--hardtrim5/3 N --phred64 <bam>`

Rejection is still correct (the flag is meaningless for BAM input in every mode), but it is the only place the change removes working behaviour.

## Verification discipline note

Every reviewer factual claim was independently re-verified before adoption. One corrected Reviewer B in passing: B's "`truncated.fq.gz` has zero references" is right *for the fixture* — a broader `grep truncated` matches unrelated error-message text in `src/`.

## Next action

Felix reviews PLAN v2. Implementation remains blocked on an explicit trigger.
