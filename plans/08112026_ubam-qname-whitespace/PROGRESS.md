# PROGRESS — #415 uBAM QNAME whitespace drops preserved aux tags

**Issue:** [#415](https://github.com/FelixKrueger/TrimGalore/issues/415)
**Status:** 📝 Planning — **v1 written, awaiting Felix's manual review**
**Base:** `dev` @ `1fc00f6`

## Pipeline

| Step | Artifact | State |
|---|---|---|
| Reproduction | inline (probe deleted) | ✅ defect confirmed on both output paths at `1fc00f6` |
| Plan | `PLAN.md` | ✅ v1, **awaiting Felix's review** |
| Plan review (dual) | `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` | ⏸️ blocked on Felix |
| Implementation | code | ⏸️ blocked on the exact trigger |
| Coverage audit | `COVERAGE.md` | ⏸️ |
| Code review (dual) | `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md` | ⏸️ |

## Decisions taken (Felix, 2026-08-11, before drafting)

| Question | Answer |
|---|---|
| How hard should the refusal be? | Hard error, always — matches #409's treatment, makes #408's A5 enforced |
| Which characters? | All ASCII whitespace — space, tab, `\n`, `\r`, `\v`, `\f`. Not the full SAM regex, which would over-reject UTF-8 and `@` |

## Decided in the plan, flagged for review

- The check is **read-side**, so it also refuses FASTQ-output runs where nothing is
  currently lost (PLAN §2.3 / O1). The one decision most likely to be overturned.
- The secondary defect (#406's notice saying "FASTQ" on a uBAM run) becomes **unreachable**
  once whitespace QNAMEs are refused, so it is resolved by deleting the possibility and
  pinning it with a test, not by rewording (§2.6).

## Findings from reproduction

1. **uBAM output:** exit 0, QNAME truncated, `--preserve-tags` tag silently gone.
2. **FASTQ output:** tag text survives in the header, but the header is unreadable to any
   whitespace-splitting consumer — so the paths differ and the issue only describes the first.
3. **`fastq::id_prefix` truncates the paired sync key** at the first whitespace — a third
   symptom of the same root cause, not in the issue.
4. **A newline in QNAME would break the FASTQ record structure** — why the check covers all
   ASCII whitespace.
5. **`noodles` refuses to write a whitespace QNAME** (`InvalidInput`), so the fixture must be
   built with samtools and committed — and our own writer can never emit this shape, making
   the defect strictly read-side.
