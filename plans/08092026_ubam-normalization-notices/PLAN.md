# Plan: make the uBAM path's normalization messages true (#406)

**Issue:** [#406](https://github.com/FelixKrueger/TrimGalore/issues/406) — filed from #400's code review.
**Scope (Felix, 2026-08-09):** fix the IUPAC direction; add a once-per-run notice for dropped header descriptions; uppercasing stays silent. **r2 additions decided 2026-08-09:** echo the actual dropped text; keep once-per-process; fold in the third documentation surface.
**Base:** `dev` @ `90dda11`.

## Revision history

- **r2 (2026-08-09):** dual plan review (PLAN_REVIEW_A.md, PLAN_REVIEW_B.md). Both reviewers reproduced their findings on the built binary rather than reasoning. **Three of r1's claims were wrong** and are corrected below: the "nearly unreachable" shared-`OnceLock` argument (reachable, both reviewers demonstrated it), the "read-side text is already correct" claim (false on uBAM→uBAM), and one clause of the drafted notice (false in exactly the case that fires it). Scope grows by the third documentation surface and shrinks by one false promise. A separate, more serious bug found during review is filed as **#408** and deliberately not fixed here.
- r1 (2026-08-09): initial plan. Superseded.

## Goal

Three messages on the uBAM path are either misdirected or missing. After this change every statement the uBAM path makes about its own normalizations is true in every direction it can run.

## Context

### The IUPAC message is shared by two opposite directions, and is wrong in **both**

| Site | Function | Direction | r1's verdict | Truth |
|---|---|---|---|---|
| `bam.rs:823` | `validate_and_normalize_seq_for_write` | FASTQ → uBAM | reversed | **reversed** — says input uBAM, output FASTQ; both false |
| `bam.rs:944` | `bam_record_to_fastq` | uBAM → FASTQ **or** uBAM → uBAM | "already correct" | **correct only for uBAM→FASTQ.** On uBAM→uBAM it says "coerced to N **for FASTQ output**" when the user's output is a BAM |

Both reviewers reproduced the uBAM→uBAM case. Fixing only the write side would have shipped #406's own defect on the other side of the same sentence. **Remedy: make the text output-neutral — name *where* the coercion happens, not what it is for.** "coerced to N on read" is true for both read-side runs; the write-side message names its own direction.

### The shared `OnceLock` suppression is reachable — r1's argument was wrong

r1 claimed cross-direction suppression was "nearly unreachable" because uBAM→uBAM coerces on read first. That half is true, **but it is not the only way to reach both sides in one process.** A **mixed FASTQ + uBAM single-end multi-input run** reaches both, and such runs are legal by design:

- `any_bam` (`main.rs:240-242`) is an `.any()` over per-input formats — mixed input is tolerated by construction.
- The design *depends* on it: `main.rs:296` offers "convert **at least one** input to uBAM" as the remedy for `--preserve-tags` with all-FASTQ inputs.
- `reject_bam_format_mismatch_in_pair` is gated on `cli.paired` — it does not apply to SE.

Reviewer A reproduced both argument orders. FASTQ-first is the damaging one: the user is told the wrong direction **and** never told about the genuine uBAM coercion, because the first (wrong) message consumed the guard. So splitting into two independent `OnceLock`s **repairs a reachable suppression bug** — it is not the tidiness item r1 described, and it belongs in the CHANGELOG.

### The description drop

- `parse_name_and_data` (`bam.rs:707-719`) splits `record.id` at the first whitespace: TAB → aux-tag tail (parsed), SPACE → discarded. **Write-side only**, one production caller (`BamWriter::write_record`, `bam.rs:589`), which serves **every** uBAM output path (trim, clump, hardtrim).
- The drop is **structurally required** — BAM QNAME forbids whitespace. The defect is that it is silent.
- **r1's frequency claim was overstated.** Reviewer A surveyed all 20 FASTQ fixtures: **9 have no whitespace at all** and would never fire the notice, including the classic non-SRA Illumina `#45/1` style. So this is common, not near-universal — which is an argument *for* lazy per-record detection over eager startup emission, since roughly half of real headers have nothing to drop.
- **There is no flag that preserves a description.** Only TAB-separated `TAG:TYPE:VALUE` tails survive.
- **A sibling message already exists on this exact path**: `NOTE: input is FASTQ; output uBAM records will have no aux fields.` (`main.rs:2021`). It establishes both the register and the precedent that an informational caveat here is acceptable. Note it fires **per input file** while this notice is specified per process (Felix's call) — accepted divergence, recorded so it is not mistaken for an oversight.

### Two edge cases r1's Behavior 6 missed (both reviewers, independently)

1. **Trailing space → false-positive notice.** `FastqReader` strips only `\n`/`\r` (`fastq.rs:446`), so `@name ` reaches the parser with a trailing space, takes the SPACE branch, and yields `rest == ""`. Keying naively on "SPACE branch taken" — literally what r1 specified — would announce a loss that did not occur. **Gate on a non-empty description.**
2. **SPACE before TAB → aux tags are silently lost too, so r1's own clause was false.** For `@readX desc<TAB>CB:Z:ATCG` the first-whitespace split discards the description *and* the tag tail; the resulting BAM record carries no `CB`. r1's drafted sentence "Tab-separated TAG:TYPE:VALUE fields are preserved" is therefore **false on a record that fires the notice** — the same defect class #406 exists to remove. **Drop the promise.** (Not reachable from uBAM input, since BAM QNAME cannot contain a space; needs a hand-built or third-party FASTQ.)

### Out of scope, filed: #408

`--rename --output-format ubam` silently discards the `:clip5:`/`:clip3:` annotation whenever the header has a space, because `append_to_id` (`fastq.rs:173-181`) splices before the first TAB but has no SPACE equivalent, and `--rename` is not in §3.4a's uBAM rejections. Independently reproduced (`@withspace 1:N:0:ACGTAC` → FASTQ keeps `:clip5:ACG`, uBAM QNAME is bare `withspace`). Filed as **#408**; fixing it is a behaviour decision, not this plan's business. **Consequence inside this plan:** the notice must not imply that only *user* metadata is at stake, or a `--rename` user who does not care about instrument IDs will dismiss it and never learn their annotation is gone. Echoing the actual dropped text (Felix's decision) handles this by construction — the `:clip5:` suffix appears in the echoed string.

### The third documentation surface (folded in per Felix)

The description drop / uppercasing / IUPAC coercion are documented **only** on the `--clump_only` surface (`cli.rs:181-193`, `docs/.../modes/clump-only.md`) after #400. Meanwhile `docs/.../guide/outputs.md`'s uBAM section states *"The output records themselves round-trip losslessly"* — false for FASTQ input on all three counts — and `--output-format`'s own help (`cli.rs:254-258`) says nothing. Both get one sentence.

## Behavior

1. **IUPAC write side (FASTQ→uBAM):** names its own direction — input FASTQ, coerced in the uBAM output.
2. **IUPAC read side (uBAM→FASTQ *and* uBAM→uBAM):** output-neutral — "coerced to N on read".
3. **Each direction fires at most once per process, independently**; neither can suppress the other.
4. **New description notice**, `NOTE:`-prefixed, once per process, fired only when a **non-empty** description is discarded, and **echoing the actual dropped text** (truncated to ~60 chars). Makes no promise about aux tags.
5. **No change to what is written** — bytes, names, tags, files, exit codes all identical. Only stderr changes.
6. **Prefix is `NOTE:`, decided from the codebase, not taste** (Reviewer B): `WARNING:` ×10, `NOTE:` ×5, `Note:` ×0 in `src/*.rs`. `WARNING:` here is reserved for "you may have made a mistake"; `NOTE:` for format-mandated, non-actionable facts — which is what this is, and what the three sibling uBAM-path messages use.
7. **Edge cases:** no whitespace → no notice. TAB tail → no notice (tags parsed). Trailing/empty description → **no notice** (Behavior 4's gate). TAB tail containing spaces in a tag *value* → no notice, tags intact (verified by Reviewer A: `CB:Z:has space here` round-trips). SPACE before TAB → notice fires and the tags are lost, which the wording no longer contradicts.

## Implementation outline

1. `bam.rs` — replace `emit_iupac_warning_once` with two independently-guarded emitters (separate `OnceLock`s). Update `:823` (write) and `:944` (read).
2. Write-side text: `WARNING: input FASTQ contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V); coerced to N in the uBAM output. This is legitimate in PacBio HiFi, ONT Dorado, and 10x cellranger data.`
3. Read-side text: as today but output-neutral — `…coerced to N on read.` (keeps the PacBio/ONT/10x reassurance).
4. `parse_name_and_data` SPACE branch — when `!rest.is_empty()`, call `emit_description_dropped_once(rest)`. Message states that BAM QNAME forbids whitespace so text after the first space is not carried into the output, and echoes the dropped text truncated to ~60 chars (ellipsis if truncated). **No aux-tag promise.**
   - **Emission stays in the parser** rather than returning a flag to `write_record` (Reviewer B's option 10). Rationale: the function already returns `Result` and performs validation, so it is not pure today; threading a value through would touch ~16 test call sites for no testable gain, because per B's own item 11 once-ness must not be asserted in unit tests anyway (the `OnceLock` is process-global and `validate_seq_coerces_iupac_to_n` already consumes the write-side guard in that binary). The notice is therefore tested at the **integration** level, which is the right altitude. Recorded as a deliberate deviation from B's suggestion.
5. Doc comments — `bam.rs:796-800` and `:865` describe the coercion without naming direction; one clause each.
6. **Third surface** — one sentence in `docs/.../guide/outputs.md`'s uBAM section replacing/qualifying "round-trip losslessly", and one in `--output-format`'s help (`cli.rs:254-258`) pointing at the per-record normalizations.
7. **Tests** — binary-driven, in the integration suite (not unit): (a) FASTQ with IUPAC → uBAM shows the write-side direction and not the read-side text; (b) uBAM with IUPAC → uBAM shows "on read" and not "for FASTQ output"; (c) a mixed FASTQ+BAM SE run shows **both** messages (this is validation 3, now known constructible); (d) a space-bearing header → notice present, echoing the dropped text; (e) `@name ` trailing-space and a TAB-tail header → **no** notice. Never assert once-ness.
8. `CHANGELOG.md` — one bullet under `### Unreleased` → `#### Bug fixes`, covering **both** the reversed text and the reachable suppression bug (r1 planned to mention only the text). Avoid the `#### Fixes` duplication hazard (#401).
9. `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings`, `cargo test`.

## Efficiency / Integration

Both notices are `OnceLock`-guarded: after the first hit the cost is one **Acquire** load per coercion/split (r1 said "relaxed" — corrected per Reviewer B). `parse_name_and_data` already computes the split, so detection adds nothing; the echo allocates one truncated `String` **once per process**, inside the guard. No new I/O ordering, no behavioural coupling.

## Assumptions

- **A1:** nothing pins any of these messages — verified across `src/`, `tests/`, `.github/`, `docs/` by both reviewers.
- **A2:** the description drop is write-side only, one production caller, serving all three uBAM output paths — verified.
- **A3 (corrected):** the read-side text is **not** already accurate; it is wrong on uBAM→uBAM.
- **A4 (corrected):** the notice fires **commonly, not near-universally** — 9 of 20 fixtures have whitespace-free headers. Lazy per-record detection is therefore the right mechanism; eager startup emission would over-disclose.
- **A5:** adding a stderr line cannot break an existing assertion — both reviewers verified no test asserts on stderr exhaustively or on line counts.
- **A6 (new):** echoing user data into stderr is acceptable here — a read name is not sensitive, and it is truncated. Flagging so the choice is explicit rather than incidental.

## Validation

| # | What | How | Expected |
|---|------|-----|----------|
| 1 | Write-side direction | IUPAC FASTQ → `--output-format ubam` | names input FASTQ + uBAM output; no "input uBAM"/"for FASTQ output" |
| 2 | Read-side output-neutral | IUPAC uBAM → uBAM, **and** → FASTQ | "on read" in both; the phrase "for FASTQ output" appears nowhere |
| 3 | Both directions, no suppression | **mixed FASTQ + uBAM SE run, both IUPAC-bearing, distinct stems** (constructible — confirmed by both reviewers; run it in **both argument orders**) | both messages present in both orders |
| 4 | Notice fires and echoes | `BS-seq_10K_R1.fastq.gz` → uBAM | notice once, containing the real dropped tail |
| 5 | Notice does **not** over-fire | `@name ` trailing space; a TAB-tail header; `BS-seq_10K_I1` (whitespace-free) | no notice in any; aux tags still parsed |
| 6 | Output bytes unchanged | the **existing** harness — `assert_ubam_eq` against `test_files/ubam_out_{se,pe}_REFERENCE.bam`, already run by `cargo test` (both reviewers; r1 proposed a manual diff) | identical |
| 7 | Expected-fail control | run validations 1, 2 and 3 against the **unpatched** build first | 1 shows reversed text; 2 shows "for FASTQ output"; 3 shows one message suppressed |
| 8 | No collateral | full `cargo test` + fmt + clippy `-D warnings` | green |

## Questions or ambiguities

- **[Resolved r2 — Felix]** Echo the actual dropped text; once per process; fold in the third surface.
- **[Resolved r2 — evidence]** Prefix is `NOTE:` (codebase frequency + sibling register).
- **[Open, minor]** Truncation length (~60 chars) is a guess; a reviewer may prefer a different budget or no truncation.

## Self-Review (r2)

- **What r1 got wrong, owned:** three claims, all load-bearing, all reproduced-wrong by both reviewers independently — the unreachability argument, the read-side accuracy claim, and a notice clause that was false precisely where it fires. The last is the second consecutive round in which my own corrective copy contained the defect class it was correcting; the countermeasure adopted here is to make the message echo observed data instead of asserting a generalisation.
- **Frequency claim corrected** from "nearly every real run" to "commonly" on the strength of a 20-fixture survey I had not done.
- **Deliberate deviations from review advice**, both recorded with reasons: emission stays in the parser (B option 10), and cardinality stays per-process despite the sibling differing (A §1.8, Felix's call).
- **Traps checked:** expected-fail control now covers three validations, not one; the byte-identity check uses the existing harness rather than an invented diff; once-ness is explicitly excluded from unit tests because the guard is process-global and already consumed by a sibling test in the same binary.
- **Remaining risk:** only the truncation budget. The change cannot alter output bytes, and validation 6 pins that with a harness that already exists.
