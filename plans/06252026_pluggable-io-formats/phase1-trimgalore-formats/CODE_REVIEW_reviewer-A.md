# Code Review — Phase 1 uBAM output (`--output-format ubam`) — Reviewer A

**Target:** working-tree diff (`git diff HEAD` + 3 untracked: `tests/integration_ubam_out.rs`, `test_files/ubam_out_se_REFERENCE.bam`, `test_files/ubam_out_pe_REFERENCE.bam`)
**Spec:** PLAN v2.1 (`.../phase1-trimgalore-formats/PLAN.md`)
**Reviewer:** A (independent; ran in parallel with B, no coordination)
**Date:** 2026-06-26

## Verification gate (re-run locally)

- `cargo build --release` — clean (only sandbox `xcrun` noise, not code).
- `cargo test --release` — **336 pass** (317 lib + 2 + 8 + 9 + 0 doc). Matches the implementer's claim.
- `cargo clippy --release --all-targets -- -D warnings` — clean.
- `cargo fmt --all -- --check` — clean.

## Summary

A well-structured, faithful implementation of PLAN v2.1. The architectural pivot (no trait, FASTQ path untouched, serial uBAM branch) is respected; the `_to_bam` trimmer entry points mirror the FASTQ ones line-for-line on stats/filtering; the `open_inner` refactor is safe at all four call sites; header propagation (`@HD`/`@PG`/`@CO` + appended `@PG PP:`) works as specified; empty-sequence and divide-by-zero edge cases are handled.

**One CRITICAL silent-data-corruption bug** undermines the headline feature (aux-tag preservation): the `--rename` clip-annotation collides with the §3.6 tag-tail encoding, corrupting the last preserved tag. This is reproducible on fully-allowed flag combinations and is not caught by any existing test. **Verdict: NEEDS_REWORK** (single targeted fix + a guard test; everything else is APPROVE-grade).

---

## CRITICAL

### C1 — `--rename` corrupts the last preserved aux tag on the uBAM-output path (silent data corruption)

**Where:** `src/trimmer.rs:261-274` (normal trim path) and `src/specialty.rs:130-138` / `:166-176` (hardtrim) call `FastqRecord::append_to_id(":clip5:…")` / `:clip3:…`. `src/fastq.rs:140-144` appends the suffix to the **end of the whole id**. `src/bam.rs:230-280` (`parse_name_and_data`) then splits the id on the first `\t` and treats *everything after* as the `\tTAG:TYPE:VALUE` tail.

**Effect:** when the input is uBAM with `--preserve-tags`, the id is `@NAME\tCB:Z:…\tUB:Z:GCTAGCTA`. `append_to_id` produces `@NAME\tCB:Z:…\tUB:Z:GCTAGCTA:clip5:AATTA`. Because `parse_tag_value` for `Z` keeps the rest via `splitn(3, ':')`, the **last tag's value is silently corrupted** with the rename annotation. Worse, this is written to the output BAM with no warning and looks like a legitimate (if odd) tag value downstream.

**Reproduced empirically** (built binary, committed fixture `test_files/ubam_test_with_tags.bam`):

```
$ trim_galore --clip_R1 5 --rename --output-format ubam --preserve-tags CB,UB \
      -o OUT test_files/ubam_test_with_tags.bam
$ samtools view OUT/ubam_test_with_tags_trimmed.bam | head -1
  ... CB:Z:ATCGATCG-1   UB:Z:GCTAGCTA:clip5:AATTA      <-- UB corrupted
```

Same corruption via `--hardtrim5 20 --rename --preserve-tags CB,UB` (the `:clip5:<45bp>` lands inside UB).

**Why it slips past tests:** every `parse_name_and_data` unit test exercises *either* a tag tail *or* a name with trailing text, never the *combination* of a real tag tail followed by an appended rename suffix. The integration tests never combine `--rename` with `--preserve-tags` on uBAM input.

**Why CRITICAL not IMPORTANT:** it is silent (exit 0, valid-looking BAM), it hits the exact feature this phase exists to deliver (tag fidelity), and `--clip_R1/R2 --rename` + uBAM + `--preserve-tags` is a fully-allowed combination — no CLI rule blocks it.

**Fix options (recommend, multiple valid approaches):**
1. **Insert the rename annotation into the *name*, before the tag tail.** Make `append_to_id` (or a uBAM-aware variant) splice `:clip5:…` immediately after the name token and before the first `\t`. Cleanest semantically — the annotation belongs to the read name, and round-trips through FASTQ identically.
2. **Reject `--rename` + `--output-format ubam`** at `Cli::validate()` (cheap, in line with the existing §3.4a rejections), deferring rename-into-BAM to a follow-up. Lowest-risk for v1.
3. Append the rename token only when there is no tag tail (detect `\t` in id first). Narrowest, but leaves rename silently dropped/unsupported when tags are present — least satisfying.

Whichever path: add an integration test combining `--rename` + `--preserve-tags CB,UB` on `ubam_test_with_tags.bam` and assert the tag value is exactly `GCTAGCTA`.

---

## IMPORTANT

### I1 — Re-running on trim_galore's own uBAM output produces a duplicate/self-referential `@PG` ID

**Where:** `src/bam.rs:177-188` (`build_output_header`) — `h.programs_mut().as_mut().insert("trim_galore".into(), pg)` and `build_pg_map` taking `last_pg_id(source)` as PP.

**Effect:** if the input BAM already has `@PG ID:trim_galore` (i.e. the user re-trims a trim_galore-produced uBAM), `insert` **overwrites** the existing entry rather than appending a new uniquely-ID'd one, and the replacement's `PP:` points at `trim_galore` — the ID it just replaced. Result is a self-referential `@PG` chain with a duplicate-by-intent ID. Reproduced:

```
$ trim_galore --output-format ubam -o OUT <a-trim_galore-output.bam>
$ samtools view -H OUT/...bam | grep ID:trim_galore
  @PG  ID:trim_galore  PN:trim_galore  VN:2.2.0  CL:...  PP:trim_galore   <-- dangling self-PP
```

The SAM spec requires unique `@PG` IDs. samtools tolerates it on read (it disambiguates), but emitting a duplicate ID with a dangling PP is malformed and loses the prior trim_galore step's provenance. PLAN §3.5 promised "provenance is preserved by adding to history, not preserving silence" — overwriting violates that exact intent.

**Fix:** when `source` already contains the chosen ID, disambiguate (`trim_galore.1`, `.2`, … — mirrors samtools' scheme) and keep PP pointing at the genuine last program ID. Low risk; isolated to `build_output_header`.

### I2 — `H:` (hex) tags round-trip asymmetrically: emitted on read, rejected on write

**Where:** read side `src/bam.rs:922-925` (`append_tag_type_and_value` emits `H:{s}` into the id tail) vs write side `src/bam.rs:307-310` (`parse_tag_value` **errors** on `H`). `B:` (array) is symmetric — rejected on both sides (read: `:926-933`, write: `:307`). But `H` is not.

**Effect:** a uBAM input carrying a hex tag listed in `--preserve-tags` reads fine (the `H:…` lands in the id tail), but `write_record` then **errors the whole run** at the first such record on the uBAM-output path. PLAN §3.6 documents B/H as not round-tripping, but the documented behavior is "stripped at read time" — for `B` that's true, for `H` it is not (it's emitted on read, then fatal on write). The mismatch turns a documented limitation into a hard run failure mid-stream.

This is unlikely in practice (hex tags are rare; CB/UB/RX are all `Z`), hence IMPORTANT not CRITICAL. **Fix:** either reject `H` on the read side too (symmetry, fail-fast at input), or have the write side skip/warn on `H` rather than abort. Recommend matching the read side's `B`-array treatment for consistency.

### I3 — Two PLAN §9 validation rows have no test

PLAN §9 lists both as test rows; neither exists:

- **`--cores N>1` + uBAM-out warns + still works.** No test asserts exit 0 + stderr "ignored". The only `--cores`-bearing uBAM test (`integration_ubam_out.rs:339`, `ubam_out_clumpify_rejected_at_cli`) expects a *rejection*, so it doesn't cover the warn-and-proceed path. (I verified the runtime behavior is correct — `main.rs:1369-1376` emits the NOTE and proceeds — so this is a coverage gap, not a bug.)
- **`--preserve-tags` mixed-batch allowed** (`preserve_tags_mixed_batch_allowed`, named in PLAN §3.4b note + §9). No such test. The §3.4b guard (`main.rs:148`, keyed on `any_bam`) is correct and I verified a `BS-seq_10K_R1.fastq.gz + ubam_test_with_tags.bam` SE batch runs and produces both BAMs — but the A-O1 loosening is untested.

**Fix:** add the two integration tests. Both pass against current behavior, so this is purely closing the promised matrix.

---

## NIT

### N1 — `_preserve_tags` field on `BamWriter` is dead

`src/bam.rs:76` stores `_preserve_tags: Vec<String>` (leading underscore, never read — the doc comment says the tail is the source of truth). It's an allocation + clone (`:109`) that serves only documentation. Consider dropping the field and keeping the rationale as a comment, or leave it as a deliberate symmetry marker. Cosmetic.

### N2 — `command_line` 4th param is an undocumented-in-signature PLAN deviation, but it IS documented inline

`BamWriter::create`'s 4th `command_line` param (PLAN §4 listed 3) is justified in the doc comment (`src/bam.rs:88-91`). Fine — flagging only because the review brief asked. No action.

### N3 — `@CO` line reorders relative to `@PG` block on propagation

uBAM-in → uBAM-out moves `@CO` from before the `@PG` block to after it. SAM header inter-group ordering is not semantically meaningful and `assert_ubam_eq` compares `@CO` independently, so harmless. Noting for completeness only.

### N4 — `build_pg_map` `.expect("…infallible…")`

`src/bam.rs:205-206` — acceptable; `Map::<Program>::build()` has no required fields. The expect message documents the invariant. Fine.

---

## Spec-conformance checks that PASSED (for the record)

- **§3.3 flag bits** — `0x04` SE / `0x4D` R1 / `0x8D` R2 verified by `bam_writer_flag_bits_correct` and re-confirmed via raw `samtools view -f`.
- **§3.3 qual** — Sanger→raw-Phred via `saturating_sub(33)`; `saturating_sub` is a sound defensive choice (FastqReader rejects sub-33; BamReader synthesises `!`=33). Round-trips byte-exact (`bam_writer_se_fastq_input_round_trip`).
- **§3.6 type fidelity** — `i`/`f`/`A`/`Z` land as real BAM types, not downgraded to `Z` (`bam_writer_aux_typed_int_and_float_round_trip`).
- **`open_inner` refactor** — all 4 call sites (`bam.rs:98,164,343,782`) destructure `(reader, header)` correctly; read path discards header, `peek_header` keeps it. No invariant lost; threaded closures still compile and pass.
- **Header propagation** — `@HD`/`@PG`/`@CO` preserved, `SO:`/`GO:` ordering tags preserved, new `@PG PP:` points at the genuine last program ID (verified on `ubam_test_with_tags.bam`).
- **§3.4a/b rejections** — all 6 rules (clumpify/passthrough/clock/implicon/demux/retain_unpaired at CLI; preserve-tags+all-FASTQ at format-detection) wired and tested (`cli.rs` unit tests + 3 integration rejection tests). `any_bam`-keyed §3.4b correctly allows mixed batches.
- **Empty sequence** — fully-trimmed read with `--length 0` writes a valid `*`-seq/`*`-qual FUNMAP record; no panic.
- **`pct()`** — divide-by-zero safe (`main.rs:1919`).
- **`@SQ`-bearing ("aligned in disguise") input** — rejected per-record by the pre-existing `bam_record_to_fastq` `is_unmapped()` check (`bam.rs:794`); the uBAM-output path reads through the same `RecordSource`, so aligned input is rejected before any BamWriter write. Good.
- **`parse_name_and_data` corners** — empty name (`@`, ``), missing type (`CB`, `CB:`), trailing tab, 2-char tag / 1-char type enforcement, space-vs-tab split all covered by lib tests and behave correctly.

---

## Verdict

**NEEDS_REWORK** — gated solely on **C1** (silent aux-tag corruption under `--rename`), which is a correctness defect in the phase's headline feature on an allowed flag combination. **I1** (duplicate `@PG` on re-processing) and **I2** (`H:` asymmetry) should be fixed in the same pass; **I3** closes the promised §9 test rows. Absent C1, this would be APPROVE_WITH_FIXES — the architecture, refactor, header handling, and the bulk of the test suite are solid.
