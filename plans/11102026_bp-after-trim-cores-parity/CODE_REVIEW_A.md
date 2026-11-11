# CODE REVIEW A — #541 `total_bp_after_trim` parity between the serial and parallel paths

**Reviewer:** A (independent; no contact with Reviewer B or the coverage audit)
**Target:** uncommitted working tree at `/Users/fkrueger/Github/TrimGalore`, branch `dev`, base `c0bbee4`
**Diff:** 3 files, +222/−5 — `src/parallel.rs`, `CHANGELOG.md`, `docs/src/content/docs/reference/changelog.md`
**Date:** 2026-11-10

---

## Summary

**The fix is correct and the two new tests are genuinely discriminating. No Critical or High findings.**

The three changed expressions are right, they are exhaustive, and nothing else in the crate is
perturbed. I confirmed the implementer's claims by execution rather than accepting them, including
the two the brief flagged as weakest — the hand-derived anchors and R2's discriminating power — and
both hold. The remaining findings are one Medium follow-up (a named state the suite still cannot
cover) and four Low items, three of which are wording.

Nothing in this change should block. The Medium is a recommendation for later, not a defect here.

### Verified by execution

| Check | Result |
|---|---|
| `cargo test`, isolated `CARGO_TARGET_DIR` | exit 0, **844 passed / 0 failed** (sum across all 21 targets) — matches the claim exactly |
| Both new tests by name | `..._parity_with_clipping` ok, `..._parity_paired` ok |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0, **0 warnings**; log line 176 confirms `Compiling trim-galore v2.4.0` — a real check, not a cache hit |
| `cargo fmt --all -- --check` | exit 0 |
| **SE anchor, independently** | Built a 2-read fixture with the test's exact sequences and ran the release binary at `--cores 1` (serial path, untouched by this diff, so a valid oracle): `Total written (filtered): 64 bp` for `43 + 34` input. Each clean sequence run alone: R1-clean **34 bp, 0 adapters**; R2-clean **33 bp, 0 adapters**. Therefore the adapter-bearing read's post-adapter length is **exactly 30**, and `15×30 + 15×34 = 960` is right. |
| **R2 discriminating power, independently** | Copied the crate to a scratch dir (`tar`, verified byte-identical), mutated `:543` alone to `r2.seq.len()`, ran `cargo test --lib bp_after_trim_parity`: the **paired test FAILED on R2 (990 vs 750)** while the SE test still passed. The shared tree was never touched — `git diff --stat` still reads 3 files / +222/−5 and `:543` is intact. |
| Changelog mirror | `diff` of the two bodies (temp files, not process substitution) → **exit 0, zero lines**: byte-identical |

The mutation run's full `TrimStats` dump confirms every derived R2 figure at once:
`total_reads: 30`, `total_reads_with_adapter: 0`, `total_bp_processed: 990`, `total_bp_written: 750`,
`reads_written: 30`. So R2 is 30 × 33 bp, picks up **no** adapter match, loses 8 bp to clipping
(33 − 5 − 3 = 25, × 30 = 750), and nothing is filtered.

### Verified by reading

Exhaustiveness, collateral scope, comment and changelog style, and the remaining coverage gaps —
detailed below with file:line.

---

## Logic

### L1 — The three changed expressions are correct. Confirmed.

`src/parallel.rs:542`, `:543`, `:1157` now add `res_r1.bp_after_cutadapt` / `res_r2.bp_after_cutadapt` /
`result.bp_after_cutadapt`. Each binding is in scope and already owned (`res_r1`/`res_r2` bound at
`:495`/`:496`, `result` at `:1133`); `bp_after_cutadapt` is a `usize`, so reading it does not disturb
the later read of `res_r2.clip_5prime_applied` at `:514`. Placement is unchanged: still after the
`--discard_untrimmed` early return and before `filters::filter_*`, matching the serial sites.

**Confidence: high** (read + compiled + tested).

### L2 — The three sites are exhaustive. Confirmed, with the clumpify and passthrough routes traced.

`grep total_bp_after_trim src/ tests/` returns exactly nine accumulation sites — `src/trimmer.rs:350`,
`:536`, `:537`, `:672`, `:780`, `:781` and the three changed ones — and all nine now add
`bp_after_cutadapt`. `src/parallel.rs` has exactly three `trim_read` call sites (`:495`, `:496`,
`:1133`), all inside `classify_paired` / `classify_single`.

The two routes the brief asked about do not bypass them:

- **Clumpify** differs only in *reader dispatch* — `read_pairs_clumpy` (`:812`) / `read_single_clumpy`
  (`:1237`) in place of the round-robin readers (`:669`, `:1202`). Both feed
  `process_paired_batch` (`:338`) / `process_single_batch` (`:1088`), which call `process_pairs`
  (`:602`) / `process_reads` (`:1187`) → `classify_*`. The two `process_pairs` call sites at `:386`
  and `:415` are the gzip and plain arms of one function, not two routes.
- **Passthrough** adds a sink argument to `process_pairs`; the passthrough slice is immutable and
  carries no stats of its own. Same `classify_paired`.
- `src/specialty.rs`, `src/demux.rs` and the `--clump_only` path contain no `total_bp_after_trim` at
  all (`--clump_only` writes `_clumping_report.txt`, and `ci.yml:1027` asserts the trimming report is
  *absent*).

**Confidence: high.**

### L3 — `total_bp_written` and the filter arms are untouched. Confirmed.

`total_bp_written` is written only in `FilterResult::Pass` / `PairFilterResult::Pass`
(`src/parallel.rs:552-553`, `:1170`), which the diff does not touch. The mutation run's dump shows
`total_bp_written: 750` identical on both sides, and the end-to-end JSON stayed at 594,545 per the
implementer.

One property worth stating because it is what keeps the report well-formed: `bp_after_cutadapt` is
captured at `src/trimmer.rs:165`, after quality trim and adapter trim, both of which only *truncate*.
So `total_bp_after_trim ≤ total_bp_processed` per read, and the percentage at `src/report.rs:636`
cannot exceed 100%. The fix raises the parallel figure; it cannot push it over the top.

**Confidence: high.**

### L4 — No golden fixture, CI md5 or test asserts the affected line.

`grep -rn "Total written" tests/ .github/` finds nothing. The four `*_trimming_report.txt` files under
`test_files/` are unreferenced leftovers — no test or CI step compares against them. The only CI steps
running `--cores >= 2` are the `--clump_only` family (ci.yml:986–1095, which md5 records, not reports)
and the clumpify memory job (ci.yml:439–449). So the brief's premise holds: CI could not have caught
this and will not now notice it changed.

**Confidence: high.**

### L5 — `run_paired_end` and `run_paired_end_parallel` treat an empty `adapters_r2` identically. Confirmed.

Both call `trimmer::trim_read(r2, config, true)` (`src/trimmer.rs:484`, `src/parallel.rs:496`), and the
`adapters_r2.is_empty() → adapters` fallback lives *inside* `trim_read` (`src/trimmer.rs:113-117`), so
the two paths share it by construction rather than by agreement. The only other `is_r2` branch before
the `bp_after_cutadapt` capture at `:165` is the directional-RRBS arm, and `rrbs` is `false` in this
config.

Note for the record: R2 *is* scanned with the Illumina adapter here (`adapters_r2` empty → falls back),
it simply does not match. `total_reads_with_adapter: 0` in the mutation dump is the proof, and it is
what makes the R2 side a clean single-variable test of clipping.

**Confidence: high.**

---

## Efficiency

No findings. The change swaps a `Vec`-length read for a `usize` field read on an owned struct, per
read. No allocation, no branch, no extra pass. `clippy --release -D warnings` is clean.

---

## Error handling

No findings. The change introduces no fallible operation. Both new tests propagate with `?` and return
`Result<()>`, matching the module's idiom; `fs::remove_dir_all(&dir).ok()` on the success path only
(so a failing run leaves its fixture for inspection) is the idiom already used by every neighbouring
test.

`fresh_tmpdir` (`:1414`) uses a fixed path under `std::env::temp_dir()` and `remove_dir_all`s it on
entry, so two tests sharing a slug would race under cargo's parallel runner. I checked all 17 slugs in
the module: both new ones (`tg_parallel_serial_bp_after_trim`, `tg_parallel_serial_bp_paired`) are
unique and neither nests inside another. **Not a finding** — just confirming the hazard was avoided.

---

## Structure

### S1 — The discrimination guard is stronger than the plan asked for, and satisfiable one way it shouldn't be. (Low)

`src/parallel.rs:1609-1613` and `:1715-1720` guard with
`total_bp_after_trim > total_bp_written`, where PLAN §7 suggested `< total_bp_processed`. The
implementer's choice is the better one: adapter trimming alone satisfies the plan's version, so it
would have been vacuous, whereas removing the clips makes this one fire (960 > 960 is false).

The residual weakness: `total_bp_written` counts only *passing* reads, so the inequality is also
satisfied by any read being filtered out, with no post-capture trimming at all. Today that cannot
happen — `length_cutoff: 0` is inherited from `parity_test_config()` and the dump confirms
`reads_written: 30`, nothing filtered — so the guard means what its message says. But a future edit
that raised `length_cutoff` enough to drop reads *and* removed the clips would leave the guard green
while the test went vacuous, which is the exact failure class the guard exists to prevent.

A guard that cannot drift: assert the config's clip fields are `Some` directly, or compare against the
per-read written length of passing reads.

**Confidence: high** that the loophole exists; **low** that anyone walks into it.

### S2 — The new paired test omits `PairValidationStats`, which the test it mirrors treats as load-bearing. (Low)

`run_paired_end` returns `(TrimStats, TrimStats, PairValidationStats)` (`src/trimmer.rs:396`). The new
test asserts `.0` and `.1` only (`src/parallel.rs:1693`, `:1699`). The existing paired parity test
asserts `.2` **first**, with a comment recording it as a prior reviewer's critical catch
(`src/parallel.rs:2523-2531`), and `PairValidationStats` derives `PartialEq` specifically for that
assertion (`src/report.rs:107-110`).

Nothing in this change touches pair counters, so this costs nothing today. It matters because the new
test is now the module's cheaper, more focused paired serial-vs-parallel harness — the one a future
author is likelier to extend — and it is one line weaker than the one it was modelled on. Same for
output records: the test writes four FASTQ files and never reads them back, where `read_records` is
already in the module.

Recommendation: add `assert_eq!(serial.2, parallel.2, ...)`.

**Confidence: high** (read).

### S3 — Fixture duplication between the two new tests and the existing parity test. (Low, no action recommended)

The SE fixture loop at `:1551-1565` is character-for-character the loop at `:1480-1494`, and the PE
test duplicates the harness again. A shared builder would be DRY-er, but the duplication is what makes
each hand-derived anchor auditable where it sits, which is the point of PLAN §6. I would leave it. I
mention it only so the next reviewer does not read the duplication as an oversight.

---

## Comments and changelog style

### C1 — The new source comment is compliant, but names one of five post-capture stages. (Low)

`// Post-adapter, pre-RRBS length, so the cutadapt section matches Perl v0.6.x (#232).`
(`src/parallel.rs:541`, `:1156`)

One line, current state, no narration of the change, invariant stated on the line it constrains — this
is exactly what CLAUDE.md asks for, and it is a clear improvement on what it replaced.

The nit is accuracy. `bp_after_cutadapt`'s own doc (`src/trimmer.rs:64-71`) lists **five** things the
capture precedes: the RRBS 2 bp adjustment, poly-A/G trimming, N-trimming and fixed clipping. The
comment names only RRBS — while this change's own tests and assertion messages say "post-adapter,
**pre-clip**" (`:1542`, `:1600`, `:1604`), and the defect that prompted the fix is clip-visible. So one
quantity carries three phrasings inside one change, and a maintainer who reads `pre-RRBS` at `:542`
and then finds a *clipping* test guarding it has to go read `trimmer.rs` to reconcile them.

`// Post-adapter length, so the cutadapt section matches Perl v0.6.x (#232).` would be accurate for all
five stages and one word shorter. Cosmetic; worth doing while the lines are being touched anyway.

For context, not as a finding: the three serial sites (`src/trimmer.rs:346-348`, `:533-535`) still
carry the older three-line comments that do narrate ("isn't double-subtracted vs Perl v0.6.x"), so the
two families now differ in both length and wording. PLAN Step 2 deliberately left them alone, and I
agree — widening the diff to reword correct serial comments is not this change's job.

### C2 — The CHANGELOG entry matches the file but runs past the one-sentence rule. (Low, Felix's call)

The entry (`CHANGELOG.md:6-13`) is a bold title plus one long sentence. Against the file's own
established format it is unremarkable — the v2.4.0 entries above it are all bold-title + multi-sentence
paragraphs, and several are longer. Against CLAUDE.md's "one sentence per entry, no exceptions" it is
two, and the second carries the internal field name `total_bp_after_trim` and the serial/parallel
mechanism, which is audit-trail material: a reader deciding whether an upgrade affects them needs the
*flag list*, not the field that was wrong.

The flag enumeration earns its place — it is the only way a user knows whether their runs moved. The
mechanism clause ("The parallel path accumulated the final record length into `total_bp_after_trim`
while the serial paths used the post-adapter, pre-clip length introduced for #232") is what the commit
message is for. Cutting it leaves roughly one sentence and loses nothing a user needs.

I am flagging the tension, not asserting the entry is wrong: the global rule says to follow the file's
existing format exactly, and this file's format is paragraphs. That conflict is Felix's to resolve.

### C3 — The docs mirror is correct and was necessary. Confirmed.

The two changelog bodies are now byte-identical (`diff` exit 0, zero lines). There is **no** automated
sync: `.github/workflows/docs.yml:8` lists `CHANGELOG.md` only as a *rebuild trigger* for the Astro
site, and `docs/scripts/` holds only chart and logo generators. The page's own preamble
(`docs/src/content/docs/reference/changelog.md:6-8`) declares itself a synced copy. So PLAN deviation
#2 is sound and the root-only edit really would have left the page desynced.

---

## What the tests still cannot catch

Named states, not generic gaps.

### G1 — `rrbs: true` is set by no test in `src/parallel.rs`. (Medium — follow-up, not a defect here)

`grep rrbs src/parallel.rs` returns only `rrbs: false` at `:1443`. So:

- The `--rrbs` row of the issue's own table (490,527 bp vs 483,967 bp) is verified only by the
  implementer's manual end-to-end run. Nothing locks it.
- `rrbs_r2_clipped_5prime` (`src/parallel.rs:515`) and `rrbs_trimmed_3prime` / `_5prime`
  (`:502-512`, `:1136-1141`) have **no** serial-vs-parallel parity coverage at all — and
  `rrbs_r2_clipped_5prime` is one of the two fields commit `3996fc5` was written to fix, per the
  comment at `:1470`.

Clipping is a structurally equivalent discriminator (same capture point, same accumulation line), so
the #541 fix itself is adequately locked. What is missing is coverage of the regime #232 was actually
about, in the paired RRBS shape where `clip_r2` is auto-set to 2 inside `trim_read`.

Cheapest closure: a third variant of the new test with `rrbs: true, is_paired: true` over a fixture
whose reads are adapter-bearing (the directional 3' trim is gated on `had_adapter`,
`src/trimmer.rs:195`) — roughly the existing paired test with two config fields changed and one
recomputed anchor.

**Confidence: high** that the gap is real; it is pre-existing and out of this change's scope.

### G2 — `--discard_untrimmed` has no behavioural coverage anywhere in the repo. (Medium — pre-existing, arguably its own issue)

No unit test sets `discard_untrimmed: true`; no integration test or CI step passes the flag; nothing
asserts `discarded_untrimmed`. The only reference outside the CLI definition is the `--clump_only`
rejection at `src/cli.rs:886`.

Why it belongs in this review: the accumulation this change edits sits *immediately after* the
`--discard_untrimmed` early return (`src/parallel.rs:535-540`, `:1152-1155`), and PLAN §3 lists
"contributes to neither field" as an edge case. An edit that moved the accumulation above that early
return would change `total_bp_after_trim` on both paths **equally**, so every `serial == parallel`
assertion in the file stays green — and the new absolute anchors stay green too, because this fixture
discards nothing. That is precisely the "test that can only ever pass" class this plan was written to
eliminate, and it is the one remaining way to break this field silently.

Per the repo's issue-filing bar this is a review spin-off with no user-visible defect behind it, so it
belongs in the plan doc rather than a new issue unless Felix wants it tracked.

**Confidence: high** (grepped `src/`, `tests/`, `ci.yml`).

### G3 — The other three post-capture stages. (Low)

`--polyA`, `--polyG` and `--trim-n` are the remaining things the capture precedes, and none is
exercised on this axis. They share the capture point with `clip_*`, so the marginal value is small —
except that poly-G now runs ahead of poly-A per #504, and that reordering moves bases between the two
counts while leaving `bp_after_cutadapt` alone. A parity test over a poly-X config would pin that
interaction. Low priority.

### G4 — Reads fully consumed by the adapter. (Low)

`bp_after_cutadapt == 0` adds zero; the fixture has no such read (every adapter-bearing read keeps
30 bp). Arithmetically trivial, but it is the one input shape where the two expressions coincide at
zero rather than diverging, so a fixture of only such reads would be vacuous. Worth a sentence in the
test doc comment at most.

### G5 — The uBAM output path at `--cores >= 2`. (Low, out of scope)

The uBAM *serial* sites (`src/trimmer.rs:672`, `:780`, `:781`) are already correct and this change does
not touch them. Per the brief, #529 exposes the same divergence on the uBAM output path; since its
`process_{single,paired}_batch_bam` route through `classify_*`, it inherits this fix on rebase — but
that inheritance is an assertion about an unmerged branch, and nothing here tests it. Worth one line on
#529 rather than work in this change.

---

## Recommendations

**Critical** — none.

**High** — none.

**Medium**
1. **G1** — add an `rrbs: true, is_paired: true` variant of the new paired parity test, which would also
   give `rrbs_r2_clipped_5prime` its first serial-vs-parallel coverage. Follow-up, not a blocker.
2. **G2** — give `--discard_untrimmed` one behavioural test with an absolute anchor, so the accumulation
   *placement* is pinned the way this change pinned the accumulation *value*. Pre-existing; record in the
   plan doc rather than filing.

**Low**
3. **C1** — reword the two new comments to `// Post-adapter length, so the cutadapt section matches Perl
   v0.6.x (#232).` — accurate for all five post-capture stages, and consistent with the "pre-clip"
   wording used elsewhere in the same change.
4. **S2** — add `assert_eq!(serial.2, parallel.2, ...)` to the new paired test, matching
   `src/parallel.rs:2525`.
5. **S1** — make the discrimination guard independent of filtering (assert the clip fields are `Some`, or
   compare against passing reads' written length).
6. **C2** — consider cutting the CHANGELOG entry's mechanism clause to the commit message, keeping the
   bold title, the issue link and the affected-flag list. Judgment call given the file's established
   paragraph format.

None of 3–6 changes behaviour; all are one-or-two-line edits if Felix wants them before commit.

---

## Housekeeping

- I edited no source file. The only file I created is this report.
- No commit, push, branch switch or GitHub write.
- The one mutation check I ran was performed on a `tar` copy under my scratchpad, never on the shared
  tree. Re-verified afterwards: `git diff --stat` = 3 files / +222/−5, `src/parallel.rs:543` intact,
  branch `dev`.
- Builds and tests used `CARGO_TARGET_DIR` under
  `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/aecc6671-429c-462c-ad9c-38f8e4d35fbe/scratchpad`;
  the repository's `target/` was read but not written.
- Logs: `revA_cargo_test.log`, `revA_clippy.log`, `revA_fmt.log`, `revA_mut543.log` in that scratchpad
  directory.
