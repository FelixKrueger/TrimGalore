# PLAN — `total_bp_after_trim` parity between the serial and parallel paths

**Issue:** [#541](https://github.com/FelixKrueger/TrimGalore/issues/541)
**Branch base:** `dev` @ `c0bbee4`
**Date:** 2026-11-10
**Revision:** 2 — see §12. Revision 1 contained three false factual claims and specified a test that
could not fail; both were caught by dual plan review.

---

## 1. Goal

Make the trimming report's cutadapt-section `Total written (filtered)` line independent of
`--cores`, by having the parallel path accumulate `TrimResult::bp_after_cutadapt` — as the serial
paths already do — instead of the final record length.

Secondary, and where the real risk lies: close the test gap that let this survive for six months, in
a way that cannot silently reopen.

Out of scope: #542 (`MultithreadedWriter`).

## 2. Context

`TrimResult::bp_after_cutadapt` (`src/trimmer.rs:71`, documented `:64-71`, captured `:165`) holds the
read length after quality and adapter trimming but **before** the RRBS 2 bp adjustment, poly-A/G
trimming, N-trimming and fixed clipping. It exists so the cutadapt-section line reports what Perl
v0.6.x's cutadapt subprocess wrote rather than what the wrapper then truncated — issue #232.

Six sites honour the contract, on **two distinct serial paths** (revision 1 conflated them):

| Path | Sites | Function |
|---|---|---|
| FASTQ serial | `src/trimmer.rs:350` | `run_single_end` |
| FASTQ serial | `src/trimmer.rs:536`, `:537` | `run_paired_end` |
| **uBAM serial** | `src/trimmer.rs:672` | `run_single_end_to_bam` (`:638`) |
| **uBAM serial** | `src/trimmer.rs:780`, `:781` | `run_paired_end_to_bam` (`:708`) |

The parallel path uses the final record length at three sites:

| Site | Current expression | Binding in scope |
|---|---|---|
| `src/parallel.rs:542` | `r1.seq.len()` | `res_r1` (bound `:495`) |
| `src/parallel.rs:543` | `r2.seq.len()` | `res_r2` (bound `:496`) |
| `src/parallel.rs:1157` | `record.seq.len()` | `result` (bound `:1133`) |

These are the **only** three sites: `parallel.rs` has exactly three `trim_read` call sites, and the
`--clumpify` and `--passthrough` paths funnel through the same `classify_single` / `classify_paired`
functions rather than accumulating independently.

**`TrimResult` is already in scope at all three sites**, so this is a value change, not plumbing.

**Accumulation placement is already correct and must not move.** Both path families accumulate after
the `--discard_untrimmed` early exit and before `filters::filter_*`, so discarded reads contribute to
neither field while length-filtered reads contribute to `total_bp_after_trim` but not
`total_bp_written`. `TrimStats::merge` (`src/report.rs:70`) is a plain field-wise sum, so a per-read
change is at the right level and batch merging cannot double-count.

**`total_bp_written`** (`src/report.rs:21`) is a separate field, written only in the filter `Pass`
arms, already correct on every path, and must come out of this change unchanged. That is why the JSON
report shows no divergence.

### Root cause — an orphaned call site

| Date | Commit | Effect |
|---|---|---|
| 17 Apr 2026 | `82d1e34` | Added `total_bp_after_trim` to the parallel SE/PE paths using `seq.len()`. **Correct at the time** — `bp_after_cutadapt` did not yet exist. |
| 27 Apr 2026 | `f4b544d` | Introduced `bp_after_cutadapt` for #232 and switched the serial sites. Touched `CHANGELOG.md`, `src/report.rs`, `src/trimmer.rs` — **not `src/parallel.rs`**. |

The parallel code was right when written; the definition moved underneath it ten days later. This is
an orphaned call site, not a dropped comment clause.

`f4b544d`'s own message says it uses `bp_after_cutadapt` "in both single-end and paired-end
pipelines". That phrase is ambiguous — SE/PE is a different axis from serial/parallel — but a reader
would reasonably take it as full coverage, and it has not held for `--cores >= 2` since April.

### Why the tests missed it

**Two** serial-vs-parallel `TrimStats` parity tests exist, and both are vacuous for this field:

| Test | Line | Covers |
|---|---|---|
| `test_parallel_serial_trim_stats_parity` | `:1461` | single-end, whole-struct |
| `test_passthrough_serial_parallel_parity` | `:2271` | **paired, R1 and R2**, asserted at `:2338-2339` |

Both build their input from `parity_test_config()` (`:1424`), which sets `trim_n: false`,
`clip_r1: None`, `clip_r2: None`, `three_prime_clip_r1: None`, `three_prime_clip_r2: None`,
`rrbs: false`, `poly_a: false`, `poly_g: false`. `passthrough_test_config()` (`:2168`) inherits it
via `let mut c = parity_test_config()`, changing only `length_cutoff` and `is_paired`. Under that
configuration `bp_after_cutadapt == record.seq.len()` by construction —
`test_bp_after_cutadapt_equals_seq_len_when_no_rrbs` (`src/trimmer.rs:1102`) asserts exactly that —
so neither test's assertion on this field can fail.

Separately, the CI Perl byte-identity matrix runs `--cores 1` (one step is labelled `--cores 1
(sequential path)`), so the invariant is only ever checked against Perl on the paths that implement
it.

### Measured, on `c0bbee4`

| Run | `--cores 1` | `--cores 4` | delta |
|---|---|---|---|
| `--clip_R1 5`, `BS-seq_10K_R1.fastq.gz` | 644,525 bp (99.2%) | 594,545 bp (91.5%) | 49,980 bp |
| `--rrbs`, `SRR24766921_RRBS_R1.fastq.gz` | 490,527 bp (96.2%) | 483,967 bp (94.9%) | 6,560 bp |

Decompressed output byte-identical in both cases (9,996 and 9,729 records) — reporting-only.

**The `--clip_R1` delta is `5 × 9,996 = 49,980`**: five bases per *written* read, with four reads
unclipped because `clip_5prime` is a no-op at or below the clip length (`src/fastq.rs:111`).
Revision 1 said "10,000 × 5 less the reads filtered out", which is the right number by the wrong
mechanism — filtered reads contribute to this field on both paths, since accumulation precedes
filtering.

### Affected population

`--clumpify` **requires `--cores >= 2`** (`src/cli.rs:745`), so every clumpify run takes the parallel
path and reports the wrong figure. That is likely the largest affected group, ahead of users who pass
`--cores` explicitly.

## 3. Behavior

1. For every read reaching the accumulation point, the parallel path adds
   `result.bp_after_cutadapt`, the same quantity the serial paths add.
2. `Total written (filtered)` is identical at `--cores 1` and `--cores N` for all N, for every
   combination of `--rrbs`, `--clip_R1/R2`, `--three_prime_clip_R1/R2`, `--polyA`, `--polyG` and
   `--trim-n`.
3. `total_bp_written`, `total_bp_processed`, every other `TrimStats` field, and the trimmed records
   are unchanged.

### Edge cases

| Case | Expected |
|---|---|
| `--discard_untrimmed` discards a read | Contributes to neither field, all paths |
| Read filtered `TooShort` / `TooManyN` | Contributes to `total_bp_after_trim`, not `total_bp_written` |
| Read at or below the clip length | `clip_*` is a **no-op** (`src/fastq.rs:111`, `:126`), so the two expressions coincide and the read adds no divergence |
| Read fully consumed by adapter trim | `bp_after_cutadapt == 0`; adds 0 |
| `--rrbs` on R2 where the 5' clip applies | Captured before the clip; `rrbs_r2_clipped_5prime` unaffected |
| Empty input / single record | Both paths equal |
| `--cores 1` | Worker pool bypassed (`src/cli.rs:358`); unaffected |

## 4. Signature

No signature changes. Three expressions change; one test helper is added.

## 5. Implementation outline

### Step 1 — fix the three accumulation sites

`src/parallel.rs:542-543`:

```rust
stats_r1.total_bp_after_trim += res_r1.bp_after_cutadapt;
stats_r2.total_bp_after_trim += res_r2.bp_after_cutadapt;
```

`src/parallel.rs:1157`:

```rust
stats.total_bp_after_trim += result.bp_after_cutadapt;
```

### Step 2 — state the invariant in the comment

The parallel comments currently read "Track bp after trimming but before length filtering
(Cutadapt-compatible stat)", which was accurate when written and became incomplete when `f4b544d`
redefined the quantity. Replace each with one line naming the invariant, so the next reader sees what
constrains the expression:

```rust
// Post-adapter, pre-RRBS length, so the cutadapt section matches Perl v0.6.x (#232).
```

One line, current state, no narration of the change. Leave the serial comments alone.

### Step 3 — add a discriminating config, without touching the shared one

`parity_test_config()` has **12 callers** (`:1483, 1553, 1590, 1620, 1663, 1664, 1758, 1826, 1909,
1980, 2025, 2169`). Do **not** enable trimming inside it.

Add a sibling helper beside it. **It must set the R2 values as well as R1** — revision 1 set only R1
in the helper and mentioned R2 in prose, which would have left `:543` undiscriminated and the paired
test passing however `:543` was written:

```rust
/// Clipping on both reads and both ends, so `bp_after_cutadapt` and the final
/// record length differ for every read longer than the clip lengths.
fn parity_test_config_post_adapter() -> TrimConfig {
    TrimConfig {
        clip_r1: Some(5),
        clip_r2: Some(5),
        three_prime_clip_r1: Some(3),
        three_prime_clip_r2: Some(3),
        ..parity_test_config()
    }
}
```

`clip_*` is the right discriminator: it applies to every read above the clip length regardless of
sequence content, so the divergence is arithmetically predictable rather than data-dependent.

### Step 4 — strengthen the single-end parity test

Extend `test_parallel_serial_trim_stats_parity` (`:1461`), or add a sibling, to repeat its comparison
under `parity_test_config_post_adapter()`. Keep the existing clean-config run — it still guards the
other fields.

### Step 5 — extend the paired parity test

`test_passthrough_serial_parallel_parity` (`:2271`) already drives serial `run_paired_end` against
`run_paired_end_parallel` and asserts R1 and R2 `TrimStats` at `:2338-2339`. Add a run under a
post-adapter config (`parity_test_config_post_adapter()` with `length_cutoff` and `is_paired` set as
`passthrough_test_config()` does) rather than building a new harness.

Assert **both** R1 and R2 — `:542` and `:543` are independent sites writing to different `TrimStats`
instances, so a fix applied to one must not pass.

### Step 6 — anchor the expected value, not just the agreement

Every assertion in Steps 4 and 5 is `serial == parallel`, which pins the two paths together but not
to anything true. A later change flipping the **serial** line to `seq.len()` would keep the whole
suite green and silently revert #232 everywhere.

Add at least one absolute assertion — on a small fixture with known lengths, assert the literal
expected `total_bp_after_trim`, derived by hand from post-adapter lengths and independent of either
implementation.

### Step 7 — guard the fixture's discriminating power

Both replacement tests must assert that their configuration still makes the two quantities differ,
or a later fixture edit returns them to "cannot fail" — the exact defect being repaired. Follow the
idiom already present in the test being strengthened (`src/parallel.rs:1518-1524`), which
sanity-checks that the fixture exercised adapter trimming "otherwise we'd be passing trivially":

```rust
assert!(
    stats_serial.total_bp_after_trim < stats_serial.total_bp_processed,
    "fixture must exercise post-adapter trimming, or this test cannot fail"
);
```

### Step 8 — CHANGELOG

Add a one-sentence entry under `### Unreleased` (currently empty), in the file's existing format.
One sentence, no migration prose, per the project's changelog rule.

### Step 9 — confirm the tests are discriminating

Revert each changed line **separately** and confirm a test fails each time (§9, validations 5–6).
Do not skip this: a test that cannot fail is the defect this plan exists to repair, and revision 1
of this plan specified one.

## 6. Efficiency

No change. `bp_after_cutadapt` is a `usize` already computed and stored in `TrimResult` for every
read on every path; reading it costs nothing that `seq.len()` did not. No allocation, no extra pass,
no complexity or memory change.

## 7. Integration

**Reads:** `TrimResult::bp_after_cutadapt`, already populated.
**Writes:** `TrimStats::total_bp_after_trim` only.
**Downstream:** `src/report.rs:635-636` formats the value and its percentage of
`total_bp_processed` into the text report. The JSON report carries `total_bp_written`, untouched.

**Perl parity:** this moves the parallel path *toward* v0.6.x parity. The CI matrix runs `--cores 1`
and stays green either way, so **it proves nothing about this change** — Steps 4–7 are the only guard.

**#529:** it adds **no accumulation of its own**. The uBAM output path already exists on `dev` and
accumulates correctly on its serial route (`:672`, `:780`, `:781`); #529's
`process_{single,paired}_batch_bam` call `classify_single` / `classify_paired`, so it **inherits this
fix on rebase**. Revision 1 claimed the opposite and would have sent someone hunting a nonexistent
site. The note on #541 needs the same correction.

**No golden-file or md5 fixture changes.** No test asserts a `total_bp_after_trim` value, so nothing
needs updating.

## 8. Assumptions

1. `bp_after_cutadapt` is the correct target semantics, grounded in #232 and the contract at
   `src/trimmer.rs:64-71`. **Fixed rule.**
2. `trim_read` mutates the record in place, so `seq.len()` after the call is the final length.
   Verified.
3. Accumulation placement is already correct on all paths and is not changed. Verified.
4. `total_bp_written` is intentionally the final length and stays as it is. **Fixed rule.**
5. `parity_test_config()` stays byte-for-byte as it is; 12 tests depend on it.
6. `clip_*` **skips** reads at or below the clip length rather than clamping them
   (`src/fastq.rs:111`, `:126`), so the divergence is `clip_total × (reads above the clip length)`,
   not `clip_total × all reads`. Revision 1 assumed unconditional application, which was false.
7. `--cores 2` suffices to exercise the worker pool for this field. **Configurable.**

## 9. Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| 1 | Single-end figure no longer varies with `--cores` | `--clip_R1 5` on `BS-seq_10K_R1.fastq.gz` and `--rrbs` on `SRR24766921_RRBS_R1.fastq.gz`, at `--cores 1` and `4` | 644,525 bp (99.2%) and 490,527 bp (96.2%) at **both** core counts |
| 2 | **Paired** figure no longer varies | `--paired --clip_R1 5 --clip_R2 5` on `BS-seq_10K_R{1,2}.fastq.gz` at `--cores 1` and `4`; compare **both** report files | Identical across core counts for R1 *and* R2 |
| 3 | Records untouched | `gzip -dc`, compare SHA-256 against pre-fix output | Identical (9,996 / 9,729 records single-end) |
| 4 | `total_bp_written` unchanged | Compare `basepair_processing.total_bp_written` before and after | 594,545 in all cases |
| 5 | Single-end site genuinely covered | Revert `:1157` alone, run the suite | Step 4's test fails |
| 6 | **Each** paired site genuinely covered | Revert `:542` alone, then `:543` alone | Step 5's test fails **each** time |
| 7 | The anchor catches a serial regression | Flip `src/trimmer.rs:350` to `record.seq.len()` | Step 6's absolute assertion fails |
| 8 | No collateral damage | `cargo test` from the crate root | 842 tests pass, plus the new ones |
| 9 | CI-clean | `cargo fmt --all -- --check`; `cargo clippy --all-targets --release -- -D warnings` | Clean — local rustc 1.95 vs CI 1.99, so necessary but not sufficient |

Validations 6 and 7 are the ones most likely to be skimped and the ones that catch, respectively, a
half-applied fix and a future silent reversion.

## 10. Questions or ambiguities

**Open (assumption taken, does not block):**

- Whether Steps 4 and 5 extend the existing tests or add siblings. Implementer's choice; the existing
  clean-config runs must survive either way.
- Which fixture carries the Step 6 absolute anchor. A small purpose-built one with hand-computed
  lengths is likelier to stay correct than a derived figure from a 10,000-read file.
- Whether MultiQC parses `Total written (filtered)`. Unverified, noted as such on #541. It changes
  priority, not the fix.
- Whether the #232 CHANGELOG wording warrants its own correction, given the SE/PE vs serial/parallel
  ambiguity. Taken: mention in this change's entry, do not rewrite history.

**No critical ambiguities.**

## 11. Self-Review

**Efficiency.** No complexity, allocation or memory change.

**Logic.** All three parallel sites traced against their serial counterparts; placement confirmed
already identical, so changing the value alone is sufficient. `TrimStats::merge` confirmed a plain
sum. Confirmed the three sites are exhaustive, including the clumpify and passthrough routes.

**Edge cases.** §3 covers discarded reads, filtered reads, reads at or below the clip length,
zero-length reads, the RRBS R2 clip interaction, empty input, single record, and `--cores 1`.

**Integration.** Only consumer is `src/report.rs:635-636`; JSON path unaffected; no fixture moves.

**Remaining risks.**

1. CI cannot catch a regression here — it runs `--cores 1`. Steps 4–7 are the whole guard.
2. The absolute anchor in Step 6 must be computed by hand. Deriving it from a run of either path
   would reproduce the defect it exists to prevent.
3. The Step 7 guards are themselves assertions that could be weakened by a later edit. This bottoms
   out somewhere; the guard makes the next such edit visibly wrong rather than silently inert.

## 12. Implementation notes

Implemented 2026-11-10 against `dev` @ `c0bbee4`. All nine steps done. Diff confined to three files,
per `git diff --numstat`: `src/parallel.rs` (+202/−5), `CHANGELOG.md` (+10),
`docs/.../reference/changelog.md` (+10) — 222 insertions, 5 deletions in total.

| Step | Outcome |
|---|---|
| 1 — three accumulation sites | `src/parallel.rs:542`, `:543`, `:1157` now add `bp_after_cutadapt`. All nine sites across both files verified consistent afterwards. |
| 2 — comments | One-line invariant at both sites. |
| 3 — helper | `parity_test_config_post_adapter()` with **all four** clip values. |
| 4 — single-end test | New sibling `test_parallel_serial_bp_after_trim_parity_with_clipping`. |
| 5 — paired test | New sibling `test_parallel_serial_bp_after_trim_parity_paired` (see deviation below). |
| 6 — absolute anchor | Hand-derived `15×30 + 15×34 = 960` (SE), `960` / `30×33 = 990` (PE R1/R2). |
| 7 — discrimination guard | `total_bp_after_trim > total_bp_written` in both tests. |
| 8 — CHANGELOG | One sentence under `### Unreleased`, mirrored to the docs copy (see deviation below). |
| 9 — mutation checks | All four performed; see below. |

### Deviations

1. **Step 5 added a focused paired test rather than extending
   `test_passthrough_serial_parallel_parity`.** §10 delegated extend-vs-sibling to the implementer.
   Extending that test would have coupled #541 coverage to passthrough, doubled a 12,000-row
   fixture's runtime, and made Step 6's hand-derived anchor impractical against generated multiome
   rows. The new test drives the same entry points (`run_paired_end` /
   `run_paired_end_parallel`) with `None` for the passthrough slots. The existing passthrough test
   is untouched and still passes.
2. **Step 8 also updated `docs/src/content/docs/reference/changelog.md`**, which the plan did not
   mention. That page is a hand-maintained byte-identical copy and states so in its own preamble;
   the root-only edit had desynced it (12 lines vs 2). Mirrored rather than left to drift.
3. **Step 7's guard asserts `total_bp_after_trim > total_bp_written`, not the plan's
   `< total_bp_processed`.** The plan's snippet would have passed under the *clean* config too
   (960 < 1155), so it would not actually have guarded discrimination — the one job Step 7 exists
   to do. Comparing against `total_bp_written` is false exactly when no post-adapter trimming
   occurred. Sound here because `length_cutoff: 0` filters nothing; a filtering fixture would need
   a different formulation.
4. **Validation 7 was run as a two-line mutation**, flipping the serial *and* parallel single-end
   sites, rather than the plan's one-line serial flip. The one-line version fails on the parity
   assertion, which proves nothing about the anchor; flipping both makes parity pass so that only
   the anchor can fire. Strictly more informative, and it is what demonstrates Step 6's purpose.
5. **The paired test leaves `length_cutoff` at 0**, inherited from `parity_test_config()`, so no
   read is filtered and the anchors are exact. Incidental rather than chosen.

### Mutation checks (validations 5–7)

Each revert was applied alone and reverted immediately after.

| Mutation | Result |
|---|---|
| `:1157` → `record.seq.len()` | SE test **failed**: one field differs, `960` vs `720`, all others identical |
| `:542` → `r1.seq.len()` | PE test **failed on R1**: `960` vs `720` |
| `:543` → `r2.seq.len()` | PE test **failed on R2**: `990` vs `750` — the site revision 1 could not have caught |
| **`trimmer.rs:350` *and* `:1157`** → `seq.len()` | Parity assertion **passed** (both agree at 720); the **absolute anchor** fired. Confirms Step 6 catches the silent #232 reversion no `serial == parallel` assertion can. |

The fourth is the one that justifies Step 6's existence, and it behaved exactly as intended.

### Verification

| Check | Result |
|---|---|
| Validation 1 (SE figure) | `644,525 bp (99.2%)` and `490,527 bp (96.2%)` at **both** `--cores 1` and `4` |
| Validation 2 (PE figure) | R1 `644,525`, R2 `645,073`, identical across core counts; R1 ≠ R2 so both sides are genuinely exercised |
| Validation 3 (records) | Byte-identical to pre-fix output — `53d619cd…` (9,996 records), `f49dea91…` (9,729) |
| Validation 4 (`total_bp_written`) | `594,545` before and after, at both core counts |
| Validation 8 (suite) | `cargo test` exit 0 — **844 passed, 0 failed** (842 + 2 new) |
| Validation 9 (CI gates) | `clippy --all-targets --release -D warnings` exit 0; `fmt --check` exit 0 after one `cargo fmt` pass that touched only the new test |

All three figures computed by hand in §2 were confirmed by the mutation output:
`total_bp_processed: 1155`, `total_bp_written: 720`, `total_bp_after_trim: 960`.

**Not done, deliberately:** no commit, branch or PR — not requested. Local rustc is 1.95 against
CI's 1.99, so the clippy pass is necessary but not sufficient.

### Post-review copy pass

Dual code review (`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`) returned **no Critical or High findings**;
both recommended the change as-is. Four changes were applied from their Low findings, then the full
gated set was re-run (844 pass, clippy and fmt clean, changelogs byte-identical):

| Change | Source |
|---|---|
| Comment reworded to "Post-adapter length, before RRBS/poly-X/N-trim/clipping" | **Both** reviewers: "pre-RRBS" understated a quantity that precedes four further stages, and clipping is what the tests discriminate on |
| CHANGELOG tightened to one sentence, in both copies | Reviewer A, against the project's one-sentence rule |
| `assert_eq!(serial.2, parallel.2)` added to the paired test | Reviewer A: the model test treats `PairValidationStats` parity as load-bearing. Added and passes — the fixture has no short reads, so it does not trip the unpaired divergence below |
| Final diff | 224 insertions, 5 deletions: `src/parallel.rs` +208/−5, both changelogs +8 |

### Follow-ups found during review — recorded here, not filed

Neither is a defect in this change; both are pre-existing and out of its scope.

1. **`r1_unpaired` / `r2_unpaired` diverge with `--cores`** (Reviewer B, reproduced independently).
   `--paired --clip_R1 5 --clip_R2 5` with no `--retain_unpaired`: `--cores 1` omits the "Unpaired
   read" lines entirely, while `--cores 4` prints `Unpaired read 1 kept: 0` and **`Unpaired read 2
   kept: 4`** with no unpaired file written, and `"r2_unpaired": 4` in the JSON. Same class as #541,
   five lines below the edited ones. Already documented as known at `src/parallel.rs:2266-2276`,
   whose comment states the parallel path increments on the filter's `r1_ok`/`r2_ok` regardless of
   whether the rescue writer is open while the serial path increments only when it is open — and
   that `gen_multiome_rows` keeps **both** mates short specifically so neither path increments and
   the parity test passes. That is the same vacuous-fixture pattern this plan was written to remove,
   deliberately left in place. **Worth an issue; Felix's call.**
2. **Coverage gaps adjacent to this field** (Reviewer A). No test in `src/parallel.rs` sets
   `rrbs: true`, so #541's own `--rrbs` row is pinned only by a manual run and
   `rrbs_r2_clipped_5prime` has no parity coverage at all. `--discard_untrimmed` has no behavioural
   coverage anywhere in the repo — which matters here because moving the accumulation above its
   early return would shift **both** paths equally and leave every assertion in the file green.

## 13. Revision History

**Revision 2 — 2026-11-10.** Dual plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`). Both
reviewers verified the claims in source; both found revision 1 factually wrong in three places. All
corrections below were re-verified against the source before being applied.

| Change | Why |
|---|---|
| §2: paired parity test **does** exist (`:2271`, asserting R1+R2 at `:2338-2339`) | Revision 1 said none existed. Both reviewers found it. Step 5 is now an extension, not new work. |
| §2: serial sites split into FASTQ and **uBAM** path families | Revision 1's single "six serial sites" table hid that `:672`/`:780`/`:781` are a separate, already-correct uBAM path. |
| §2: root cause rewritten as an **orphaned call site** | `82d1e34` wrote the parallel sites correctly for their time; `f4b544d` created `bp_after_cutadapt`, switched the serial sites, and never touched `parallel.rs`. Revision 1's "dropped comment clause" story was backwards. |
| §2: `--clip_R1` arithmetic corrected to `5 × 9,996`, four reads unclipped | Revision 1 blamed filtered reads, but accumulation precedes filtering. Right number, wrong mechanism — inherited verbatim into #541, which still needs correcting. |
| §2, §3: `--clumpify` requires `--cores >= 2`, so every clumpify run is affected | Unmentioned in revision 1; likely the largest affected population. |
| **Step 3: helper now sets `clip_r2` and `three_prime_clip_r2`** | **Revision 1 set only R1 in the helper and put R2 in Step 5's prose, contradicted by the adjacent code block. As written, `:543` was undiscriminated and the paired test would pass however `:543` was written — the "test that can only ever pass" failure mode, inside a plan written to fix it.** |
| Step 6 added: absolute-value anchor | Every assertion was `serial == parallel`, pinning the paths to each other but not to anything true. A future flip of the *serial* line would keep the suite green and silently revert #232. |
| Step 7 added: fixture discrimination guard, following `:1518-1524` | Nothing asserted the config still discriminates, so a later fixture edit would silently restore vacuity. |
| Step 8 added: CHANGELOG entry | Both reviewers; `### Unreleased` is empty. |
| §7: #529 **inherits** this fix rather than needing its own | Revision 1 said it adds its own accumulation. It routes uBAM through `classify_*`. |
| Assumption 6 corrected: `clip_*` skips rather than clamps | Revision 1 assumed unconditional application. |
| §9: validations 2 and 7 added; 6 split per site | Revision 1's artifact checks were all single-end, leaving `:542`/`:543` with no end-to-end coverage. |
| §11: removed the claim that the paired test was "the plan's second-biggest correction" | That self-review claim was itself false, resting on the false premise that no paired test existed. |

Steps 1 and 2 — the fix itself — were confirmed correct by both reviewers and are unchanged.
