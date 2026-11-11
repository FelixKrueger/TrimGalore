# Plan Coverage Report

**Mode:** B (no `IMPL.md` — the design plan stands in as the implementation spec, per the
audit brief)
**Plan(s):** `plans/08162026_clumpify-memory-accounting/PLAN.md` (v3)
**Scope:** **Phase 1 only.** Phases 2 and 3 are explicitly future work in the plan's Staging
section and are not counted as gaps.
**Codebase:** `dev`, base commit `2ab0918`, changes uncommitted
**Date:** 2026-08-16
**Verdict:** **INCOMPLETE — 3 items unresolved** (1 MISSING, 2 PARTIAL; all three are
test-side guard strength, none is a behaviour or data-integrity gap)

## Summary

- In-scope items (Phase 1): **20**
- DONE: **15**
- PARTIAL: **2**
- MISSING: **1**
- DEVIATED (both documented with a reason): **2**
- Out of scope (Phases 2–3), listed but not counted: **4**

**Suite state, verified independently:** `cargo test --release --no-fail-fast` →
**729 passed, 0 failed** across 16 test targets. This matches the implementer's recorded
figure exactly. `cargo fmt --all -- --check` clean. `cargo clippy --all-targets --release
-- -D warnings` clean on a **forced** re-check (`Checking trim-galore`, 0 warnings) — the
first clippy invocation returned in 0.33 s, which is a cache hit and not evidence of
anything, so the source files were touched to force a genuine compile.

**One transient failure, discounted with cause.** My first run reported
`integration_output_collision` → 56 passed / 1 failed, in
`hardtrim3_ubam_accepts_distinct_stems_and_names_output`, with
`failed to run trim_galore: Os { code: 2, kind: NotFound }`. That is the test harness
failing to *find* the binary, not the binary misbehaving: two code-review agents were
building in the same working tree and `target/release/trim_galore` was momentarily absent.
The binary is present (13:25, 8,052,288 bytes) and the test passes in the clean re-run.
The first run also stopped at that target (no `--no-fail-fast`), which is why it showed 57
tests rather than 729.

## Coverage ledger — Implementation outline (Phase 1)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | Lever A: work-channel depth 2 → 1, **gated on clumpy mode**, at both construction sites | Outline 1 | DONE | `src/parallel.rs:117` (paired) and `src/parallel.rs:943` (SE): `let queue_depth = if clump_layout.is_some() { 1 } else { 2 };`, consumed at `:121` and `:947`. Grep confirms these are the only two `PairedWork`/`SingleWork` channels; the `result_tx` channels (`:130`, `:959`) are deliberately untouched at `cores * 2`. Plain mode keeps depth 2. |
| 2 | Lever B: reserve output buffers from `estimated_record_bytes` instead of `len * 300`; size R2 from R2 | Outline 2 | **DEVIATED** (documented) | Arithmetic is right: `src/parallel.rs:352-353` size `buf_r1`/`buf_r2` from their **own** mate's sum, and `:1074` sizes the SE buffer from its batch. `grep '\* 300'` over `src/*.rs` returns nothing, so no straggler sites. **The deviation is the gating:** Staging says "Levers A and B, **gated on clumpy mode**"; lever B shipped **ungated**. Documented in the notes ("Built" item 2) with a reason — sizing from the batch's own bytes is correct in every mode. See Gaps. |
| 3 | `ClumpLayout` constructibility — resolve the E0063 hazard | Outline 3 | DONE | `ClumpLayout::new(n_bins, bin_byte_budget, cores)` at `src/clump.rs:148`; the test helper at `src/parallel.rs:1615-1616` now calls it instead of a struct literal. The plan sanctioned either remedy ("updating the struct literal … or adding a constructor"); the constructor was taken. Compiles and the five clumpy tests that use the helper pass. |
| 4 | Drop the drift surface: store the sizing `cores`, remove the parameter from `predicted_peak_bytes` | Outline 4 | DONE | Private `cores` field at `src/clump.rs:142`; `predicted_peak_bytes(&self)` at `:158` reads `self.cores`; `resolve_layout` populates it at `:230`. The only production caller is `src/main.rs:206`, on a layout that `resolve_layout` produced — so the banner cannot be evaluated against a different core count than the sizing used. `src/clump_only.rs` prints no predicted peak, so there is no second drift site. |
| 5 | `--clump_only` fallback decision (warn-and-continue vs documented hard failure) + which `cores` it passes | Outline 5 | **DEVIATED** (documented deferral) | Not implemented. `src/clump_only.rs:275, 413, 748, 905` still call `clump::resolve_layout(memory_budget_bytes, cores.max(1))?` and propagate with a bare `?`. Deferred to Phase 2 in the notes with a stated reason, and **the reason checks out** (see Gaps) — allowed as DEVIATED, not MISSING. Also covers the plan's "`--clump_only`'s abort path and single-threadedness" edge case. |
| 6 | Banner unit label: prints MiB, labelled `MB` | Outline 6 | DONE | `src/main.rs:203` → `"clumpify: {} bins × {} MiB; predicted peak ≈ {} MiB (gzip level {})"`, both divisors `1024 * 1024`. Four more at `src/clump_only.rs:281, 420, 760, 999`. A sweep of `src/` for `MB` finds no remaining unit label on a divide-by-1024² line (the surviving `MB` hits are unrelated batch-size comments in `fastq.rs`, `bam.rs`, `parallel.rs`, `cli.rs`). Live confirmation in the test log: `16 bins × 21 MiB budget`. |

## Coverage ledger — Tests

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 7 | `resolve_layout_default` — affected | Tests | DONE | `src/clump.rs:457`. Needed no edit: Phase 1 changed no sizing arithmetic. Passes. |
| 8 | `resolve_layout_scales_with_cores` — affected | Tests | DONE | `src/clump.rs:467`. Same; passes. |
| 9 | `resolve_layout_bails_when_too_small` — affected | Tests | DONE | `src/clump.rs:495`. Same; passes. |
| 10 | `resolve_layout_bails_below_static_overhead` — affected; **panics** rather than fails if the reservation drops | Tests | DONE | `src/clump.rs:506`, passes unchanged. The fragility the plan named is still there — `res.unwrap_err()` at `:509` with no `is_err()` guard — but it is latent: it only bites when the reservation moves, which is Phase 2. Nothing was owed in Phase 1. |
| 11 | `resolve_layout_predicted_peak_fits_budget` — needs any new argument | Tests | DONE | `src/clump.rs:482` updated to `layout.predicted_peak_bytes()`. Passes. |
| 12 | Same test — "**Retire it or move the assertion somewhere it bites**" | Tests | **MISSING** | Neither was done, and the notes do not mention it. See Gaps. |
| 13 | `"static overhead"` substring assertion that a bail-message rewrite breaks | Tests | DONE | `src/clump.rs:502` intact and passing; the bail message at `:217-218` was not rewritten, so nothing was owed. |
| 14 | `clump_only.rs`'s `small_memory_budget()` helper and its 13 call sites | Tests | DONE | Helper at `src/clump_only.rs:1187` returns 1 GiB, unchanged; 14 occurrences = 1 definition + **13** call sites, matching the plan's count. All pass — consistent with the floor being untouched. |
| 15 | Both `tests/integration_clump_only*.rs` files, default-budget invocations | Tests | DONE | `integration_clump_only.rs` 13 passed, `integration_clump_only_ubam.rs` 16 passed. Their live banners in the log resolve at the default budget (`16 bins × 21 MiB budget`), confirming no invocation flipped to aborting. |
| 16 | Banner guard in `tests/integration_preflight_tripwire.rs` | Tests | **PARTIAL** | New assertion added at `:597-602` and it genuinely bites — it would have failed on the old banner. But it pins only the **bins** label; a regression of the **peak** label alone would slip through. See Gaps. |
| 17 | **New:** banner-present assertion in every measured cell | Tests ("New:") | DONE | In the measurement harness, which is where "measured cell" points: `spikes/spike_memory-levers.py:116-118` raises `SystemExit` unless `"clumpify:"` appears in the run's stderr, per cell. The test suite carries the same guard for the shipped binary at `integration_preflight_tripwire.rs:591-595`. One note on the harness in Gaps. |
| 18 | **New:** regression test pinning resolved `B` at a named cell, so a coefficient edit cannot move bin size unnoticed | Tests ("New:") | **PARTIAL** | No new test was added. The protective effect exists via the two pre-existing range assertions (items 7 and 8), which pin `bin_byte_budget` at two named cells inside ~±2% windows. The requested addition is absent and the omission is undocumented. See Gaps. |

## Coverage ledger — Edge cases explicitly called out (Phase-1 relevant)

| # | Edge case | Source | Status | Notes |
|---|-----------|--------|--------|-------|
| 19 | Single-end's one bin pool — "lever A patches the SE channel too" | Self-Review / Validation 8 | DONE | SE site patched (`src/parallel.rs:943, 947`). Five clumpy tests drive `Some(small_clump_layout())` through the depth-1 queue — `test_clumpy_single_end_is_a_permutation`, `test_clumpy_paired_preserves_pair_lockstep`, `test_clumpy_stats_parity_with_plain`, `test_clumpy_compresses_clusterable_input`, `test_clumpy_tiny_input` — so the shallower queue is shown not to deadlock or drop records, including on tiny input. |
| 20 | Silent fallback masquerading as a good ratio | Self-Review | DONE | Guarded on both sides: harness (item 17) and `se_trim_clumpify`'s `"clumpify did not engage, so this case proves nothing about its writers"` assertion. |

## Out of scope — Phases 2 and 3 (listed for completeness, not counted)

| Item | Plan location | Why out of scope |
|---|---|---|
| `--fastqc_args` alone → use `fastqc_requested()` | Self-Review edge cases | Only arises once the FastQC share of `STATIC_OVERHEAD_BYTES` is split out. Phase 1 left the reservation unconditional at `src/clump.rs:179`, so the case cannot fire. `Cli::fastqc_requested()` exists at `src/cli.rs:592` for Phase 2 to use. |
| FastQC-path floor under mimalloc | Staging Phase 2, Validation 7 | Measurement, not code. |
| Multi-pair per-process peak | Validation 10 | Measurement. |
| Read-length dependence of σ; fixture rebuild on distinct reads; `sort -u` distinct-count logging | Staging Phase 2 | Explicitly Phase 2. The harness logs no distinct counts yet, as expected. |
| Docs arithmetic, CI check, calibrated model | Staging Phases 2–3 | Explicitly future. **Worth flagging for Phase 3:** `docs/src/content/docs/performance/clumpy.md:102` still reads "**32 bins × 12 MB**" and "32 bins × 66 MB", and the same sentence says "The binary prints those resolved values at startup" — which now prints MiB. The docs and the banner disagree until Phase 3 lands. |

## Gaps (detail)

### Item 12: `resolve_layout_predicted_peak_fits_budget` was neither retired nor moved — MISSING

**Expected:** the plan states the test "**passes for any coefficient** — it compares the
formula against its own inverse. Retire it or move the assertion somewhere it bites."

**Found:** the test still exists at `src/clump.rs:477-492`, still asserting
`predicted_peak <= budget` over four `(memory, cores)` cells. The only change is the
signature update from item 11. Because `predicted_peak_bytes()` is the algebraic inverse of
`resolve_layout`'s division, the assertion remains a tautology — item 4 arguably tightened
the tautology, since the layout now carries the very `cores` the prediction divides by.

**Gap:** neither prescribed action was taken, and the implementation notes do not mention
the item at all, so this is a silent skip rather than a documented deviation.

**One nuance for the reader to rule on:** the same plan sentence continues "Phase 3's CI
check is the real replacement", which can be read as pairing the retirement with Phase 3.
The item sits in the Tests section with no phase label, alongside the signature change that
*was* made, so I have scored it in scope. If Felix reads it as Phase 3, this item drops and
the verdict moves to INCOMPLETE — 2 items.

### Item 16: the new banner guard pins only one of the two labels — PARTIAL

**Expected:** a banner guard covering the unit-label fix, which is a fix to **two** labels
in one line.

**Found:** `tests/integration_preflight_tripwire.rs:598-602`:

```rust
out.stderr.contains("MiB; predicted peak ≈") && !out.stderr.contains(" MB;"),
```

The banner is `"clumpify: {} bins × {} MiB; predicted peak ≈ {} MiB (gzip level {})"`. Both
clauses key on the **bins** label: `"MiB; predicted peak ≈"` matches the `MiB` that precedes
the semicolon, and `" MB;"` is the pre-fix form of that same label. The **peak** label is
followed by `" (gzip"`, not `";"`, so a revert of the peak label alone yields
`"… MiB; predicted peak ≈ 884 MB (gzip level 1)"` — first clause true, second clause false —
and the assertion **passes** on a wrong label.

**Gap:** extend the guard to the peak figure, e.g. also require
`"predicted peak ≈ "` to be followed by a `MiB` token and reject `" MB ("`. The claim in the
notes that the assertion "would have failed before the change" is correct as stated; the
shortfall is forward coverage, not the historical check.

### Item 18: no new test pins resolved `B` at a named cell — PARTIAL

**Expected:** "a regression test pinning resolved `B` at a named cell so a future
coefficient edit cannot move bin size unnoticed", listed under "New:".

**Found:** no new test. `resolve_layout_default` (4 GiB / 2 cores → `bin_byte_budget` in
[150, 156) MiB) and `resolve_layout_scales_with_cores` (4 GiB / 8 cores → [65, 70) MiB) do
pin `B` at two named cells with windows tight enough that a coefficient edit would break
them, so the intent is substantially covered by incumbents.

**Gap:** the requested new test is absent and the notes do not record why. Since Phase 2 is
where coefficients move, this is the guard that Phase 2 will lean on — worth resolving
before Phase 1 is called finished, either by adding the test or by recording that items 7
and 8 are deemed to discharge it.

## Deviations (both documented — informational, not counted as gaps)

### Item 2: lever B shipped ungated, where Staging said "gated on clumpy mode"

The plan's Staging reads "**Phase 1 — structural reduction.** Levers A and B, gated on
clumpy mode." Lever B is unconditional: `src/parallel.rs:352-353` and `:1074` apply to every
run. The notes document this explicitly, with a reason ("sizing from the batch's own bytes
is correct in every mode, not a memory lever conditional on clumpy"), so it scores DEVIATED.

Two consequences a reader should see, both already visible in the notes' own numbers:

- **Validation row 6 expects "plain path untouched".** With lever B ungated the plain path
  *is* touched on memory: the notes' table records plain `c4` at 96.8 → 91.2 MiB, −5.8%.
  Output is byte-identical (407,440,179 bytes) and CPU is +0.47%, so the row's intent —
  no throughput or output regression — holds; its literal wording no longer does.
- The deviation is not listed in the notes' "**Deviation from the plan**" section, which
  names only item 5. It is documented inline in the "Built" list instead. Reason given, so
  the substance is covered; the notes' own deviation register is incomplete.

### Item 5: `--clump_only` fallback deferred to Phase 2

The notes defer it with this reason: "The abort path only becomes reachable when the floor
moves, and Phase 1 changes no sizing arithmetic — `clumpify_min_memory_bytes` is untouched,
so every budget that resolved before still resolves."

**I verified the reason rather than taking it.** In the diff, `src/clump.rs` gains only the
`cores` field, `ClumpLayout::new`, and the signature change. `resolve_layout`'s arithmetic,
`STATIC_OVERHEAD_BYTES` (512 MiB, `:179`), `MIN_BIN_BYTES` (1 MiB, `:38`) and
`clumpify_min_memory_bytes` (`:168-172`) are all byte-for-byte unchanged. Empirically, the
29 `clump_only` integration tests pass and their banners resolve at the default budget
(`16 bins × 21 MiB`), well clear of the 1 MiB floor. The deferral is therefore safe for
Phase 1 and correctly reasoned, and the open behaviour question it names (warn-and-continue
vs documented hard failure) is genuinely Felix's call. Accepted as DEVIATED.

## Test verification

| Test / target | File | Status |
|---|---|---|
| Whole suite, 16 targets | — | **729 PASS, 0 FAIL** (`--no-fail-fast`) |
| `resolve_layout_default` | `src/clump.rs:457` | PASS |
| `resolve_layout_scales_with_cores` | `src/clump.rs:467` | PASS |
| `resolve_layout_predicted_peak_fits_budget` | `src/clump.rs:477` | PASS (tautological — item 12) |
| `resolve_layout_bails_when_too_small` | `src/clump.rs:495` | PASS |
| `resolve_layout_bails_below_static_overhead` | `src/clump.rs:506` | PASS |
| `se_trim_clumpify` (+ new MiB assertion) | `tests/integration_preflight_tripwire.rs:581` | PASS (guard partial — item 16) |
| `test_clumpy_single_end_is_a_permutation` | `src/parallel.rs:1637` | PASS |
| `test_clumpy_paired_preserves_pair_lockstep` | `src/parallel.rs:1691` | PASS |
| `test_clumpy_stats_parity_with_plain` | `src/parallel.rs:1794` | PASS |
| `test_clumpy_compresses_clusterable_input` | `src/parallel.rs:1840` | PASS |
| `test_clumpy_tiny_input` | `src/parallel.rs:1914` | PASS |
| `integration_clump_only` (13 tests) | `tests/integration_clump_only.rs` | PASS |
| `integration_clump_only_ubam` (16 tests) | `tests/integration_clump_only_ubam.rs` | PASS |
| Regression test pinning resolved `B` | — | **MISSING as a new test** (item 18) |
| `hardtrim3_ubam_accepts_distinct_stems_and_names_output` | `tests/integration_output_collision.rs` | PASS on clean re-run; transient NotFound in a concurrent-build race |

## Other observations (no plan item owed)

- **`spikes/spike_memory-levers.py:119` will now fail to parse the shipped banner.** It does
  `int(text.split("predicted peak ≈ ")[1].split(" MB")[0])`. `" MB"` is not a substring of
  `" MiB"`, so the split no longer separates the number and the `int()` raises. Iteration
  log #3 records fixing "the verification harness" regex to accept both forms; the **spike**
  harness in the plan directory was not updated. It is a spike artifact, so nothing was owed
  in Phase 1 — but Phase 2 measures with this family of scripts, so it is a live trap.
- **`ClumpLayout::new` reopens, in principle, the drift the plan closed.** Item 4 removed the
  ability for a caller to evaluate the prediction against a different core count; the new
  constructor lets a caller supply any `cores`. Production never does — `main.rs` sizes and
  predicts from one `resolve_layout` result. The one live example is the test helper
  (`src/parallel.rs:1616`), which declares `cores = 2` and is then handed to runs at
  `cores = 4`; harmless because `predicted_peak_bytes()` is never called on it. The
  constructor was one of the two remedies the plan itself sanctioned, so this is a note, not
  a gap.
- **The `CHANGELOG.md` entry is beyond the plan**, which specified no changelog work. Its
  presence is not a gap in either direction.

## Verdict

**INCOMPLETE — 3 items unresolved.** All three are test-side; no behaviour item in the
Phase 1 outline is missing, output byte-identity is unaffected, and the suite is green at
729/0.

To close Phase 1:

1. **Item 12 (MISSING)** — retire `resolve_layout_predicted_peak_fits_budget`, or move its
   assertion somewhere it is not the formula's own inverse. If this is intended as Phase 3
   work alongside the CI check, record that and the item drops.
2. **Item 16 (PARTIAL)** — extend the `se_trim_clumpify` banner guard to the *peak* label; as
   written it only pins the *bins* label, and a revert of `predicted peak ≈ {} MB` alone
   passes.
3. **Item 18 (PARTIAL)** — add the requested regression test pinning resolved `B` at a named
   cell, or record that `resolve_layout_default` and `resolve_layout_scales_with_cores`
   discharge it.

Two further decisions for Felix, not defects: whether lever B may stay **ungated** (it makes
Validation row 6's "plain path untouched" literally false while keeping its intent), and
whether item 5's deferral of `--clump_only`'s fallback to Phase 2 stands (its stated reason
verifies).
