# Plan Coverage Report

**Mode:** B (code vs. plan §4/§3/§9)
**Plan(s):** plans/a2-ignored-with-presets/PLAN.md (v2)
**Date:** 2026-07-29
**Verdict:** INCOMPLETE — 3 items unresolved (one of them a live `cargo fmt` failure introduced during this audit)

## Summary

- Total items: 53
- DONE: 47 / PARTIAL: 2 / MISSING: 1 / DEVIATED: 3

Every behaviour requirement in §3 is satisfied and 24 of the 25 §9 validation rows pass, verified independently rather than from §13's table. Unresolved: §4 Step 9 (the #369 comment, MISSING — self-reported as deferred), §4 Step 5 (test coverage, PARTIAL — three named sub-requirements absent), and **V14, which now fails**: `cargo fmt --all -- --check` exits 1 on the current working tree.

### The tree changed during this audit

`src/adapter.rs` was **not** modified when I started (`git status` showed five modified tracked files: `ci.yml`, `CHANGELOG.md`, `flags.md`, `cli.rs`, `main.rs`). At 16:03:42, mid-audit, a sixth appeared: `src/adapter.rs`, adding `parse_adapter_spec_quiet` / `parse_adapter_specs_quiet` and rewiring `cli.rs:941` to the silent variant. I did not make this edit and did not revert it.

Its effect is to resolve what I had recorded as Observation O2 — `Cli::validate` no longer re-announces an `A{N}` expansion that the trimming run announces itself. Re-verified after the change: the notice now prints **once** (was twice), the release binary is newer than the sources so it includes the change, the full suite is still **454 passed / 0 failed**, and `cargo clippy --all-targets --release -- -D warnings` still exits 0.

But it also **breaks `cargo fmt`** — see the V14 gap below. Everything else in this report was measured before 16:03:42; because the change is stderr-only, none of the output-identity findings (V4–V7, V13, V16) is affected.

### Verification independence

§13's own results table was treated as a claim, not evidence. I rebuilt the **pre-fix** binary from `dev` @ `4276810` in a throwaway worktree with a separate `--target-dir`, which let me generate the "before" baselines myself instead of trusting captured ones. That converts V4, V5, V6, V7b, V7c, V13 and V15 from unverifiable-in-principle into directly measured. The worktree was removed afterwards; `git status` shows exactly the five expected modified tracked files and no repo file was altered except this report.

## Coverage ledger

### §4 Implementation outline

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `apply_adapter2_override` added, called from `setup_trimming`; `resolve_adapter` **not** restructured; redundant `-a2` parse deleted from the `-a` branch | §4.1 | DONE | `src/main.rs:821` (call), `:994` (helper). All seven branches of `resolve_adapter` (`:1026`–`:1088`) verified byte-unchanged apart from the `-a` branch's r2 slot becoming `Vec::new()` (`:1029`) — equivalent, since `parse_adapter_specs(&[])` returns an empty vec. Read 1 guarantee is structural, not merely tested |
| 2 | `-a2` validated in `Cli::validate` so malformed input fails before the 1 M-read scan | §4.2 | DONE | `src/cli.rs:940` calls `parse_adapter_specs(&self.adapter2)?` ahead of any I/O |
| 3 | Override NOTE emitted after the adapter display; `Adapter 2 (Read 2):` gated on `cli.paired` | §4.3 | DONE | `main.rs:847` (NOTE), `:836` (`cli.paired &&` gate). Measured: SE `--small_rna` prints 1 `Adapter 2` line pre-fix, 0 post-fix |
| 4 | Unusable-`-a2` warning for single-end **and** all four specialty modes with a per-mode reason; suppression carve-out + warning | §4.4 | DONE | `cli.rs:941-956` (four-arm reason chain), `main.rs:1004-1011` (carve-out). All five modes verified to warn with a mode-true reason |
| 5 | New `tests/integration_adapter2.rs`: per-branch e2e for all seven branches, plus V16–V23; **and** a unit test of the helper | §4.5 | **PARTIAL** | 13 tests, all passing; seven branches and V16–V22 covered. Missing: (a) no unit test of `apply_adapter2_override` anywhere, (b) V23 (uBAM) not automated, (c) `--implicon` and `--hardtrim3` warning branches not automated. See Gaps |
| 6 | CI `validation` matrix case: PE `--illumina -a2 <seq>` md5-compared against Perl 0.6.11 | §4.6 | DONE | `.github/workflows/ci.yml:334-356`, inside the `validation` job. Covers `--illumina` **and** auto-detect, md5s both mates, plus a grep that R2 got the `-a2` sequence |
| 7 | CHANGELOG appended to the existing `#### Fixes` block under `### Unreleased` | §4.7 | DONE | `CHANGELOG.md:86-135`. States all-v2.x scope, the `-a`-only exception, and that Read 1 *per-read trimming* is unchanged while `_val_1.fq.gz` differs — the §1 correction is carried through correctly |
| 8 | Docs: `flags.md:12` hedge; the five other cited spots become correct as written | §4.8 | DONE | Hedge applied at `flags.md:12`, extended to `--bgiseq` (documented deviation 3). All six cited examples executed and verified to now do what they claim |
| 9 | Comment on #369 linking the PR; close on merge | §4.9 | **MISSING** | #369 is still `OPEN` with 2 comments (`FelixKrueger`, `MathieuUm`) — none links this fix. §13 self-reports this as deferred pending go-ahead. Not actioned by me, per instruction |

### §3 Behaviour

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 10 | Resolution order: candidate from `resolve_adapter`, override applied in the caller; `-a2` parsed early; `parse_adapter_specs` on every path | §3.1 | DONE | Single decision point confirmed — `resolve_adapter` has exactly one caller. `A{N}`, repeated `-a2` and `file:` all verified working on preset **and** auto-detect paths |
| 11 | Truth table, all four rows | §3.2 | DONE | Row 1 (`--illumina`/`--nextera`/`--stranded_illumina`) → `-a2`. Row 2 → `-a2` + NOTE, verified for `--small_rna`, `--bgiseq` **and** auto-detected BGI. Row 3 (`--illumina`, no `-a2`) → R1 fallback. Row 4 (`--small_rna`/`--bgiseq`, no `-a2`) → default stands. Read 1 unaffected in all four |
| 12 | Override NOTE, only on real displacement | §3.3 | DONE | `NOTE: Read 2 adapter taken from -a2; the smallRNA default (GATCGTCGGACT) is not used.` Absent for `--illumina`/`--nextera`/`--stranded_illumina`/auto-detected-Illumina; present for `--small_rna`, `--bgiseq`, auto-detected BGI (label `BGI/DNBSEQ`, correct) |
| 13 | `-a2` where it cannot apply warns, with a reason true in every firing case | §3.4 | DONE | Reason-neutral prefix + specific cause, matching the plan's prescribed wording. Sited in `Cli::validate` before the deprecation block, as specified. `--clump_only` left alone (still hard-errors) |
| 14 | Unchanged behaviour, all six bullets | §3.5 | DONE | R1 selection/precedence/trimming; `-a`+preset permissive with no new warning; `--small_rna` cutoff 18; `-a2` absent → default intact; byte-identity on every non-`-a2` path (measured pre vs post across 25 invocations) |
| 15 | Newly-fatal malformed `-a2`; `--consider_already_trimmed` suppression wins | §3.6 | DONE | `ZZZQQQ`, `""`, `file:missing.fa` all exit 1 with **the identical message** the `-a` path gives. Suppression: warning fires, both mates `Reads with adapters: 0 (0.0%)`, no `Adapter 2` line, no NOTE |

### §9 Validation (V1–V23, verified independently)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 16 | V1 — reported case fixed | §9 | DONE | `--paired --illumina -a2 AAATCAAAAAAAC` → R2 report `-a AAATCAAAAAAAC` |
| 17 | V2 — all five presets honour `-a2` | §9 | DONE | `--nextera`, `--small_rna`, `--stranded_illumina`, `--bgiseq` all report `-a AAATCAAAAAAAC` for R2 |
| 18 | V3 — auto-detect honours `-a2` | §9 | DONE | No adapter flag → R2 gets `-a2`; R1 per-read figures unchanged from a no-`-a2` run |
| 19 | V4 — `-a` + `-a2` unchanged | §9 | DONE | **Pre vs post binary: byte-identical** across FASTQ, `.txt` and `.json` (paths normalised) |
| 20 | V5 — `-a2` absent → preset default intact | §9 | DONE | `--small_rna` → `GATCGTCGGACT`, `--bgiseq` → BGI R2. Both byte-identical pre vs post |
| 21 | V6 — `-a2` absent, no R2 default → R1 fallback | §9 | DONE | `--illumina` → R2 uses `AGATCGGAAGAGC`; byte-identical pre vs post |
| 22 | V7a — R1 per-read trimming never moves | §9 | DONE | `Reads with adapters: 3,775 (37.8%)` / `Reads written: 10,000 (100.0%)` identical pre vs post, and identical between `-a2` and no-`-a2` runs on both `--illumina` and auto-detect |
| 23 | V7b — R1 **file** identity where `-a2` does not newly apply | §9 | DONE | Verified for V4, V5, V6 and all V13 cases via the pre-fix binary. Correctly **not** asserted for V1–V3 |
| 24 | V7c — reports too, not just FASTQ | §9 | DONE | `*_trimming_report.txt` and `.json` compared alongside every `.fq.gz`; identical |
| 25 | V8 — `--small_rna` cutoff unaffected by `-a2` | §9 | DONE | `length cutoff of 18 bp` in the report with **and** without `-a2`; A4 holds |
| 26 | V9 — NOTE only on real displacement | §9 | DONE | Present for `--small_rna`; absent for `--illumina`, `--nextera`, `--stranded_illumina`, auto-detected Illumina |
| 27 | V10 — SE `-a2` warns, output unchanged | §9 | DONE | Warning on stderr, exit 0, trimmed FASTQ byte-identical to omitting `-a2` |
| 28 | V11 — `-a` + preset still permissive | §9 | DONE | `--illumina -a AAAACCCCGGGG` → exit 0, `-a` wins, **zero** warnings or errors (A6 preserved) |
| 29 | V12 — `A{N}` / `file:` / repeated `-a2` on preset paths | §9 | DONE | `A{10}` → `AAAAAAAAAA`; `file:` → both FASTA entries; two `-a2` → both sequences. Also verified with **no** `-a` at all |
| 30 | V13 — validation matrix unmoved | §9 | DONE (surrogate) | Cannot run the job itself (needs a Linux runner + conda Perl 0.6.11). Substituted a **stronger local check**: every one of the job's `trim_galore` invocations run under both the pre-fix and post-fix binary — 25 cases, **all identical**. See Gaps for what this does and does not establish |
| 31 | V14 — full suite + fmt + clippy | §9 | **PARTIAL** | `cargo test --release`: **454 passed, 0 failed** (matches §13, and still 454/0 after the mid-audit change). `cargo clippy --all-targets --release -- -D warnings` exit 0, forced re-check. But **`cargo fmt --all -- --check` exits 1** on the current tree — regression introduced at 16:03:42 by the `src/adapter.rs` edit. See Gaps |
| 32 | V15 — the new CI gate fails pre-fix | §9 | DONE | **Reproduced.** Pre-fix binary reports `-a AGATCGGAAGAGC` for R2 under both `--illumina` and auto-detect, so the step's `grep -q 'a AAATCAAAAAAAC'` fails. The gate is non-vacuous |
| 33 | V16 — differential oracle | §9 | DONE | Post-fix `--illumina -a2 SEQ` **and** post-fix auto-detect `-a2 SEQ` are byte-identical, both mates, to the **pre-fix** `-a AGATCGGAAGAGC -a2 SEQ` output. Also asserted post-vs-post by `oracle_equivalence` |
| 34 | V17 — malformed `-a2` fatal on preset paths | §9 | DONE | All three forms exit 1. Messages verified **character-identical** to the `-a` path, as V17 requires |
| 35 | V18 — suppression wins over `-a2` | §9 | DONE | Warning fires; R1 and R2 both `Reads with adapters: 0 (0.0%)` (symmetric); no `Adapter 2` line; no NOTE |
| 36 | V19 — specialty modes warn with a mode-true reason | §9 | DONE | All five verified, arguments passed individually to avoid the zsh word-splitting trap: `--clock`, `--implicon`, `--hardtrim5`, `--hardtrim3`, single-end. Each names its own reason; all exit 0 |
| 37 | V20 — SE output not self-contradictory | §9 | DONE | `--illumina -a2 SEQ single.fq.gz`: warning present, zero `Adapter 2 (Read 2)` lines. See Observation O1 for a residual case outside V20's specified invocation |
| 38 | V21 — poly-G piggyback survives | §9 | DONE | `Scanning for poly-G content` = **1** for `--illumina`, **0** for auto-detect; identical pre vs post. Poly-G verdict line identical in both |
| 39 | V22 — repeated `-a2` sizes R2 stats correctly | §9 | DONE | R2 report contains **2** `=== Adapter N ===` blocks for two `-a2` values, and 3 for three. `r2_adapter_count` switches correctly |
| 40 | V23 — uBAM output inherits the fix | §9 | DONE | `--paired --illumina -a2 SEQ --output-format ubam` → one interleaved `_val.bam`, R2 report `-a AAATCAAAAAAAC`, `Adapter 2 (Read 2)` line present. Behaviour correct; **no automated test** (folded into Step 5's PARTIAL) |

### §11 Enumerated edge cases

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 41 | `A{N}` shorthand for `-a2` on preset paths | §11 | DONE | Also with no `-a` (`-a2 'T{12}'` → 12 T's) |
| 42 | `file:` spec for `-a2` on preset paths | §11 | DONE | Both with and without a paired `-a file:` |
| 43 | Repeated `-a2` on preset paths | §11 | DONE | Two and three sequences both reach R2 |
| 44 | `-a2` empty after parsing | §11 | DONE | `-a2 ""` → `Error: Empty adapter sequence`, exit 1 |
| 45 | `-a2` on single-end | §11 | DONE | Warns, exits 0, output unchanged, no `Adapter 2` line |
| 46 | `-a2` with `--clump_only` still rejected | §11 | DONE | `Error: --clump_only does not trim; -a2/--adapter2 is not compatible`, exit 1 — unchanged |
| 47 | Preset **with** an R2 default and no `-a2` | §11 | DONE | V5; byte-identical pre vs post |
| 48 | Preset **without** an R2 default | §11 | DONE | V6; byte-identical pre vs post |
| 49 | `-a` + preset + `-a2` (the #369 workaround) | §11 | DONE | V4/V11; byte-identical pre vs post, no new warning |
| 50 | Auto-detected BGI (R2 default via detection, not a flag) | §11 | DONE | Verified on a purpose-built fixture: detection picks `BGI/DNBSEQ`, `-a2` displaces its R2 default, NOTE names the right preset and sequence. Without `-a2` the detected default stands. **Nothing in the test suite reaches this path** — the standard fixture always detects Illumina |

### §13 self-reported deviations (verified as claims)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 51 | Suppression detected from the R1 adapter shape, not `DetectionResult.suppressed` | §13 dev. 1 | DEVIATED (documented) | **Claim verified sound.** `seq: ""` occurs exactly once in `src/adapter.rs` (`:222`, the synthetic suppressed preset); `parse_adapter_spec` rejects every empty form (`:314-330`, `:328-330`, `:335-340`); all real presets have non-empty `seq`. So `len()==1 && seq.is_empty()` is a faithful discriminator, not a proxy. Also confirmed suppression is reachable **only** on the auto-detect path — clap rejects preset + `--consider_already_trimmed` (exit 2, identical pre and post), so the carve-out's placement is complete |
| 52 | No `md5` dev-dependency; decompressed bytes compared directly | §13 dev. 2 | DEVIATED (documented) | Equivalent strength for V16's purpose, and consistent with the project's no-external-deps posture. `flate2::read::MultiGzDecoder` is already in the tree |
| 53 | One extra docs edit (`--bgiseq` hedge alongside `--small_rna`) | §13 dev. 3 | DEVIATED (documented) | Present at `docs/src/content/docs/guide/flags.md:12`. Justified — `--bgiseq` has the identical auto-set behaviour |

## Gaps (detail)

### Item 9: §4 Step 9 — comment on #369 and close on merge

**Expected:** A comment on issue #369 linking the fix/PR, and a manual close on merge (auto-close does not fire on `dev`).
**Found:** #369 is still `OPEN`. It has 2 comments, by `FelixKrueger` and `MathieuUm` — neither links this work. No PR exists yet; nothing is committed or pushed.
**Gap:** Post the comment and close the issue after merge. §13 self-reports this accurately as "Not done — awaiting go-ahead", and it is outward-facing, so deferral is a defensible sequencing choice rather than an oversight. I did not comment or close anything, per instruction.

### Item 5: §4 Step 5 — test coverage has three named sub-requirements unmet

**Expected:** `tests/integration_adapter2.rs` covering all seven branches "plus V16–V23", **and** "a unit test of the helper as well, but not as a substitute".
**Found:** 13 tests, all passing, covering the seven branches and V16–V22.
**Gap:** three specifics.

1. **No unit test of `apply_adapter2_override`.** `grep -rn apply_adapter2_override src tests` returns only the definition (`main.rs:994`), the call site (`:821`) and a doc-comment mention; `src/main.rs` has no `#[cfg(test)]` module. The plan's own parenthetical sanctions the integration-binary route ("the helper can be tested via the integration binary or moved to the library if a unit test is wanted"), so this is close to plan-permitted — but the sentence "Keep a unit test of the helper as well" was an instruction, and §13 does not list its omission among the deviations. That undisclosed omission is why this is PARTIAL rather than DONE.
2. **V23 is not automated.** No test file combines `-a2` with `--output-format ubam` (`grep -rn "output-format" tests | grep -c a2` → 0). The behaviour is correct — I verified it directly — but §13's table reports "V23 Pass" without recording that the evidence is manual only, so nothing guards it against regression.
3. **Two warning branches are not automated.** `unusable_a2_warns_with_a_mode_accurate_reason` covers `clock`, `hardtrim5` and `single_end`. `--implicon` has its own conditional arm (`cli.rs:947`) with no automated coverage; `--hardtrim3` shares `hardtrim5`'s arm, so its risk is lower. The plan's V19 named `--clock`, `--hardtrim5` **and** `--hardtrim3`. All five modes verified correct manually.

**Severity:** low. No behaviour is wrong; the shortfall is regression protection, and the primary bug class (a branch not calling the override) is fully covered by `a2_wins_on_every_resolution_path`.

### Item 31: V14 — `cargo fmt --all -- --check` currently fails

**Expected:** "`cargo test`; `cargo fmt --all -- --check`; `cargo clippy --all-targets --release -- -D warnings` → All pass, zero warnings."
**Found:** tests pass (454/0) and clippy exits 0, but `cargo fmt --all -- --check` **exits 1**. The offending hunk is the new function signature in `src/adapter.rs:375`, which rustfmt wants collapsed onto one line:

```
-fn parse_adapter_specs_inner(
-    specs: &[String],
-    announce: bool,
-) -> Result<Vec<(String, String)>> {
+fn parse_adapter_specs_inner(specs: &[String], announce: bool) -> Result<Vec<(String, String)>> {
```

**Gap:** run `cargo fmt --all`. This is CI-blocking, not cosmetic — the `lint` job (`ci.yml:153`) runs `cargo fmt --all -- --check` verbatim, and `CLAUDE.md` records "CI fails on unformatted code". It was introduced at 16:03:42 by the mid-audit `src/adapter.rs` edit, after §13's validation table was written; §13's "fmt --check clean" claim was true of the tree it described. One command fixes it.

### Item 30: V13 — what the surrogate does and does not establish

**Expected:** Run the whole `validation` job (~30 steps, 24 byte-identity assertions) on a branch push.
**Found:** Not runnable locally — the job needs a Linux runner plus conda-installed Perl Trim Galore 0.6.11. I substituted a differential run: **every** `trim_galore` invocation in the job (extracted from `ci.yml:243-832`), executed under both the pre-fix and post-fix binary, comparing all outputs — 25 cases including SE, PE, multi-pair, `--cores 1`, `-j 4`, `--retain_unpaired`, hardtrim5/3, clock, implicon, demux, `--fastqc`, all four `--clump_only` variants, uBAM output, all five presets, and `--consider_already_trimmed`.

**Result: all 25 identical.** Two apparent diffs resolved to artefacts of my own harness: reports and the uBAM `@PG CL:` line embed the invoking binary's absolute path, and my pre-fix binary sat at a longer path (which also shifted the BAM byte count, 723 vs 672). After normalising the binary path and output directory, uBAM headers are identical and record multisets match under `samtools view`.

**Can the diff affect a non-`-a2` invocation?** Structurally, no: `apply_adapter2_override` returns `Ok(None)` immediately when `cli.adapter2.is_empty()`, and the `Cli::validate` block is inside `if !self.adapter2.is_empty()`. The `-a` branch's `Vec::new()` is equivalent to `parse_adapter_specs(&[])`. The **one** non-`-a2` change is stderr-only: SE runs with an R2 candidate (`--small_rna`, `--bgiseq`, auto-detected BGI) no longer print `Adapter 2 (Read 2):` — measured, 1 line → 0. `grep -n "Adapter 2" .github/workflows/ci.yml` returns nothing, so no assertion depends on it.

**Not confirmed:** that the outputs still match **Perl 0.6.11**. My check proves the fix moved nothing; it cannot prove the pre-existing baselines were correct — but they were green on `dev`, so this is the right invariant. The new step at `ci.yml:334` establishes a genuinely new Perl baseline and must be watched on its first run. A5 re-verified independently: the anchored grep now returns 7 hits, all inside the new step.

## Observations (not plan-item gaps)

Recorded because they are user-visible and were found while auditing, not because the plan requires otherwise.

**O1 — single-end `--small_rna`/`--bgiseq` with `-a2` emits contradictory stderr.**

```
WARNING: -a2/--adapter2 was given but is not used in this mode (it applies to Read 2 of a pair, and this is a single-end run). Ignoring.
Adapter: smallRNA (TGGAATTCTCGG)
NOTE: Read 2 adapter taken from -a2; the smallRNA default (GATCGTCGGACT) is not used.
```

The warning says the flag is ignored; the NOTE says the Read 2 adapter was taken from it. Both are printed. Trimming is unaffected (single-end never consults an R2 adapter). This is **not** a plan-item gap: §3.3 does not gate the NOTE on `cli.paired`, and V20 specifies only `--illumina -a2 SEQ single.fq.gz`, which has no R2 default and therefore no NOTE — that case passes. Flagged because V20's stated purpose was "SE output is not self-contradictory", and this is the residue of exactly that class. Step 3 gated the `Adapter 2 (Read 2):` line on `cli.paired` but not the NOTE.

**O2 — `A{N}` expansion announced twice — RESOLVED mid-audit.** `--illumina -a2 'A{10}'` originally printed `Adapter sequence A{10} expanded to AAAAAAAAAA` twice, because `parse_adapter_specs` ran once in `Cli::validate` and again in the override helper. The plan prescribed exactly that double parse (§4 Step 2: "then parse once in the override helper"), so it was a consequence of the specified design rather than a deviation; §6's cost analysis covered the second parse but not its stderr side effect. The 16:03:42 `src/adapter.rs` edit added silent `_quiet` parse variants and pointed `cli.rs:941` at one. Re-measured after the change: the notice prints **once**. The same edit is what broke `cargo fmt` (Item 31).

**O3 — the NOTE repeats per pair in multi-pair runs.** Two pairs with `--small_rna -a2 SEQ` produce two NOTEs, since `setup_trimming` runs per pair. Consistent with per-pair adapter resolution; the plan states nothing either way.

**O4 — housekeeping.** A worktree from the implementation session is still registered: `/private/tmp/claude-501/tg_dev_369` (detached at `4276810`). Left in place — not mine to remove. My own audit worktree was removed.

## Test verification

| Test name | File | Status |
|-----------|------|--------|
| `a2_wins_on_every_resolution_path` | tests/integration_adapter2.rs | PASS (7 branches) |
| `oracle_equivalence` | tests/integration_adapter2.rs | PASS |
| `read1_trimming_is_unaffected_by_the_r2_adapter` | tests/integration_adapter2.rs | PASS |
| `preset_r2_defaults_survive_when_no_a2_given` | tests/integration_adapter2.rs | PASS |
| `r2_falls_back_to_r1_when_no_default_and_no_a2` | tests/integration_adapter2.rs | PASS |
| `note_fires_only_on_a_real_displacement` | tests/integration_adapter2.rs | PASS |
| `suppression_wins_over_a2` | tests/integration_adapter2.rs | PASS |
| `unusable_a2_warns_with_a_mode_accurate_reason` | tests/integration_adapter2.rs | PASS (3 of 5 modes) |
| `single_end_does_not_announce_a_read2_adapter` | tests/integration_adapter2.rs | PASS |
| `malformed_a2_is_rejected_on_preset_paths` | tests/integration_adapter2.rs | PASS (3 forms) |
| `repeated_a2_applies_on_a_preset_path` | tests/integration_adapter2.rs | PASS |
| `brace_expansion_works_for_a2_on_a_preset_path` | tests/integration_adapter2.rs | PASS |
| `polyg_scan_is_not_duplicated_or_lost` | tests/integration_adapter2.rs | PASS |
| unit test of `apply_adapter2_override` | — | **MISSING** (§4 Step 5) |
| V23 uBAM + `-a2` regression test | — | **MISSING** (§4 Step 5) |
| `--implicon` unusable-`-a2` warning | — | **MISSING** (§4 Step 5) |
| Full workspace suite | all | **454 passed, 0 failed** (re-run after the mid-audit change: still 454/0) |
| `cargo fmt --all -- --check` | — | **FAIL, exit 1** — `src/adapter.rs:375`, introduced 16:03:42 |
| `cargo clippy --all-targets --release -- -D warnings` | — | exit 0, forced re-check (re-confirmed after the change) |

## Verdict

**INCOMPLETE — 3 items unresolved.** The fix itself is complete and correct: all six §3 behaviour subsections hold, all 10 §11 edge cases are covered, 24 of 25 §9 rows pass under independent verification, and all three §13 deviations are both documented and sound — deviation 1's discriminator in particular is a faithful reading of the suppressed-preset semantics, not a fragile proxy. The Read 1 guarantee is structural: zero branches of `resolve_adapter` were touched.

To close out, in priority order:

1. **V14 (PARTIAL) — CI-blocking, one command.** `cargo fmt --all -- --check` exits 1 on `src/adapter.rs:375`. Run `cargo fmt --all`. Introduced at 16:03:42 by an edit made during this audit, so it postdates §13's (then-accurate) "fmt --check clean" claim. The `lint` job runs this command verbatim and will fail the build.
2. **§4 Step 9 (MISSING)** — comment on [#369](https://github.com/FelixKrueger/TrimGalore/issues/369) linking the PR, and close the issue manually after merge (auto-close does not fire on `dev`). Deliberately deferred, and outward-facing; needs a go-ahead, not a fix.
3. **§4 Step 5 (PARTIAL)** — three sub-requirements absent: no unit test of `apply_adapter2_override`; V23 (uBAM + `-a2`) verified manually but not automated; `--implicon`'s warning arm not automated. No behaviour is wrong — this is regression protection only. The first is close to plan-permitted by §4 Step 5's own parenthetical; it counts against coverage because §13 lists three deviations and does not disclose it.

Only item 1 blocks merge, and it is trivial. Two further things worth a maintainer's eye, neither a plan gap: **O1**, the contradictory single-end stderr for `--small_rna`/`--bgiseq` with `-a2` (a one-line `cli.paired &&` on the NOTE would settle it, matching what Step 3 already did for the `Adapter 2` line); and the new `ci.yml:334` step, which establishes a **new** Perl 0.6.11 baseline and so should be watched on its first CI run — my local differential proves the fix moved no existing baseline, but only the job itself can confirm the new one.

**Caveat on freshness.** This audit describes the working tree as of 16:04, and the tree was being edited concurrently. If further edits have landed since, re-run `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings` and `cargo test --release` before trusting the V14 row.
