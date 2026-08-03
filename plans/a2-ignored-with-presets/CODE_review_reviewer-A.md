# CODE review — Reviewer A — `-a2` must always win over preset/auto-detected Read 2 adapters (#369)

**Reviewed:** working tree on `fix/a2-honoured-with-presets`, based on `dev` @ `4276810`. Nothing committed.
**Files:** `src/main.rs`, `src/cli.rs`, `tests/integration_adapter2.rs` (new), `.github/workflows/ci.yml`, `CHANGELOG.md`, `docs/src/content/docs/guide/flags.md`.
**Read-only:** no repo file was modified except this report. Scratch output went to the session scratchpad.

**Verification method key** — every finding below is labelled:
- **[ran]** — I executed the binary, the test suite, Perl 0.6.11, or a bash transcription of the CI step and observed the result.
- **[read]** — I traced the code and am confident from the source alone.
- **[suspect]** — flagged but not confirmed.

---

## Summary

The design is right and the implementation is small, correct on the central contract, and well tested. The caller-side override is the correct shape: Read 1 provably cannot move, and I confirmed that by reading the diff. The full suite (454 tests), `cargo fmt --all -- --check` and `cargo clippy --all-targets --release -- -D warnings` are all clean **[ran]**.

Two of the reviewer-brief's open questions came back **stronger than the plan claims**:

- **§13 deviation 1's suppression-detection claim is confirmed**, and the coupling is tighter than §13 states — see "Verified clean" §A.
- **The new CI step is not just plausible, it passes for real.** `cutadapt 5.2` (the exact CI pin) and `perl` are installed on this machine, so I ran both legs of the new validation step against Perl Trim Galore 0.6.11 extracted from the `0.6.11` tag. Both reads are **byte-identical** in both modes **[ran]**. My prior suspicion — that a homopolymer-rich adapter (`AAATCAAAAAAAC`) would expose alignment-DP tie-breaking differences against real cutadapt — is **refuted**.

Against that, one real user-visible defect ships in this state:

- **High:** on single-end input with `--small_rna` or `--bgiseq`, stderr now contradicts itself — a WARNING saying `-a2` is ignored, followed by a NOTE saying `-a2` displaced the preset's Read 2 default. This is exactly the V20 defect class the `cli.paired` display gate was added to close, re-introduced through the NOTE, which was not gated **[ran]**.
- **Medium:** the `A{N}` expansion notice is now printed twice (and N+1 times for N pairs), including on the previously-correct `-a` + `-a2` path **[ran]**.
- **Medium:** §13's headline claim "Zero branches of `resolve_adapter` were touched" is false as written — the `-a` branch *was* modified. The Read 1 guarantee survives, but it needs restating **[read]**.

---

## Issues by area

### 1. Logic

#### 🔴 High — single-end + a preset with an R2 default prints a WARNING and a contradictory NOTE

`apply_adapter2_override` (`src/main.rs:994`) has **no `cli.paired` guard**. It applies the override and returns the displaced default unconditionally. The display gate added at `src/main.rs:837` covers only the `Adapter 2 (Read 2):` line; the NOTE at `src/main.rs:847-851` is **not** gated.

Observed **[ran]** — `./target/release/trim_galore --small_rna -a2 AAATCAAAAAAAC test_files/BS-seq_10K_R1.fastq.gz -o …`:

```
WARNING: -a2/--adapter2 was given but is not used in this mode (it applies to Read 2 of a pair, and this is a single-end run). Ignoring.
Adapter: smallRNA (TGGAATTCTCGG)
NOTE: Read 2 adapter taken from -a2; the smallRNA default (GATCGTCGGACT) is not used.
```

The two messages assert opposite things about the same flag. Reachable via SE + `--small_rna`, SE + `--bgiseq`, and SE + auto-detected BGI.

This is not cosmetic in the plan's own terms. Plan §4 Step 3 and V20 exist *because* "SE output is not self-contradictory" was judged a real defect worth fixing; the fix gated the display line and left a message that says the same wrong thing more emphatically. `tests/integration_adapter2.rs:409` (`single_end_does_not_announce_a_read2_adapter`) asserts only the absence of `Adapter 2 (Read 2)` and uses `--illumina`, which has no R2 default — so no test covers this.

A secondary consequence **[read]**: in SE mode the override still writes the user's `-a2` into `TrimConfig::adapters_r2` (`src/main.rs:960`), so the single-end trimmer config carries a Read 2 adapter that nothing consumes. Harmless today — `trimmer.rs:112` keys on `is_r2`, `r2_adapter_count()` is only called from paired code paths (`trimmer.rs:400`/`711`, `parallel.rs:221`/`317`), and the SE report hardcodes `adapters_r2: Vec::new()` (`src/main.rs:1200`) — but it is dead state that a future reader could reasonably act on.

**Recommended fix** — one line, closes both, and keeps every other property:

```rust
fn apply_adapter2_override(…) -> Result<Option<String>> {
    if cli.adapter2.is_empty() || !cli.paired {
        return Ok(None);
    }
```

`Cli::validate` has already warned by this point, so nothing is silently lost. Malformed-`-a2` rejection is unaffected (that happens in `Cli::validate`, `src/cli.rs:941`). `--paired` is set on the interleaved-uBAM paired path too, so uBAM output is unaffected **[read]**. Add an SE + `--small_rna` case to `single_end_does_not_announce_a_read2_adapter` asserting the NOTE is absent.

#### 🟡 Medium — §13's "zero branches were touched" is false as stated

`resolve_adapter`'s `-a` branch **was** modified (`src/main.rs:1026-1030`): the `parse_adapter_specs(&cli.adapter2)` line was deleted and the returned R2 slot changed to `Vec::new()` **[read]**.

The Read 1 guarantee itself is intact, and I verified it by reading the whole diff: no branch's `adapters_r1` or `adapter_label` computation changed, `autodetect_poly_g` is untouched, and `apply_adapter2_override` takes `adapters_r1` by shared reference (`src/main.rs:996`) so it cannot mutate it. But the claim the PR narrative rests on should be the accurate one:

> No branch's Read 1 computation changed; the only edit inside `resolve_adapter` is the Read 2 slot of the `-a` branch.

Worth fixing in §13 before this becomes a PR description, since "structural rather than tested" is doing load-bearing work in the argument that R1 is safe.

#### 🟢 Low — the suppression check's coupling is invisible at the other end

`src/main.rs:1004` — `adapters_r1.len() == 1 && adapters_r1[0].1.is_empty()` — is correct (see §A), but nothing at `src/adapter.rs:219-227`, where the synthetic empty preset is constructed, says that `main.rs` reads emptiness as a suppression signal. A one-line note next to `seq: ""` would keep the two in sync. `adapters_r1.iter().all(|(_, s)| s.is_empty())` would also be marginally more robust than the `len() == 1` form, though I could not construct an input that distinguishes them.

#### 🟢 Low — `--consider_already_trimmed` cannot combine with any preset

`src/cli.rs:293` gives `consider_already_trimmed` `conflicts_with_all = ["illumina", "nextera", "small_rna", "stranded_illumina", "bgiseq"]`, so the suppression carve-out is reachable **only** on the auto-detect path **[ran — `--illumina --consider_already_trimmed 10000` exits with a clap conflict error]**. `-a` + `--consider_already_trimmed` *is* accepted and does bypass suppression **[ran]**, which is the pre-existing asymmetry §3.6 records as out of scope — correctly. No change needed; noting it because it narrows the carve-out's blast radius and is worth stating in §13.

### 2. Efficiency

#### 🟡 Medium — `A{N}` expansion notice printed twice, and once per pair thereafter

`-a2` is parsed in `Cli::validate` (`src/cli.rs:941`) and again in `apply_adapter2_override` (`src/main.rs:1014`). `parse_adapter_spec` emits `eprintln!("Adapter sequence {} expanded to {}", …)` (`src/adapter.rs:341`) as a side effect, so the notice repeats. Measured **[ran]**:

| Invocation | expansion notices | pre-fix |
|---|---|---|
| `--paired --illumina -a2 'A{10}' R1 R2` | 2 | 0 (flag discarded) |
| `--paired -a AGATCGGAAGAGC -a2 'A{10}' R1 R2` | **2** | **1** |
| `--paired --illumina -a2 'A{10}' R1 R2 R1' R2'` (2 pairs) | 3 | 0 |

The middle row is a regression on a path the plan asserts is unchanged (§3.5, V4). It is cosmetic — output bytes are unaffected — but it lands on exactly the feature the CHANGELOG advertises as newly working, and it scales with pair count.

**Recommended fix:** give `Cli::validate` a non-printing validation entry point. Splitting `parse_adapter_spec`'s body from its `eprintln!` (a `parse_adapter_specs_quiet`, or a `bool` parameter threaded through the three call sites) keeps the early-failure property §3.1 wanted and prints once. Do **not** drop the eager validation — that would re-introduce the 1 M-read scan before the failure.

The other cost of the double parse — a `file:` FASTA re-read once in `validate` plus once per input pair — is negligible and already true for `-a`.

### 3. Errors

Nothing found. Specifically checked and clean:

- **Double validation does not double-report.** `--paired --illumina -a2 ZZZQQQ` prints `Error: Adapter sequence must contain only DNA characters (A, C, G, T, N, X), got: 'ZZZQQQ'` exactly once, exit 1 **[ran]**. `Cli::validate` is `&self` and is called once from `src/main.rs:166` **[read]**.
- **`--clump_only` ordering is correct.** `--clump_only -a2 SEQ` prints only `Error: --clump_only does not trim; -a2/--adapter2 is not compatible`, exit 1, with **no** preceding warning **[ran]**. The hard error at `src/cli.rs:716` sits inside the `clump_only` block at ~line 700, well before the new block at `src/cli.rs:938`.
- **All four unusable-mode reasons fire and are true for their mode** **[ran]**: `--clock`, `--hardtrim5`, `--hardtrim3`, `--implicon` (note `--implicon` uses `require_equals`, so it is `--implicon=4`, not `--implicon 4` — my first attempt was a bad invocation, not a missing warning), and single-end.
- **`-a` + `-a2` is unchanged for every form.** Single sequence (via `oracle_equivalence`, passing), repeated `-a2`, `A{N}`, and `file:` all behave identically on the `-a` and `--illumina` paths — `Adapters R2 (2 sequences):` and a two-entry per-adapter breakdown in the R2 report on both **[ran]**. Structurally the pre- and post-fix expressions are the same `parse_adapter_specs(&cli.adapter2)` applied to the same input, and the new suppression check cannot fire on the `-a` path because that branch returns before auto-detection **[read]**.

### 4. Structure and style

#### 🟢 Low — comment style

`CLAUDE.md`'s `Code comments: one line, state the fact` rule:

| Site | Verdict |
|---|---|
| `src/main.rs:820` — `// #369 — a user -a2 wins over the preset/auto-detected Read 2 candidate.` | Compliant — one line, states the fact. |
| `src/main.rs:1003` — `// Suppression is the only way an R1 adapter reaches here with an empty sequence.` | Compliant, and the single most valuable comment in the diff. |
| `src/main.rs:836` — `// Read 2 adapters are only used in paired mode.` | Compliant. Becomes fully true once the High finding is fixed; today the *value* is still populated in SE. |
| `src/cli.rs:938-939` — two lines | At the stated maximum. Acceptable. |
| `src/main.rs:988-993` — the `apply_adapter2_override` doc comment | The third sentence ("trimming R2 while R1 is left alone would be asymmetric, and the mode announces that only quality trimming will happen") is a reasoning chain the rule pushes to the commit message. Existing doc comments in this file are similarly discursive and predate the rule, so this is a judgement call, not a defect. |
| `.github/workflows/ci.yml:335-338` — four comment lines | Over the limit for committed config, though `src/cli.rs:700-711` shows a twelve-line precedent nearby. |

#### 🟢 Low — no unit test of the helper, contrary to plan Step 5

Plan §4 Step 5: *"Keep a unit test of the helper as well, but not as a substitute."* There is none, and this is not listed among §13's three deviations. `apply_adapter2_override` is a private `fn` in `main.rs`, which has no `#[cfg(test)]` module — so the plan's own escape hatch ("or moved to the library if a unit test is wanted") applies. I agree with the priority (the integration matrix is what catches this bug class), but the deviation should be recorded rather than silent.

### 5. Test quality

The 13 tests are genuinely good — the oracle, the per-branch loop, and the poly-G scan-line count are all non-vacuous, and the `r2_report_line` rewrite of `suppression_wins_over_a2` correctly sidesteps the `split_whitespace` trap (`Reads with adapters` containing ` 0 (0.0%)` is the real behaviour, and I confirmed the report formats it that way) **[ran]**. Remaining gaps, all Low:

1. **`malformed_a2_is_rejected_on_preset_paths` (`:424`) asserts only `!ok`.** It would pass if the binary failed for an unrelated reason (a missing fixture, an output-collision preflight hit). Add a stderr assertion — `stderr.contains("Adapter sequence must contain only DNA")` for `ZZZQQQ`, `contains("Empty adapter sequence")` for `""`, `contains("Cannot open adapter FASTA")` for the missing file. I confirmed all three messages are the real ones **[ran]**.
2. **`unusable_a2_warns_with_a_mode_accurate_reason` (`:369`) covers 3 of 5 modes** — `--implicon` and `--hardtrim3` are absent. I verified both manually **[ran]**; adding them is two array rows (with `--implicon=8`, per `require_equals`).
3. **`repeated_a2_applies_on_a_preset_path` (`:448`) asserts on stderr only.** V22 promised the R2 report's per-adapter breakdown. I confirmed it is two entries **[ran]** — worth asserting, since the stderr line only proves the display saw two sequences, not that the stats were sized from `adapters_r2`.
4. **`read1_trimming_is_unaffected_by_the_r2_adapter` (`:213`) discards both `run()` results.** It fails loudly via the `expect("no R1 report")`, so this is a readability point, not a soundness one.
5. **No test cleans up its temp dir.** `tempdir()` (`:32`) removes a stale dir on entry but leaves the final one behind. Matches `integration_clump_only.rs`, so it is house style.

### 6. CI step

**Bash semantics: correct.** I transcribed the step verbatim into a `bash -e` script (GitHub Actions' default shell for `run:` on Linux) and instrumented the argv **[ran]**:

```
[mode='--illumina']  argc=4  argv=(--paired --illumina -a2 AAATCAAAAAAAC)
[mode='autodetect']  argc=3  argv=(--paired -a2 AAATCAAAAAAAC)
```

The empty iteration collapses to nothing under bash word-splitting, as intended. The zsh trap that produced the false result during implementation does not apply here.

**The grep cannot pass vacuously.** `grep -q 'a AAATCAAAAAAAC' …/BS-seq_10K_R2.fastq.gz_trimming_report.txt` targets the correct file (`naming::report_name` → `<input filename>_trimming_report.txt`) and matches exactly one line **[ran]**:

```
 9: Optional adapter 2 sequence (Read 2): 'AAATCAAAAAAAC'   ← no match (quote before A)
18: Command line parameters: … -O 1 -a AAATCAAAAAAAC …      ← the match
33: Sequence: AAATCAAAAAAAC; Type: regular 3'; …             ← no match ("e: " not "a ")
```

Neither `-a2 AAATCAAAAAAAC` nor `--adapter2 AAATCAAAAAAAC` matches the pattern `a AAATCAAAAAAAC` (a `2` intervenes), so an echoed command line cannot satisfy it. The R1 report does not match at all **[ran]**.

**The Perl comparison passes — verified for real, not inferred.** `cutadapt 5.2` (the exact CI pin at `ci.yml:278`) and `perl` are on this machine, so I extracted `git show 0.6.11:trim_galore` and ran both legs:

| mode | file | Perl 0.6.11 md5 | Rust md5 | |
|---|---|---|---|---|
| `--illumina` | `_val_1.fq.gz` | `d9246013dd91b9628b02e1d782ba2218` | same | MATCH |
| `--illumina` | `_val_2.fq.gz` | `fd2935a73a9673102ff51db1eb3df0a3` | same | MATCH |
| autodetect | `_val_1.fq.gz` | `d9246013dd91b9628b02e1d782ba2218` | same | MATCH |
| autodetect | `_val_2.fq.gz` | `fd2935a73a9673102ff51db1eb3df0a3` | same | MATCH |

The step's premise also checks out against the Perl source **[read]**: 0.6.11 dies on `-a` + preset (lines 3151-3164) but has no such guard for `-a2` + preset, and `$a2` is assigned only at lines 518/522 (`--polyA`), 563 (smallRNA), 573 (BGISEQ) and 585 (default `''`) — never by auto-detection. So both legs are legitimately comparable, and `if ($validate and $a2)` at line 1347 puts `-a $a2` on R2's cutadapt call exactly as we do.

#### 🟢 Low — the two legs are md5-redundant on this fixture

Auto-detection is *inconclusive* on `BS-seq_10K_R1` — Perl reports `count Nextera: 0, count Illumina: 0, count smallRNA: 0` and defaults to Illumina **[ran]** — so the autodetect leg produces byte-identical output to the `--illumina` leg. It still earns its place (it exercises a different resolution branch and its own grep), but a reader could mistake it for coverage of *detection*. One comment line would prevent that.

### 7. Invariants

- **`validation` job byte-identity: safe.** No accepted path that omits `-a2` can change. `Cli::validate`'s new block is guarded by `!self.adapter2.is_empty()` (`src/cli.rs:940`) and `apply_adapter2_override` returns `Ok(None)` immediately on the same condition (`src/main.rs:999`), so with no `-a2` the only reachable delta is the `-a` branch's R2 slot changing from `parse_adapter_specs(&[])` — which returns `Vec::new()` — to a literal `Vec::new()`. Identical **[read]**. The display gate at `src/main.rs:837` only affects stderr, which the validation job does not hash **[read]**. §13's V13 result (24 baseline files identical) is consistent with this.
- **No external runtime dependency added.** Byte comparison via `flate2::read::MultiGzDecoder` in the test's `gunzip` (`tests/integration_adapter2.rs:97`) is the right call — `flate2` is already a dependency, and full-content comparison is strictly stronger than an md5 for a test that can print the mismatch. Deviation 2 is sound.
- **`-D warnings` clean; `fmt --check` clean; 454 tests pass** (369 lib + 13 new + 72 other integration) **[ran]**.

### 8. User-facing text

Grammatical and accurate in every case I could fire, with the one exception in the High finding (the NOTE is accurate about the override but is emitted in a context where the override does not matter).

- `WARNING: -a2/--adapter2 was given but is not used in this mode (<reason>). Ignoring.` — plain register, matches the deprecation-warning house pattern at `src/cli.rs:959+`, no second person. The reason-neutral framing plus a per-mode cause does what §3.4 asked; I confirmed all five reasons are true for their mode **[ran]**.
- `WARNING: -a2/--adapter2 not applied — adapter trimming is suppressed for this library (--consider_already_trimmed). Ignoring.` — accurate. Fires once per input pair in a multi-pair run, consistent with other per-pair messages.
- `NOTE: Read 2 adapter taken from -a2; the <label> default (<seq>) is not used.` — reads correctly for every label that can reach it: `smallRNA`, `BGI/DNBSEQ` **[ran]**.
- `CHANGELOG.md` — factually correct on every claim I checked, including "affects every v2 release (v2.1.0 through v2.3.0)", the pair-retention explanation, and "Two documented examples in the adapter guide were consequently no-ops" (`docs/…/guide/adapters.md:31` and `:38` both pass `-a2` with no `-a`; both work now **[ran]**). It is ~50 lines against three neighbouring fix entries of 5-15 — deliberate per §7, and I would leave it, but "are correct as written from this release" is an awkward closing clause worth a re-word.
- `docs/…/guide/flags.md:12` — the hedge copies the `--clip_R2` phrasing as §4 Step 8 asked. It calls the flag's adapter "BGISEQ-500" where `src/adapter.rs:61` names the preset "BGI/DNBSEQ"; both appear in the Perl lineage, so this is only worth aligning if the maintainer cares.

---

## Verified clean (things the brief asked me to try to break, and could not)

### A. §13 deviation 1 — the suppression-detection claim is **confirmed**

The claim is that suppression is the only way an R1 adapter reaches `apply_adapter2_override` with an empty sequence. I attacked it from four directions **[read, exhaustive over `src/adapter.rs`]**:

| Candidate producer | Result |
|---|---|
| `parse_adapter_spec` case 3 (single) | `src/adapter.rs:328` bails on empty; `:335` bails on `A{0}`; `:344` validates characters. Cannot yield empty. |
| `parse_adapter_spec` case 2 (embedded ` -a `) | `:314` skips empty parts, `:320` bails if all were empty. Cannot yield empty. |
| `parse_adapter_spec` case 1 (`file:`) | `read_fasta_adapters` bails on an empty record (`:416`, `:434`) **and** bails if the file yields none (`:449`). So it cannot return an empty vec either. |
| The five presets | `src/adapter.rs:36-64` — every `seq` is non-empty. |
| Auto-detect, inconclusive | Falls back to `ILLUMINA`, non-empty (`:206-211`). |
| Auto-detect, suppressed | `src/adapter.rs:219-227` — `seq: ""`, and `to_adapter_vec()` (`:23-25`) always produces exactly one element. **The only producer.** |

So the `len() == 1 && [0].1.is_empty()` predicate is exactly equivalent to `DetectionResult.suppressed`, with no false positives *and* no false negatives. It is in fact tighter than §13 claims, because clap forbids `--consider_already_trimmed` with all five presets (`src/cli.rs:293`) — auto-detect is the only route in.

### B. Report-level SE contradiction — refuted

I hypothesised that gating the *display* but not the *value* would leak a Read 2 adapter into the single-end trimming report and JSON. It does not: `src/main.rs:1200` hardcodes `adapters_r2: Vec::new()` for the SE report config, so the SE text report is byte-identical with and without `-a2`, and `"adapters_r2": []` in the JSON **[ran, diffed]**. The paired JSON does carry the override as §7 item 6 predicted: `"adapters_r2": [{"name": "adapter_1", "sequence": "AAATCAAAAAAAC"}]` **[ran]**.

### C. `--poly_a` interaction — none

Unlike Perl, `--poly_a` in this port is a `TrimConfig` boolean handled in `src/trimmer.rs:197-205` and never touches the adapter lists, so the override cannot clobber a poly-A R2 adapter **[read]**. §12's note on this is correct.

---

## Recommendations, prioritized

### Critical
None.

### High
1. **Guard `apply_adapter2_override` on `cli.paired`** (`src/main.rs:999`) so single-end runs stop emitting the contradictory NOTE and stop carrying a dead R2 adapter in `TrimConfig`. Extend `single_end_does_not_announce_a_read2_adapter` (`tests/integration_adapter2.rs:409`) with an SE + `--small_rna` case asserting no NOTE — the existing `--illumina` case cannot catch it.

### Medium
2. **Stop double-printing the `A{N}` expansion notice.** Add a non-printing validation entry point for `Cli::validate` (`src/cli.rs:941`) rather than dropping the eager validation, which is load-bearing for the pre-scan failure. Note this also restores the pre-fix behaviour of `-a` + `-a2 'A{N}'`, which currently prints twice.
3. **Restate §13's Read 1 guarantee accurately** — the `-a` branch was edited; the guarantee is that no branch's Read 1 computation changed. This wording is headed for a PR description.

### Low
4. Assert on stderr content in `malformed_a2_is_rejected_on_preset_paths` (`tests/integration_adapter2.rs:424`) so it cannot pass on an unrelated failure.
5. Add `--implicon=8` and `--hardtrim3` rows to `unusable_a2_warns_with_a_mode_accurate_reason` (`:369`); both work, neither is covered.
6. Assert the two-entry R2 per-adapter breakdown in `repeated_a2_applies_on_a_preset_path` (`:448`) — that is what V22 actually promised.
7. Add a one-line note at `src/adapter.rs:219-227` recording that `main.rs` treats an empty R1 sequence as the suppression signal, so the coupling verified in §A cannot be broken silently.
8. Note in `.github/workflows/ci.yml:334` that auto-detection is inconclusive on this fixture and therefore defaults to Illumina, so the second leg is branch coverage rather than detection coverage.
9. Record in §13 that plan Step 5's unit test of the helper was deliberately not written, and why (`main.rs` has no test module).
10. Consider trimming the `apply_adapter2_override` doc comment's rationale sentence (`src/main.rs:991-993`) and the four-line `ci.yml` comment per `CLAUDE.md`'s comment rule. Judgement call — nearby code predates the rule.
11. Re-word the CHANGELOG's "are correct as written from this release".

---

## Commands run

```
cargo build --release
cargo test --release                                    # 454 pass, 0 fail
cargo fmt --all -- --check                              # clean
cargo clippy --all-targets --release -- -D warnings     # clean
bash <verbatim transcription of ci.yml:334-356>         # argv instrumented, greps pass
git show 0.6.11:trim_galore                             # Perl source, $a2 + die audit
perl trim_galore_perl.pl --paired [--illumina] -a2 AAATCAAAAAAAC …   # cutadapt 5.2, both legs
./target/release/trim_galore …                          # ~15 manual invocations (SE/PE, all 5 presets,
                                                        #   autodetect, A{N}, file:, repeated -a2,
                                                        #   malformed, clump_only, clock, hardtrim3,
                                                        #   implicon, consider_already_trimmed)
```
