# PLAN v2 — `-a2` must always win over preset and auto-detected Read 2 adapters (#369)

**Issue:** [#369](https://github.com/FelixKrueger/TrimGalore/issues/369) — *`-a2` adapter ignored when `--illumina` is used*, reported by @MathieuUm
**Branch:** `fix/a2-honoured-with-presets` off `dev`
**Reported impact:** ~20% drop in mapping rate and properly-paired alignments after upgrading 0.6.11 → 2.3.0, traced by the reporter to Read 2 being trimmed with the Illumina adapter instead of their `-a2` sequence.

**Scope decisions taken with the user before writing this plan:**

- `-a` + preset stays **permissive** (`-a` wins, no error). Perl 0.6.11 rejects that combination, but it is the workaround now recommended on #369 and possibly in use elsewhere; turning a working command into a hard error is not acceptable here. The divergence is documented, not fixed.
- A user `-a2` that displaces a preset's own R2 default **announces itself**.

---

## 1. Goal

Restore the 0.6.11 contract: **a user-supplied `-a2` always determines the Read 2 adapter.** Presets and auto-detection supply a Read 2 default only when `-a2` is absent.

Secondary, same defect class: stop silently discarding `-a2` when it cannot apply (single-end input).

Non-goal: changing which adapter Read 1 gets, in any mode. Read 1 *per-read trimming* is correct today and must not move.

**Read 1's output FILE will change**, however, for every invocation where `-a2` newly applies. Measured on `test_files/BS-seq_10K_R{1,2}.fastq.gz`: trimming R2 with the user's adapter instead of `AGATCGGAAGAGC` takes the joint pair-length filter from 4 dropped pairs (0.0%) to 139 (1.4%), so `_val_1.fq.gz` differs even though R1's own cutadapt-stage figures are identical (`Reads written (passing filters): 10,000` in both). v1 of this plan asserted R1 byte-identity as its central safety property; that was wrong, and §9 V7 is re-specified accordingly. This also makes the reporter's ~20% mapping-rate drop more explicable — it is a pair-retention effect, not only R2 mis-trimming.

---

## 2. Context

### 2.1 Root cause

`src/main.rs:985` `resolve_adapter` is an early-return chain, and **only the first branch reads `cli.adapter2`**:

```rust
if !cli.adapter.is_empty() {                                     // reads -a2 ✓
    let adapters_r2 = adapter::parse_adapter_specs(&cli.adapter2)?;
    …
}
if cli.nextera           { return (…, Vec::new(), …) }            // -a2 dropped
if cli.small_rna         { return (…, SMALL_RNA.to_r2_vec(), …) } // -a2 overridden
if cli.stranded_illumina { return (…, Vec::new(), …) }            // -a2 dropped
if cli.bgiseq            { return (…, BGISEQ.to_r2_vec(), …) }    // -a2 overridden
if cli.illumina          { return (…, Vec::new(), …) }            // -a2 dropped
// auto-detect           → detection.adapter.to_r2_vec()          // -a2 dropped
```

An empty `adapters_r2` then falls back to Read 1's adapter at `src/trimmer.rs:112`:

```rust
let adapters = if is_r2 && !config.adapters_r2.is_empty() { &config.adapters_r2 } else { &config.adapters };
```

which is why the reporter saw `AGATCGGAAGAGC` on Read 2 rather than untrimmed reads. That fallback is correct for an *unspecified* `-a2` and wrong for a *discarded* one, and the current code cannot distinguish the two cases.

### 2.2 Reproduced scope — wider than reported

`-a2` is honoured **only** when `-a` is also given. Verified on `dev` @ `4276810` by reading the `Command line parameters` line of each R2 trimming report, using the reporter's `-a2 AAATCAAAAAAAC`:

| Invocation | R2 adapter actually used | |
|---|---|---|
| `-a AGATCGGAAGAGC -a2 SEQ` | `AAATCAAAAAAAC` | ✅ |
| `--illumina -a2 SEQ` | `AGATCGGAAGAGC` | ❌ reported |
| `--nextera -a2 SEQ` | `CTGTCTCTTATA` | ❌ |
| `--stranded_illumina -a2 SEQ` | `ACTGTCTCTTATA` | ❌ |
| `--small_rna -a2 SEQ` | `GATCGTCGGACT` | ❌ |
| `--bgiseq -a2 SEQ` | `AAGTCGGATCGTAGCC…` | ❌ |
| auto-detect (no adapter flag) `-a2 SEQ` | `AGATCGGAAGAGC` | ❌ |

The last row is the most likely to affect others: the **default** path discards `-a2` too.

### 2.3 Not a v2.3.0 regression

The `--illumina` branch is byte-identical at v2.1.0, v2.2.0 and v2.3.0 (`git show <tag>:src/main.rs`). Every v2.x release is affected, so downgrading within v2 is not a workaround. The reporter perceived it as a 2.3.0 change only because they upgraded from 0.6.11 directly.

### 2.4 The 0.6.11 contract, from the tagged source

`git show 0.6.11:trim_galore`, lines 560-578:

```perl
unless (defined $a2) { $a2 = 'GATCGTCGGACT' }   # small RNA, PE only
unless (defined $a2) { $a2 = 'AAGTCGGATCGTAGCC…' }   # BGISEQ, PE only
```

`unless (defined $a2)` is the whole contract: the preset fills a default, a user value wins. The same idiom appears for the small-RNA length cutoff (`# user defined length cutoff wins over auto-detection`), so this is the script's deliberate precedence convention rather than an accident.

Two further Perl behaviours, both relevant and both *divergent* from v2.x in the opposite direction:

- **`-a` + preset → `die`** (line ~3131: *"You can't supply an adapter sequence AND use the Illumina universal adapter sequence. Make your choice."*). v2.x accepts it and lets `-a` win. **Left as-is** per the scope decision — this is the #369 workaround.
- **`-a2` without `--paired` → `die`** (line ~3149: *"An optional adapter for read 2 of paired-end files requires '--paired' to be specified as well!"*). v2.x silently accepts and ignores it — verified: `trim_galore -a2 ACGTACGT single.fq.gz` runs to completion with no warning. Addressed here as a **warning**, not an error (§3.4).

### 2.5 Why the Perl-parity CI gate missed it

`grep -cE '(^|[[:space:]])(-a2|--adapter2)([[:space:]=]|$)' .github/workflows/ci.yml` returns **0**. (The unanchored `grep 'a2\|adapter2'` that v1 of this plan quoted returns 8 spurious hits inside `sha256sum` and the pinned `dtolnay/rust-toolchain@29eef336…` SHA — right conclusion, misleading command.) The `validation` job never passes `-a2`. (Its scope is also larger than `CLAUDE.md`'s "SE, PE, hardtrim5, clock, demux" summary suggests — roughly 30 `Validate …` steps and 24 byte-identity assertions, including `--clump_only` byte-identity and uBAM `@PG` parity. That summary is out of date; see §9 V13.) A parity gate that never exercises a flag cannot detect that the flag is ignored. §4 Step 6 closes this.

### 2.6 Existing display code already does most of the announcing

`src/main.rs:833-842` already prints, whenever `adapters_r2` is non-empty:

```
Adapter 2 (Read 2): <seq>
```

Once the fix populates `adapters_r2` from `-a2`, that line appears with no new code. Two consequences:

- The §3.3 NOTE only needs to state the *override* — that a preset default was displaced — not repeat the sequence.
- The **absence** of that line was itself available evidence that `-a2` had been dropped. Worth remembering as a diagnostic: if `-a2` was given and no `Adapter 2 (Read 2):` line appears, it was not applied.

---

## 3. Behavior

### 3.1 Resolution order (the fix)

There is exactly **one** place where the Read 2 adapter is decided, and it is in the **caller**, not inside `resolve_adapter`.

`resolve_adapter` keeps its seven branches exactly as they are and its return type unchanged; its `adapters_r2` becomes a documented *candidate*. `setup_trimming` — its only caller (`main.rs:817`) — then applies the override immediately before the display block that has to announce it:

1. `resolve_adapter` yields Read 1 plus a *candidate* Read 2 default, by the existing precedence: explicit `-a` → `--nextera` → `--small_rna` → `--stranded_illumina` → `--bgiseq` → `--illumina` → auto-detect.
2. In the caller: **if `cli.adapter2` is non-empty it replaces the candidate**, except under adapter-trimming suppression (§3.6).
3. If it is empty, the candidate stands (which may be empty → Read 2 falls back to Read 1's adapter at `trimmer.rs:112`, unchanged).

**Why the caller and not a restructure.** v1 proposed reshaping all seven branches from early returns into a single value, and booked the resulting "restructure touches Read 1 incidentally" as this plan's only Medium risk. That risk was self-inflicted. Applying the override one frame up is ~10 lines, touches **zero** existing branches — so Read 1 provably cannot move — and removes the §5 signature change entirely, because the caller computes the displacement itself. §3.1's guarantee is unchanged: there is still exactly one decision point, and no branch can forget, because no branch is involved.

Sketch:

```rust
let (adapter_label, adapters_r1, mut adapters_r2, autodetect_poly_g) =
    resolve_adapter(cli, input_file)?;

// user -a2 wins over preset/auto-detected R2 defaults
let displaced = apply_adapter2_override(cli, &mut adapters_r2)?;
```

The only tidy-up is deleting the now-redundant `parse_adapter_specs(&cli.adapter2)` from the `-a` branch (`main.rs:988`), a one-line deletion.

**Parse `-a2` early.** `cli.adapter2` must be validated in `Cli::validate`, not at the override site, so a malformed `-a2` fails before the auto-detect branch spends a 1 M-read scan on an input that is going to be rejected anyway.

`-a2` must go through `adapter::parse_adapter_specs` on every path, so `A{N}` shorthand, repeated `-a2`, and `file:adapters.fa` behave identically to the `-a` branch. Today only that branch parses it.

### 3.2 Truth table

| `-a2` given | Preset / detection R2 default | Read 2 adapter |
|---|---|---|
| yes | none (`--illumina`, `--nextera`, `--stranded_illumina`, most auto-detect) | **`-a2`** |
| yes | present (`--small_rna`, `--bgiseq`, auto-detected BGI) | **`-a2`**, default displaced → NOTE (§3.3) |
| no | none | empty → falls back to R1 adapter (unchanged) |
| no | present | the default (unchanged) |

Read 1 is unaffected in all four rows.

### 3.3 Override NOTE

When `-a2` displaces a non-empty preset default, emit one line after the existing adapter display:

```
NOTE: Read 2 adapter taken from -a2; the <preset> default (<seq>) is not used.
```

Only when a default was actually displaced — `--illumina -a2 SEQ` displaces nothing, so the existing `Adapter 2 (Read 2):` line alone is correct there. Perl warns in the mirror-image case (when it *sets* the default: *"Setting the Illumina smallRNA 5' adapter as adapter 2"*), so announcing the displacement matches the house register.

### 3.4 `-a2` without `--paired`

Warn whenever `-a2` is given but cannot be used. That is broader than `!cli.paired`: the specialty modes (`--hardtrim5/3`, `--clock`, `--implicon`) bypass the trimming pipeline entirely and never call `resolve_adapter`, so `--paired --clock -a2 SEQ` discards the flag too — verified, exit 0 with no mention of an adapter. A `!cli.paired` condition would stay silent there while claiming to close the silent-discard class.

The message must state a reason that is true in every case it fires. `"...is ignored without --paired"` is false for `--hardtrim5 -a2 SEQ R1 R2`, where the flag is ignored because hardtrim does no adapter trimming at all. Prefer reason-neutral phrasing plus the specific cause:

```
WARNING: -a2/--adapter2 was given but is not used in this mode (<reason>). Ignoring.
```

matching the house pattern at `cli.rs:939-963`. Sited in `Cli::validate`, which already ends with a deprecation-warning block, so this is not a new kind of side effect there.

Perl rejects instead. A warning is chosen because these invocations are otherwise valid and complete successfully — rejecting would fail runs that currently produce correct output for every flag except the ignored one. (v1 justified this as "should not convert working invocations into failures"; §3.6 shows that principle does not hold universally, so the argument stands on the narrower ground above.) `--clump_only` already hard-errors on `-a2` (`cli.rs:716`) and is left alone.

Recorded as an open item (§10) in case the maintainer prefers Perl's hard error.

### 3.6 Newly-fatal input, and the suppression carve-out

Two consequences of §3.1 that v1 did not own.

**`-a2` is now validated on every path, so malformed values become fatal.** Today `--illumina -a2 ZZZQQQ` exits 0 with the flag silently discarded, while `-a … -a2 ZZZQQQ` already errors with *"Adapter sequence must contain only DNA characters"*. After the fix both error. Affected forms: invalid characters, `-a2 ""` (empty), and `-a2 file:missing.fa`.

This is a deliberate change and is **the right one** — silently ignoring a malformed adapter is worse than rejecting it, and it makes `-a2` consistent with `-a`. But it contradicts the "do not convert working invocations into failures" argument v1 used in §3.4, so that argument is restated on its own merits there, and this becomes a listed user-visible change (§7) with its own CHANGELOG sentence and V-row.

**`--consider_already_trimmed` suppression wins over `-a2`.** When that flag decides a library is already trimmed, `adapter.rs:219-227` returns a synthetic preset with `seq: ""` and `seq_r2: None`, and both reads go untrimmed. Applying the override there would adapter-trim R2 while leaving R1 untrimmed — an asymmetry no user would predict, under a banner that says *"Only quality trimming will be carried out"*. So the override is **skipped** when `DetectionResult.suppressed` is set (`adapter.rs:78`, already in scope), and a warning says `-a2` was not applied.

Accepted cost: one documented exception to "`-a2` always wins". Noted asymmetry that remains: `-a SEQ --consider_already_trimmed N` already bypasses suppression for Read 1 today, because the `-a` branch returns before auto-detection runs. That is pre-existing and out of scope.

### 3.5 Unchanged behaviour, stated so it is testable

- Read 1 adapter selection, precedence, and per-read trimming. (Not the R1 output *file* — see §1.)
- `-a` + preset: accepted, `-a` wins, no error and no new warning.
- `--small_rna`'s length-cutoff side effect (default 18 rather than 20) — driven by the Read 1 adapter, independent of `-a2`.
- `-a2` **absent** (empty `Vec`) on a preset with an R2 default: still gets the default. Note `-a2 ""` is a *different* case and now errors — see §3.6.
- Byte-identity on every CI validation-matrix path, none of which passes `-a2`.

---

## 4. Implementation outline

1. **Add `apply_adapter2_override` and call it from `setup_trimming`** (`main.rs:817-819`), per §3.1. `resolve_adapter` is not restructured: its seven branches and its 4-tuple stay exactly as they are, and only its doc comment changes. Delete the now-redundant `parse_adapter_specs(&cli.adapter2)` from the `-a` branch (`main.rs:988`).

2. **Validate `-a2` in `Cli::validate`** so malformed input fails before the auto-detect scan (§3.1), then parse once in the override helper for the override decision and the NOTE.

3. **Emit the override NOTE** (§3.3) after the existing adapter display at `main.rs:833-842`. The displacement comes straight from the helper's return value. Also **gate the `Adapter 2 (Read 2):` display on `cli.paired`**, which fixes the pre-existing single-end case where that line prints for a flag that is never used (§3.4, V20).

4. **Add the unusable-`-a2` warning** (§3.4) in `Cli::validate`, covering single-end *and* the specialty modes, with a per-mode reason. Add the suppression carve-out and its warning (§3.6) at the override site, keyed on `DetectionResult.suppressed`.

5. **Tests** (§9). A unit test of the override helper **cannot catch this bug class** — the defect was never "the helper computes the wrong answer", it was "six branches never called it". So the primary coverage is per-branch end-to-end, in a new `tests/integration_adapter2.rs` following `tests/integration_clump_only.rs`: run the binary, assert the R2 report's `Command line parameters` line for all seven branches, plus V16-V23. Keep a unit test of the helper as well, but not as a substitute. (`main.rs` has no `#[cfg(test)]` module; the helper can be tested via the integration binary or moved to the library if a unit test is wanted.)

6. **CI matrix case** in `.github/workflows/ci.yml`'s `validation` job: a PE run with `--illumina -a2 <seq>` md5-compared against Perl 0.6.11. This is the gate that should have caught the defect, and its absence is why the defect shipped in three releases.

7. **CHANGELOG** — append to the existing `#### Fixes` block under `### Unreleased`. State that it affects all v2.x, not just 2.3.0; that `-a2` was ignored on every path except explicit `-a`; and that Read 1 output is unchanged.

8. **Docs** — the situation is worse than "check": two *recommended* examples have been no-ops for three releases. Verified edit list:
   - `guide/adapters.md:31` — presented as **"The cleanest (v2.x)"**: `--paired -a2 … -a2 … -a2 … -n 3` with no `-a`. Discards all three.
   - `guide/adapters.md:38` — same defect in the embedded-string form.
   - `guide/adapters.md:50`, `guide/flags.md:27`, `guide/flags.md:39` — claim `A{N}` expansion works for `-a2`; false today unless `-a` is also given.
   - `guide/flags.md:12` — `--small_rna` "auto-sets `--adapter2`"; needs the "unless you supply your own" hedge. `flags.md:13` already uses that phrasing for `--clip_R2` — copy it.
   - `reference/migration.md:26` — sells repeatable `-a`/`-a2` as a v2 improvement; true only with `-a`.

   All become correct *as written* once the fix lands, so most need no edit beyond `flags.md:12`'s hedge — but that they were wrong is worth a CHANGELOG sentence.

9. **Comment on #369** linking the PR, and close it on merge (manual — auto-close does not fire on `dev`).

**Comment style:** per the current `CLAUDE.md`, committed comments are **one line, stating the fact**; two is the maximum and needs a reason. Reasoning belongs in the commit message and this plan, not the source. The `-a2` override site needs at most `// user -a2 always wins over preset/auto-detected R2 defaults` — not the history.

---

## 5. Signature

**No signature change.** v1 proposed growing `ResolvedAdapter` to a 5-tuple so the caller could announce a displacement; the caller-side override (§3.1) computes it locally instead, so `resolve_adapter` keeps its existing 4-tuple and all seven branches stay untouched.

One small private helper in `main.rs`:

```rust
/// Apply a user `-a2` over the preset/auto-detected candidate. Returns the
/// displaced default, if any, for the caller's NOTE.
fn apply_adapter2_override(
    cli: &Cli,
    candidate: &mut AdapterList,
) -> Result<Option<String>>
```

Returns `Option<String>` rather than `(preset_name, seq)`: the preset name is already `adapter_label` at the call site, so a tuple would duplicate it.

`resolve_adapter`'s doc comment gains one line noting that `adapters_r2` is a *candidate*, subject to override by the caller.

## 6. Efficiency

No measurable change to the resolution code: one extra `parse_adapter_specs` over a vector that is almost always empty or length 1, once per input file, outside any per-read path. Nothing is restructured, so no I/O moves.

Two qualifications v1 stated too broadly:

- **R2 adapter length becomes user-controlled on the preset and auto-detect paths.** The adapter DP is O(read_len × adapter_len), and `alignment.rs:19`/`:272` cap the Myers bit-parallel prefilter at 64 bp — longer patterns skip it and run the full scalar DP on every read. A `-a2 T{150}` (the Perl poly-A idiom) therefore costs roughly an order of magnitude more per R2 read than the 13 bp adapter it replaces. Reachable today only via `-a`; the fix makes it reachable everywhere. Not a reason to change the design, but "no measurable change" is wrong unqualified.
- **`-a2 file:…` is re-read once per input file**, since resolution runs per input. Already true for `-a`, so no regression — but a 20-pair run opens the FASTA 20 times.

---

## 7. Integration

**Reads:** `cli.adapter`, `cli.adapter2`, the five preset flags, `cli.paired`, plus the auto-detection result. **Writes:** nothing beyond stderr.

**Order:** unchanged — `resolve_adapter` is called once per input from `setup_trimming` (`main.rs:818`).

**Behaviour changes visible to users:**

1. `-a2` is now applied on every path (the fix). **Trimmed output changes for anyone currently combining `-a2` with a preset or with auto-detection** — for *both* reads, since the joint pair-length filter retains a different set of pairs (§1). That is the point of the fix, but it means output differs from v2.1.0–v2.3.0 for those invocations. Prominent in the CHANGELOG.
2. A new `Adapter 2 (Read 2): <seq>` line appears for those invocations, from existing display code (§2.6).
3. A NOTE when a preset R2 default is displaced (§3.3).
4. A WARNING when `-a2` is given but unusable — single-end or a specialty mode (§3.4).
5. **A new non-zero exit** for malformed `-a2` on preset and auto-detect paths, which previously ran to completion with the flag discarded (§3.6).
6. `parameters.adapters_r2` in the JSON report (`src/report.rs:813`) changes from `[]` to the user's sequence on the `--illumina`/`--nextera`/`--stranded_illumina`/auto-detect paths.
7. The `--output-format ubam` paths inherit all of the above: `setup_trimming` is the shared entry point for both FASTQ and uBAM dispatch, so `--paired --illumina -a2 SEQ --output-format ubam` changes too.

No change for `-a … -a2 …`, for any invocation without `-a2`, or to Read 1 in any mode.

**Downstream:** `trimmer.rs:112`'s fallback and `TrimConfig::r2_adapter_count` (`trimmer.rs:49`) both key on `adapters_r2.is_empty()` and need no change — they simply now see the user's adapter where they previously saw the preset's or nothing.

---

## 8. Assumptions

- **A1.** Read 1 resolution and preset precedence are correct today and must not move. The fix touches only the R2 slot. Testable via V7.
- **A2.** `unless (defined $a2)` in 0.6.11 is the authoritative contract — verified in the tagged source, and corroborated by the reporter's own 0.6.11 report line showing `-a AAATCAAAAAAAC` for R2 under `--illumina`. Two independent sources.
- **A3.** `-a2`'s `A{N}` / `file:` / repeat handling should match `-a`'s exactly. Currently it only does so on the `-a` branch, because that is the only branch that parses it.
- **A4.** `--small_rna`'s length-cutoff default is a function of the Read 1 adapter, not `-a2`, so the fix cannot disturb it. Verify rather than assume (V8) — the two are adjacent in the Perl source and could plausibly have been coupled in the port.
- **A5.** No CI validation-matrix case passes `-a2` (verified by grep), so the fix cannot move any md5 in that job. The new case in Step 6 establishes a *new* baseline rather than changing an existing one.
- **A6.** `-a` + preset remains permissive. This is a deliberate, user-confirmed divergence from Perl, and the #369 reply already tells the reporter the workaround is sound — so it must keep working.

---

## 9. Validation

Run from the crate root; `cargo build --release` first (`cargo clippy` leaves the binary stale).

| # | Verify | How | Expected |
|---|---|---|---|
| V1 | The reported case is fixed | `--paired --illumina -a2 AAATCAAAAAAAC R1 R2`, read the R2 report's `Command line parameters` | `-a AAATCAAAAAAAC` |
| V2 | All five presets honour `-a2` | Repeat V1 for `--nextera`, `--small_rna`, `--stranded_illumina`, `--bgiseq` | `-a2` sequence in every R2 report |
| V3 | Auto-detect honours `-a2` | `--paired -a2 SEQ R1 R2`, no adapter flag | `-a2` sequence on R2; R1 unchanged from a no-`-a2` run |
| V4 | `-a` + `-a2` unchanged | `-a AGATCGGAAGAGC -a2 SEQ` | Identical output to `dev` before the fix — this path was already correct |
| V5 | `-a2` absent → preset default intact | `--small_rna` and `--bgiseq` with **no** `-a2` | `GATCGTCGGACT` / BGISEQ R2 respectively; byte-identical to pre-fix |
| V6 | `-a2` absent, no R2 default → R1 fallback | `--illumina`, no `-a2` | R2 uses `AGATCGGAAGAGC`; byte-identical to pre-fix |
| V7a | **Read 1 per-read trimming never moves** | For V1–V6, compare the R1 report's `Reads with adapters`, `Reads written (passing filters)` and `Total written (filtered)` bp, before and after | Unchanged in all cases. This is the R1 invariant that actually holds |
| V7b | R1 **file** byte-identity, where `-a2` does not newly apply | md5 R1 output for V4, V5, V6 and V13 only | Identical. Deliberately **not** asserted for V1–V3, where R1 changes by design (§1) |
| V7c | Reports, not just FASTQ | md5 `*_trimming_report.txt` alongside the `.fq.gz` in V7b | Identical. A reshape returning the right sequence with the wrong *label* is invisible to a FASTQ md5 |
| V8 | `--small_rna` length cutoff unaffected by `-a2` | `--small_rna` with and without `-a2`; check the reported length cutoff | 18 in both cases (A4) |
| V9 | NOTE fires only on real displacement | `--small_rna -a2 SEQ` vs `--illumina -a2 SEQ` | NOTE present for small_rna, absent for illumina |
| V10 | SE `-a2` warns | `-a2 SEQ single.fq.gz` | Warning on stderr, exit 0, output byte-identical to omitting `-a2` |
| V11 | `-a` + preset still permissive | `--illumina -a AAAACCCCGGGG` | Exit 0, `-a` wins, no error and no new warning (A6) |
| V12 | `A{N}` / `file:` / repeated `-a2` work on preset paths | `--illumina -a2 'A{10}'`; `--illumina -a2 S1 -a2 S2` | Expanded to `AAAAAAAAAA`; both sequences present for R2 |
| V13 | Validation matrix unmoved | Run the **whole** `validation` job (~30 steps, 24 byte-identity assertions) on a branch push — not "the five invocations", which undercounts it | All pass (A5) |
| V14 | Full suite + lint | `cargo test`; `cargo fmt --all -- --check`; `cargo clippy --all-targets --release -- -D warnings` | All pass, zero warnings |
| V15 | The new CI case actually fails pre-fix | Run Step 6's invocation against a pre-fix binary | Must **differ** from Perl. A gate that passes before the fix guards nothing. **Known satisfiable:** `-a2 AAATCAAAAAAAC` moves R2 `Reads with adapters` from 4,882 (48.8%) to 5,444 (54.4%) on the standard fixture, so the sequence discriminates |
| **V16** | **Differential oracle — the primary check.** Post-fix, `--illumina -a2 SEQ` and `-a AGATCGGAAGAGC -a2 SEQ` are the same invocation in everything reaching the FASTQ output: same R1 sequence, same R2 sequence, same length cutoff (keyed on the R1 *sequence* at `main.rs:848-856`), same poly-G decision | md5 both reads' output for both invocations | **Identical.** Needs no Perl, pins both reads at once, and is strictly stronger than V1–V3's report-line assertions. Same for the auto-detect path, which selects Illumina on this fixture |
| **V17** | Malformed `-a2` is fatal on preset paths (§3.6) | `--illumina -a2 ZZZQQQ`; `--illumina -a2 ""`; `--illumina -a2 file:missing.fa` | Non-zero exit with the same message the `-a` path already gives |
| **V18** | Suppression wins over `-a2` (§3.6) | `--paired --consider_already_trimmed 10000 -a2 SEQ` | Both reads untrimmed, symmetric; warning that `-a2` was not applied; **no** `Adapter 2 (Read 2):` line |
| **V19** | Specialty modes warn rather than silently discard (§3.4) | `--paired --clock -a2 SEQ`; `--hardtrim5 30 -a2 SEQ R1 R2`; `--paired --hardtrim3 5 -a2 SEQ` | Warning in each, with a reason that is true for that mode; output otherwise unchanged |
| **V20** | SE output is not self-contradictory | `--illumina -a2 SEQ single.fq.gz` | Warning present, **and** no `Adapter 2 (Read 2):` line. v1's V10 checked only the FASTQ, which is byte-identical either way, so it would have passed while stderr contradicted itself |
| **V21** | **The poly-G piggyback survives the change** | Count `Scanning for poly-G content` lines per branch, and compare the `Poly-G trimming:` verdict before/after | **0** for auto-detect (piggybacked on the adapter scan), **1** for each preset. Verified as the current behaviour. Invisible to every other check: dropping the piggyback costs a silent second 1 M-read pass with byte-identical output, and wrongly setting it disables poly-G silently — and the standard fixture is 2/10000, below threshold, so R1 md5 cannot see either |
| **V22** | Repeated `-a2` sizes the R2 stats correctly | `--illumina -a2 S1 -a2 S2` | R2 report shows a two-entry per-adapter breakdown. `TrimConfig::r2_adapter_count` (`trimmer.rs:49`) switches from `adapters.len()` to `adapters_r2.len()` once R2 is non-empty — a shape reachable today only via `-a` |
| **V23** | uBAM output inherits the fix | `--paired --illumina -a2 SEQ --output-format ubam` | R2 adapter is the `-a2` sequence; check the report, since paired uBAM output is one interleaved BAM |

V16, V21 and V15 are the ones that would quietly not hold. V21 in particular is the only check that can see the poly-G piggyback, and V15 encodes the lesson of the defect itself: a plausible-looking parity job that never exercises the flag guards nothing.

---

## 10. Questions or ambiguities

**Resolved with the user:**

- **R1.** `-a` + preset stays permissive (not Perl's hard error) — it is the #369 workaround.
- **R2.** A displaced preset R2 default is announced with a NOTE.

**Open (assumption taken, no blocker):**

1. **SE `-a2` warns rather than errors** (§3.4). Perl rejects. Warning chosen for consistency with R1's permissiveness. Escalate to a hard error if preferred.
2. **Whether to backport.** The defect affects every v2 release (v2.1.0 through v2.3.0; there is no v2.0.0 tag). Taken: fix on `dev` for the next release, no patch releases for older lines — consistent with how this project has shipped fixes so far.
3. **Whether the fix warrants a minor rather than patch bump.** It changes trimmed output for `-a2`-plus-preset users, which is a behaviour change even though it is a bug fix. Taken: note it prominently in the CHANGELOG and let the release decide, since v2.4.0 is already the next cut.

**Critical:** none outstanding.

---

## 11. Self-Review

**Efficiency.** Nothing in a per-read path; one extra parse of a near-always-empty vector per input file.

**Logic.** Traced all seven resolution paths (five presets, explicit `-a`, auto-detect) and confirmed exactly one reads `cli.adapter2` today. Checked the two downstream consumers of `adapters_r2` (`trimmer.rs:112` fallback, `TrimConfig::r2_adapter_count` at `:49`) — both key on `is_empty()` and need no change. Confirmed the display code at `main.rs:833` already announces a non-empty R2 adapter, so the fix gets that for free.

**Adjusted while writing.** Three changes: the NOTE shrank once I found `main.rs:833` already prints the sequence, so it now states only the displacement; the SE `-a2` gap was added after finding it silently accepted, since it is the same defect class; and V15 was added because a new CI case that passes before the fix would be worthless.

**Edge cases enumerated.** `-a2` with `A{N}` shorthand, `file:` spec, and repeats on preset paths (V12 — these work today only on the `-a` branch); `-a2` given but empty after parsing; `-a2` on single-end (V10); `-a2` with `--clump_only` (already rejected at `cli.rs:716`, unaffected); preset with an R2 default and no `-a2` (V5); preset without one (V6); `-a` + preset + `-a2` (V4, the workaround); auto-detected BGI, which carries an R2 default via detection rather than a flag.

**Integration.** Four user-visible changes in §7. The first is an intentional output change for affected invocations and is the reason the CHANGELOG entry needs to be explicit rather than terse.

**Remaining risks.**

- *Medium:* the restructure touches the Read 1 path incidentally, since all branches are being reshaped. V7 (md5 R1 before/after across six invocations) is the specific guard, and it is the check to run first, not last.
- *Low:* a preset's R2 default is silently lost for users who wanted it while also passing `-a2` — but that is the requested semantics, and the NOTE makes it visible.
- *Low:* the new CI case is written to pass rather than to discriminate. V15 addresses it directly.

---

## 12. Revision history

### v2 — 2026-07-29, after dual independent plan review

`PLAN_review_reviewer-A.md`, `PLAN_review_reviewer-B.md`. Both reviewers reproduced §2.2's seven-row table exactly and confirmed §2.1, §2.3, §2.4, §2.6, A2–A6. Both were read-only; the source tree was untouched.

**The Critical: v1's central safety property was unachievable.** §1, §3.5, A1 and V7 all asserted Read 1 byte-identity after the fix. Reviewer A measured that it is false on paired-end input — the joint pair-length filter drops 139 pairs instead of 4 once R2 is trimmed with the user's adapter, so `_val_1.fq.gz` changes by design while R1's own cutadapt-stage figures stay identical. Independently re-verified. v1's V7 would have failed on V1–V3 with no way for the implementer to tell "expected" from "I broke Read 1". Split into V7a (trimming invariance, all rows) / V7b (file identity, only where `-a2` does not newly apply) / V7c (reports too).

This was also the one place the reviewers **contradicted** each other: A declared A1 wrong and measured it; B judged A1 "sound as a requirement" and critiqued only the guard. A was right, and the distinction mattered — strengthening V7 without re-specifying it would have left the implementer facing a failing check.

**The design changed.** v1 restructured all seven branches of `resolve_adapter` into a single value and booked "the restructure touches Read 1 incidentally" as its only Medium risk. Reviewer A observed the risk was self-inflicted: `resolve_adapter` has exactly one caller, so applying the override one frame up is ~10 lines, touches zero existing branches, and drops the §5 signature change entirely. Adopted. Reviewer B independently proposed a named struct as the *mitigation* for the same risk — the caller-side override removes the need for mitigation, so A's shape wins.

**Corrections to v1's own claims:**

| v1 said | Actually |
|---|---|
| Read 1 output is byte-identical after the fix | False on PE — 139 vs 4 pairs dropped by the joint filter (measured) |
| `grep -n 'a2\|adapter2' ci.yml` returns nothing | Returns **8** spurious hits (`sha256sum`, the pinned toolchain SHA). Anchored form returns 0 |
| The validation job is "five invocations" (SE/PE/hardtrim5/clock/demux) | ~30 steps, 24 byte-identity assertions. `CLAUDE.md`'s summary is out of date |
| Affects "v2.0.0–v2.3.0" | There is no v2.0.0 tag; earliest v2 tag is v2.1.0. This wording was headed for a user-facing CHANGELOG |
| Auto-detect is one branch | Two — `--consider_already_trimmed` returns a synthetic suppressed preset (§3.6) |
| §6: "no measurable change", "adds no I/O" | True for the resolution code; wrong unqualified — a `-a2` over 64 bp bypasses the Myers prefilter |
| §3.5: "Empty `-a2` … still gets the default" | Ambiguous; `-a2 ""` now errors (§3.6) |
| §3.4 rationale: "should not convert working invocations into failures" | The fix does exactly that for malformed `-a2`; rationale restated on narrower ground |
| Four user-visible changes (§7) | Seven, including a new non-zero exit and a JSON-report field |
| Unit-testing the resolution matrix suffices | It cannot catch this bug class — the defect was six branches never calling the helper |

**Added:** the suppression carve-out (§3.6, user-decided); the newly-fatal input surface (§3.6); early `-a2` validation so bad input fails before the 1 M-read scan; the widened unusable-`-a2` warning with a per-mode reason (§3.4); gating `Adapter 2 (Read 2):` on `cli.paired`; V16 (the Perl-free differential oracle, now the primary check); V17–V23; the verified docs edit list (§4 Step 8).

**Two things the reviewers checked that came back clean, and are worth recording as such:** the seven-path enumeration is complete — `--poly_a` has no adapter branch in this codebase (Perl's `--polyA` does substitute `-a2`, overwriting a user value with no `unless defined` guard and applying it 5' via `-g`, so "user `-a2` always wins" is not universal even in Perl); and `--nextseq`/`--rrbs` do not touch adapter resolution. A also verified `clippy::type_complexity` would not have fired on v1's 5-tuple, so that choice was CI-safe — moot now.

**Recorded as deliberately not restored** (both reviewers found this independently): 0.6.11 keys its R2 defaults on the resolved R1 *sequence*, not the preset flag (`if ($adapter eq 'TGGAATTCTCGG')`), so Perl supplies the smallRNA/BGI R2 default even for an explicit `-a TGGAATTCTCGG`. This port does not. Out of scope for #369, but §1 claims to "restore the 0.6.11 contract", so the part not being restored is named here.

### v1 — 2026-07-29

Initial plan. Correct on the root cause, the seven-row reproduction, the not-a-2.3.0-regression finding, the Perl contract, and A4 — all independently confirmed by both reviewers.

---

## 13. Implementation notes

**Implemented 2026-07-29** on branch `fix/a2-honoured-with-presets` off `dev` @ `4276810`. All of §4 done except Step 9 (issue comment + close), which is outward-facing and awaiting go-ahead.

### What was built

| File | Change |
|---|---|
| `src/main.rs` | `apply_adapter2_override` + its call from `setup_trimming`; the NOTE; `Adapter 2 (Read 2):` display gated on `cli.paired`; redundant `-a2` parse removed from the `-a` branch; `resolve_adapter` doc line |
| `src/cli.rs` | `-a2` validated in `Cli::validate`; unusable-`-a2` warning with a per-mode reason |
| `tests/integration_adapter2.rs` | **New**, 13 tests |
| `.github/workflows/ci.yml` | `Validate -a2 with a preset and with auto-detection (#369)` — Perl md5 comparison for `--illumina` and auto-detect, plus a grep that R2 got the `-a2` sequence |
| `CHANGELOG.md` | Appended to the existing `#### Fixes` block |
| `docs/.../guide/flags.md` | `--small_rna`/`--bgiseq` `--adapter2` hedge |

**No branch's Read 1 computation changed** — the only edit inside `resolve_adapter` is the Read 2 slot of the `-a` branch (its redundant `parse_adapter_specs(&cli.adapter2)` deleted, the slot now `Vec::new()`, which is what that call already returned for an empty input). `apply_adapter2_override` takes `adapters_r1` by shared reference and cannot mutate it. That is what makes the Read 1 guarantee structural rather than tested. (v1 of these notes claimed "zero branches were touched", which both code reviewers correctly flagged as false as written.) Net: 441 → **456** tests.

### Deviations from the plan

1. **Suppression is detected from the R1 adapter, not from `DetectionResult.suppressed`.** §3.6 said `suppressed` is "already in scope" — true inside `resolve_adapter`, but *not* in the caller where the override now lives, and propagating it would have reinstated the signature change the plan had just removed. Instead the helper keys on `adapters_r1.len() == 1 && adapters_r1[0].1.is_empty()`. That is the actual semantics — suppression is the only way an R1 adapter reaches that point with an empty sequence (`adapter.rs:219-227`; `parse_adapter_specs` rejects every empty form) — and it is what `trimmer.rs` itself skips on. A label string-match would have been the fragile alternative.

2. **No `md5` dev-dependency.** The plan's V16 oracle implies hashing; there is no hashing crate in the tree and adding one to a project whose story is "no external runtime deps" is the wrong trade. The test compares decompressed bytes directly instead — same strength, and a mismatch can be inspected.

3. **One extra docs edit.** §4 Step 8 named only `flags.md:12` as needing a hedge; `--bgiseq` has the same auto-set behaviour and the same need, so the hedge covers both.

### Validation results

Baselines for V7/V13/V16 were captured **before** the first edit — V16's oracle is *today's* `-a`+`-a2` output, so it could not be reconstructed afterwards.

| Check | Result |
|---|---|
| V1/V2/V3 — `-a2` wins on all seven paths | Pass. `-a`, `--illumina`, `--nextera`, `--stranded_illumina`, `--small_rna`, `--bgiseq`, auto-detect all report `-a AAATCAAAAAAAC` for R2 |
| **V16 — differential oracle** | **Pass.** Post-fix `--illumina -a2 SEQ` *and* auto-detect `-a2 SEQ` are byte-identical to the pre-fix `-a AGATCGGAAGAGC -a2 SEQ` output, both reads |
| V4 — `-a` + `-a2` unchanged | Pass, byte-identical |
| V5/V6 — defaults intact without `-a2` | Pass, byte-identical for `--small_rna`, `--bgiseq`, `--illumina` |
| V7a — R1 per-read trimming invariant | Pass (`Reads written (passing filters): 10,000` either way) |
| V7b/V7c — R1 file + report identity where `-a2` does not newly apply | Pass |
| V8 — `--small_rna` cutoff unaffected | Pass (A4 held, as both reviewers predicted) |
| V9 — NOTE only on real displacement | Pass: fires for `--small_rna`, silent for `--illumina` |
| V10/V20 — SE warns, no contradictory R2 line | Pass |
| V11 — `-a` + preset still permissive | Pass |
| V12/V22 — `A{N}`, repeated `-a2` on preset paths | Pass (`A{10}` → `AAAAAAAAAA`; two R2 sequences) |
| **V13 — no accepted path moved** | **Pass — all 24 baseline files (FASTQ + reports) byte-identical** |
| V14 — suite + fmt + clippy | Pass: 454 tests, `fmt --check` clean, `clippy -D warnings` clean |
| **V15 — the new CI gate fails pre-fix** | **Pass.** Built `dev` in a throwaway worktree: both modes report `-a AGATCGGAAGAGC`, so the gate's grep fails. Non-vacuous |
| V17 — malformed `-a2` fatal | Pass: `ZZZQQQ`, `""`, `file:missing.fa` all exit 1 |
| V18 — suppression wins | Pass: warning fires, R2 `Reads with adapters: 0`, no `Adapter 2` line |
| V19 — specialty modes warn with a true reason | Pass for `--clock`, `--implicon`, `--hardtrim5`, `--hardtrim3`, and SE — each naming its own reason, all still exit 0 |
| V21 — poly-G piggyback intact | Pass: 1 scan line for a preset, 0 for auto-detect — unchanged from pre-fix |
| V23 — uBAM output inherits the fix | Pass |

### Iteration log

**#1 — two false negatives from my own test harness, not the code.** `zsh` does not word-split unquoted parameter expansions the way bash does, so a loop passing `m="--hardtrim5 30"` as `$m` handed clap one bogus argument `--hardtrim5 30`; V19 read as "NO WARNING" for both hardtrim modes while the direct invocation warned correctly. Re-run with arguments passed individually: all five modes warn. Same class as the earlier `mktemp` failure this session — a broken harness produces a confident wrong conclusion exactly like broken reasoning does.

**#2 — `suppression_wins_over_a2` failed on a helper bug.** The test asserted `r2_adapter(&dir) == ""`, but the suppressed adapter *is* empty, and `split_whitespace` collapses it so the helper returned the input filename. Rewritten to assert the behaviour — R2 `Reads with adapters: 0` — via a new `r2_report_line` helper, which is what actually matters and cannot be defeated by tokenisation.

**#3 — borrow errors, then a clippy lint.** `&r1()` produced temporaries dropped while borrowed; converted the fixture helpers to `const`. `cargo clippy` then flagged two `map(|s| *s)` → `.copied()`. Both mechanical; full suite re-run after each.

### Follow-ups

- **Not done: Step 9** — the comment on #369 linking the fix, and closing the issue on merge (manual; auto-close does not fire on `dev`).
- **Deviation 4 (not in the original list):** no unit test of `apply_adapter2_override`. §4 Step 5 asked for one alongside the integration matrix; `main.rs` has no `#[cfg(test)]` module and the helper is private, so the plan's own escape hatch applies — but the omission was undisclosed until the coverage audit flagged it.
- Nothing committed or pushed.
- Recorded but out of scope, per §12: Perl keys its R2 defaults off the resolved R1 *sequence* rather than the preset flag, so 0.6.11 supplies the smallRNA/BGI R2 default even for an explicit `-a TGGAATTCTCGG`. Still not restored.

---

## 14. Verification round — dual code review + coverage audit

Ran 2026-07-29 after implementation: two independent code reviewers plus a coverage audit, all in fresh contexts, all **read-only** (the #363 round had reviewers colliding on the same tree; `git status` confirmed only the five expected files modified afterwards).

`CODE_review_reviewer-A.md`, `CODE_review_reviewer-B.md`, `COVERAGE.md`. **Coverage verdict: INCOMPLETE — 2 items**, both low-severity: Step 9 (the #369 comment, deliberately deferred) and Step 5 (test coverage, PARTIAL). 53 items: 48 DONE, 1 PARTIAL, 1 MISSING, 3 DEVIATED.

### Defects found and fixed

| # | Found by | Defect |
|---|---|---|
| 1 | **A (High), B (B1), auditor (O1)** | **Single-end with a preset that has a Read 2 default printed contradictory stderr**: a WARNING that `-a2` is unused, then a NOTE that `-a2` displaced the preset default. Step 3 gated the `Adapter 2 (Read 2):` line on `cli.paired` but not the NOTE, and the override itself ran regardless. This is precisely the class V20 was written to close, re-introduced through the NOTE — and V20 could not see it, because it uses `--illumina`, which has no default to displace. Fixed by returning early from `apply_adapter2_override` when `!cli.paired`, which also stops single-end `TrimConfig` carrying a dead Read 2 adapter. Regression guard added over both `--illumina` and `--small_rna` |
| 2 | **A (Medium), B (B3), auditor (O2)** | **`A{N}` expansion announced twice**, and `1 + pairs` times on multi-pair input — because `-a2` is parsed in `Cli::validate` and again in the override, and `parse_adapter_spec` prints as a side effect. A regression on the `-a` + `-a2` path, which §3.5/V4 declare unchanged. Fixed with `parse_adapter_spec{,s}_quiet` for the validation pass, preserving §3.1's fail-before-the-scan property |
| 3 | B (B2) | **A malformed `-a2` became fatal in exactly the modes that declare it unused.** `--paired --clock -a2 ZZZQQQ` was exit 0, now exit 1 — failing on a value the tool is about to announce it will ignore. B's realistic scenario: a pipeline template reusing `-a2 file:r2_adapters.fa` across paired and single-end samples breaks on the single-end sample. Fixed by computing the unusable-reason first and skipping validation when the flag will not be used |
| 4 | B (B5) | `--help` was not hedged, only the docs page. `-a2`, `--small_rna` and `--bgiseq` doc comments now state the precedence |

### Claims corrected

- **"Zero branches of `resolve_adapter` were touched"** was false as written — the `-a` branch *was* edited. Both reviewers flagged it, and it was headed for a PR description where it does load-bearing work. Restated as "no branch's Read 1 computation changed", with the actual edit named.
- **"No change for any invocation without `-a2`"** (§3.5, §7) was false: single-end `--small_rna`/`--bgiseq` stop printing `Adapter 2 (Read 2):` whether or not `-a2` is given, because those presets set a Read 2 default regardless of pairing. B measured it against the pre-fix binary. CHANGELOG corrected.
- **Deviation 4 added:** no unit test of the helper, contrary to Step 5.

### Test coverage added

Five gaps the reviewers and auditor named, now closed: the SE-plus-`--small_rna` NOTE regression guard; stderr-content assertions on malformed `-a2` (previously a bare `!ok`, which would pass on an unrelated failure); `--implicon` and `--hardtrim3` warning rows; the R2 per-adapter breakdown V22 actually promised; and V23 (`--output-format ubam` with `-a2`), previously verified only by hand. Plus a new test that a malformed `-a2` is *tolerated* where the mode ignores it — the guard for fix 3. **15 tests in this file, 456 total.**

### What the reviewers established that I could not

- **Reviewer A ran Perl 0.6.11 for real.** `cutadapt 5.2` — the exact CI pin — plus `perl` were available on its machine, so it extracted `git show 0.6.11:trim_galore` and ran both legs of the new CI step. Both reads byte-identical in both modes. It also refuted its own hypothesis that the homopolymer-rich `-a2 AAATCAAAAAAAC` would expose alignment-DP tie-breaking differences against real cutadapt. The new gate is confirmed working, not merely plausible.
- **The auditor rebuilt the pre-fix binary and ran every `trim_galore` invocation in the `validation` job** — 25 cases extracted from `ci.yml` — under both binaries. All 25 identical. Two apparent diffs resolved to artefacts of its own harness: reports and the uBAM `@PG CL:` line embed the invoking binary's absolute path, and its pre-fix binary sat at a longer path. That converts V13 from "trust the captured baselines" to independently measured.
- **Both reviewers independently confirmed the suppression sentinel** (§13 deviation 1) by exhausting every producer of an R1 adapter in `adapter.rs`. A found it *tighter* than claimed: clap forbids `--consider_already_trimmed` with all five presets (`cli.rs:293`), so auto-detect is the only route in, and the predicate is exactly equivalent to `DetectionResult.suppressed`.

### Not adopted

- **A's suggestion to note the suppression coupling at `adapter.rs:219-227`** — worth doing, but it edits a file this change otherwise does not touch; folded into the follow-up list rather than widening the diff.
- **Trimming the `apply_adapter2_override` doc comment** per the one-line comment rule. A called it a judgement call and noted neighbouring doc comments predate the rule; the rationale sentence is the reason the early return exists, so it stays.
- **A's note that the two CI legs are md5-redundant on this fixture** (auto-detection is inconclusive on `BS-seq_10K_R1` and defaults to Illumina, so both legs produce identical output). Correct, and the second leg still earns its place as branch coverage with its own grep. Recorded here rather than as a `ci.yml` comment.
