# PLAN — `--phred64` quality-encoding fixes on the uBAM paths (issue #358)

**Issue:** [#358](https://github.com/FelixKrueger/TrimGalore/issues/358)
**Branch (proposed):** `fix/phred64-ubam`
**Base:** `dev` @ `2e3da11`
**Shape:** one PR — **a correctness constraint, not a preference** (see §5.1)
**Revision:** v2, after dual `plan-reviewer` (both APPROVE WITH REVISIONS). See §12.

---

## 1. Goal

Close two `--phred64` defects on the uBAM code paths.

| # | Combination | Current behaviour | Target |
|---|---|---|---|
| **1** | `--phred64` + `--output-format ubam`, **FASTQ** input | BAM `QUAL` stores raw `71` where the true score is `40` — off by exactly `64−33`. Inside SAM-legal `0..93`, so **no consumer errors**; quality is silently inflated. | `BamWriter` subtracts the run's actual input offset. Stored value is the true Phred score. |
| **2** | `--phred64` with uBAM **input** | `BamReader` normalises BAM raw Phred → Phred+33 ASCII; `quality.rs` then subtracts 64 from `+33` data, flooring every base to Q0. **100 % of reads discarded** when `-q > 0`. | Hard error at startup: the flag does not apply to BAM input. |

Bug 2 is not in the issue; found while reproducing Bug 1.

**These two fixes are coupled.** Fix 1 without Fix 2 introduces a *new*, silent, worse defect. See §5.1 — this is the single most important thing in this plan.

Non-goal: the FASTQ output path is untouched (§7).

---

## 2. Context

### Why FASTQ is correct and BAM is not

`FastqReader`/`FastqWriter` never interpret quality encoding. `FastqRecord.qual` is a raw ASCII `String`, sliced during trimming (`src/fastq.rs:73-135`) and written verbatim (`:61-67`). Only `quality.rs` applies an offset, and only for trim decisions. TrimGalore's FASTQ behaviour is therefore **pass-through** — Phred+64 in → Phred+64 out — matching Perl v0.6.x (which delegated to Cutadapt's `--quality-base=64`, also preserving).

BAM is categorically different: the SAM/BAM spec requires `QUAL` to hold **true Phred values** (`0..93`), not ASCII. The ASCII→raw conversion at the BAM boundary is the one place that must know the input encoding, and it hardcodes 33.

### Reproduction (release build)

```
Bug 1  --phred64 --output-format ubam, 'h'-qual FASTQ  -> stored raw 71   (correct: 40)
       control: Phred+33 'I' input, no flag            -> stored raw 40   correct
Bug 2  --phred64 -q 20 test_files/ubam_test.bam        -> 10/10 "too short", 0 written
       same input, no --phred64                        ->  0   too short, 10 written
```

`samtools view` renders `QUAL` as `ASCII(raw+33)`, so a displayed `h` (104) means stored raw 71.

### Code touchpoints

| Location | Role |
|---|---|
`src/bam.rs:60` | `const PHRED_OFFSET: u8 = 33` — **read-side only**; see §4.1 naming hazard |
`src/bam.rs:509-512` | `BamWriter` struct — gains `input_phred_offset` |
`src/bam.rs:524-543` | `BamWriter::create` — gains a parameter |
`src/bam.rs:576-579` | Comment with **two** false claims (§5 step 4) |
`src/bam.rs:580` | **Bug 1 site**: `.map(\|b\| b.saturating_sub(33))` |
`src/bam.rs:904-922` | Read side: raw → Phred+33 ASCII. Correct as-is; the reason Bug 2 exists |
`src/bam.rs:920` | Latent `u8` overflow on mixed-`0xFF` QUAL (§5 step 9) |
`src/cli.rs:144-152` | `--phred64` / `--phred33` flags |
`src/cli.rs:970-972` | `phred_offset() -> u8 { if self.phred64 { 64 } else { 33 } }` |
`src/main.rs:178-221` | "uBAM-specific format-dependent validation (T17)" — **Bug 2 guard lands here** |
`src/main.rs:237-260` | uBAM startup `NOTE:` block — O-3's target |
`src/quality.rs:32,79` | The two `saturating_sub(phred_offset)` sites — unchanged |
`src/clump_only.rs:721-726` | Documented lossless-qual invariant the guard protects (§7) |

### Where the Bug 2 guard must live

`cli.validate()` is called at `main.rs:164`; format detection not until `main.rs:181`. **`Cli::validate()` cannot see the input format**, so the guard cannot join §3.4a in `cli.rs`.

`main.rs`'s T17 block already computes `input_formats` and `any_bam`, and precedes *every* dispatch branch — specialty (`:256`), clump-only (`:349`), trim. Direct precedent: `cli.rs:544-546` records that the §3.4b rule "requires format detection and lives in main.rs." Verified by both reviewers.

### Affected write paths — 7 production, 7 test (14 total)

All uBAM output funnels through `BamWriter::write_record`, so **every** path carries Bug 1 — wider than the issue states.

| Site | Path | Status |
|---|---|---|
`src/main.rs:1797` | trim uBAM out, SE | reproduced |
`src/main.rs:1938` | trim uBAM out, PE two-file | shared writer |
`src/main.rs:2069` | trim uBAM out, PE interleaved | shared writer |
`src/clump_only.rs:766` | clump-only uBAM out, SE | shared writer |
`src/clump_only.rs:981` | clump-only uBAM out, PE | shared writer |
`src/specialty.rs:127` | `hardtrim5_to_bam` | **reproduced** (stored raw 71) |
`src/specialty.rs:175` | `hardtrim3_to_bam` | shared writer |

Test callers: `src/bam.rs:1350, 1400, 1428, 1483, 1522, 1646` and **`src/clump_only.rs:1531`** (a test caller of `clump_only_single_to_bam`, which gains a parameter). All compiler-caught.

---

## 3. Behavior

### 3.1 Bug 2 guard — reject `--phred64` with BAM input

**Implement this first.** See §5.1.

1. After `any_bam` is computed in the T17 block, `if cli.phred64 && any_bam` → `anyhow::bail!`.
2. Applies **regardless of output format and mode**.
3. Error text: the flag is meaningless for BAM input, currently corrupts results, remedy is to drop it.

**Justification — corrected.** v1 justified this as "the defect is on the read side." That is true for the trimming pipeline but **false for modes that do no quality arithmetic**. Verified: `--hardtrim5 20 --phred64 <bam>` and `--clump_only --phred64 <bam>` both succeed today with the flag entirely inert. The real justification is:

> `--phred64` declares the ASCII encoding of the input. BAM has no ASCII encoding to declare — it stores raw Phred by spec. The flag is therefore meaningless for BAM input in **every** mode. A uniform rule is safer to document, test, and refactor against than a mode-dependent one ("rejected unless you're in hardtrim, or clump-only-to-FASTQ, or…"), which is one refactor away from being wrong.

**Reconciling with O-2.** §10 O-2 declines extra rejection for `--clump_only` + FASTQ on the grounds that the flag is inert. §3.1 rejects `--clump_only` + BAM where it is *equally* inert. The governing principle, stated once:

> Reject where the flag is meaningless **because of the input format**. Leave inert-but-meaningful FASTQ cases alone — FASTQ genuinely *has* an ASCII encoding to declare, so the flag is well-formed there even when a given mode ignores it. BAM does not.

**Edge cases**
- Mixed FASTQ + BAM input → `any_bam` true → rejected. Correct; the offset cannot be per-input.
- `--phred64` + FASTQ-only → unaffected.
- BAM input without `--phred64` → unaffected.
- `--phred33 --phred64` together with BAM input → still rejected (`self.phred64` is true). Benign.

### 3.2 Bug 1 fix — honour the offset in `BamWriter`

1. `BamWriter::create` accepts `input_phred_offset: u8`, stored on the struct.
2. `write_record` subtracts `self.input_phred_offset` instead of the literal `33`.
3. All 7 production sites pass `cli.phred_offset()`; all 7 test sites pass `33`.

**Edge cases**
- `saturating_sub` is **retained, and is load-bearing rather than defensive.** The existing comment claims `FastqReader` rejects sub-33 input; it does not — there is no quality-range validation anywhere in `fastq.rs`. Reviewer B demonstrated an ASCII-32 byte reaching the writer and saturating. So it guards real malformed input.
- Post-fix, a Phred+33 file mistakenly run with `--phred64` floors to raw 0. This is the better failure direction than today's silent inflation into the legal window — but "visible" is optimistic: an all-`!` BAM is consumable as legitimate Q0. **O-3's startup `NOTE:` is the mitigation**, and that is its rationale.
- **Genuine Solexa / Illumina 1.0 data also floors.** True Solexa files contain ASCII 59–63 (negative Solexa scores), and `--phred64` is exactly the flag those users reach for. Post-fix such bytes floor to 0 rather than fabricating Q26. Still strictly better than today and identical to `quality.rs`'s existing behaviour, so no change — but it is the more plausible occurrence and belongs in the record.
- Missing-qual BAM records: `BamReader` synthesises `'!'` (33). Cannot co-occur with `--phred64` — §3.1 rejects that combination.
- Zero-length qual: reachable (a read fully consumed by adapter trimming under `--length 0`); `noodles` emits `*`/`*`, no panic. Confirmed by Reviewer B.
- **`debug_assert_eq!` at `bam.rs:581` is compiled out in CI.** `[profile.release]` in `Cargo.toml` sets no `debug-assertions`, and `.github/workflows/ci.yml:74` runs `cargo test --release`. It is **not** edge-case coverage and must not be cited as such. The fix cannot change lengths, so no action beyond not relying on it.

### 3.3 Out of scope, documented

`--phred33` / `--phred64` carry no clap `conflicts_with`, and `phred_offset()` tests only `self.phred64`, so both together silently yields 64. Pre-existing, not uBAM-specific. §10 O-1.

---

## 4. Signature

```rust
/// Open `path` for uBAM writing.
///
/// - `source_header`: input uBAM's SAM header (propagated as-is plus a
///   trim_galore `@PG` line). Pass `None` for FASTQ input — a minimal
///   `@HD VN:1.6` + trim_galore `@PG` header is synthesised.
/// - `_preserve_tags`: kept for API symmetry with `BamReader::open_*_with_tags`;
///   not consulted (the textual tail in `FastqRecord.id` is authoritative).
/// - `command_line`: verbatim invocation, for `@PG CL:`.
///   **PLAN §4 deviation:** <preserve the existing paragraph here verbatim —
///   it explains why `command_line` is a 4th parameter. Do not drop it.>
/// - `input_phred_offset`: ASCII offset of the *incoming* `FastqRecord.qual`
///   (33 or 64), i.e. `cli.phred_offset()`. BAM `QUAL` must hold true Phred
///   scores per the SAM spec, so this is subtracted per byte on write.
///   **Distinct from the module-level `PHRED_OFFSET` const** (§4.1).
pub fn create<P: AsRef<Path>>(
    path: P,
    source_header: Option<&Header>,
    _preserve_tags: &[String],
    command_line: &str,
    input_phred_offset: u8,
) -> Result<Self>
```

**Required parameter, not a defaulted setter.** A `with_phred_offset()` builder avoids touching 14 call sites but silently defaults any forgotten site to the buggy behaviour. A required parameter makes the compiler enumerate every site — the property we want, since the bug being fixed is precisely "a site that never considered it." CI is `-D warnings`; no partial-adoption escape hatch.

### 4.1 Naming hazard — `PHRED_OFFSET` already exists

`src/bam.rs:60` defines `const PHRED_OFFSET: u8 = 33`, used only at `:920` on the **read** side. Post-fix the module would carry two near-identical names with *different semantics*:

| Name | Value | Meaning |
|---|---|---|
`PHRED_OFFSET` (const) | fixed 33 | offset the **reader adds**; must stay fixed |
`BamWriter.input_phred_offset` | 33 or 64 | offset the **writer subtracts** |

"Conflating two offsets" is the exact bug class being fixed, so: name the field **`input_phred_offset`** (not `phred_offset`), and add a one-line comment at the const marking it read-side-only and unrelated to the writer field. Assumption 7's rationale must acknowledge a named constant already exists in this module.

---

## 5. Implementation outline

### 5.1 Ordering is load-bearing — guard first

> **v1 said step 8 was "independent of 1–7 and could land first."** That was wrong, and the error is dangerous. Both reviewers independently proved it by execution.

`--clump_only` never interprets quality, and `BamReader` normalises to Phred+33 — so today's hardcoded `saturating_sub(33)` is **correct** on BAM-input paths even with `--phred64` set. Verified on `dev`:

```
$ trim_galore --clump_only --phred64 --output-format ubam ubam_test.bam
$ samtools view ubam_test.bam        | cut -f11 | sort | md5   ->  31878eb2…
$ samtools view out/…_clumped.bam    | cut -f11 | sort | md5   ->  31878eb2…   identical
```

Apply the writer fix **without** the guard and it subtracts 64 from the reader's Phred+33 bytes → saturates → **an all-Q0 BAM, every read retained, zero exit status, no warning.** Strictly worse than either bug: Bug 2 fails *loudly* (100 % discarded); this fails *silently*.

Consequences:

1. **Implement the guard (step 1) before the writer fix (steps 2-8).** No compile dependency either way, so the tree builds at each step regardless — but no intermediate commit may silently corrupt a working path.
2. **This is why it is one PR.** Not a shape preference: shipping the writer fix alone would introduce a new silent-corruption bug. State this in the PR body; it pre-empts any request to split.
3. **The guard can never be relaxed to warn-and-ignore** without re-opening this path.
4. **Behaviour regression to disclose.** These work correctly today and become hard errors:
   - `--clump_only --phred64 <bam> [--output-format ubam]`
   - `--hardtrim5/3 N --phred64 <bam>`
   The rejection is still right (§3.1), but it is the only place the change removes working behaviour, so it needs a CHANGELOG *behaviour-change* line, not just a fix line.

### 5.2 Steps

1. **`src/main.rs`, T17 block (~`:190`, after `any_bam`)** — add the guard:
   ```rust
   if cli.phred64 && any_bam {
       anyhow::bail!(
           "--phred64 cannot be used with unaligned BAM input. BAM stores raw \
            Phred scores directly, so there is no ASCII encoding to declare — \
            the reader always yields Phred+33 internally. Passing --phred64 \
            here subtracts 64 from Phred+33 data, discarding every read as \
            low-quality (or, with --clump_only, silently zeroing quality). \
            Drop --phred64 for BAM input."
       );
   }
   ```
2. **`src/bam.rs:509`** — add `input_phred_offset: u8` to `struct BamWriter`.
3. **`src/bam.rs:524-543`** — add the parameter to `create`, store it, extend the doc comment per §4 (**keeping** the existing "PLAN §4 deviation" paragraph).
4. **`src/bam.rs:576-580`** — the fix plus a comment rewrite that corrects **both** falsehoods:
   ```rust
   // BAM stores raw Phred bytes; `FastqRecord.qual` is ASCII at the input's
   // declared offset (33, or 64 under --phred64). `saturating_sub` is
   // load-bearing, not merely defensive: nothing in fastq.rs validates the
   // quality range, so a malformed sub-offset byte can reach here and must
   // floor to 0 rather than wrap.
   let raw_qual: Vec<u8> = record.qual.bytes()
       .map(|b| b.saturating_sub(self.input_phred_offset))
       .collect();
   ```
   Do **not** carry forward "FastqRecord.qual is Sanger ASCII (+33)" (the bug's premise) or "FastqReader rejects sub-33 input" (also false).
5. **`src/bam.rs:60`** — one-line comment marking `PHRED_OFFSET` read-side-only (§4.1).
6. **`src/main.rs:1797, 1938, 2069`** — pass `cli.phred_offset()` (all three enclosing fns already take `cli: &Cli`).
7. **`src/specialty.rs:127, 175`** — add `phred_offset: u8` to `hardtrim5_to_bam` / `hardtrim3_to_bam` (6 → 7 params; **`specialty.rs` carries no `#[allow(clippy::too_many_arguments)]` and needs none** — 7 is at clippy's threshold, not over it). Update callers at `main.rs:275, 299`.
8. **`src/clump_only.rs:766, 981`** — same for `clump_only_single_to_bam` / `clump_only_paired_to_bam_one_pair` (both already carry the clippy allow). Update callers at `main.rs:460, 526, 557` and the test caller at `clump_only.rs:1531`.
9. **`src/bam.rs:920`** — fold in the neighbouring latent overflow: `(b + PHRED_OFFSET)` → `b.saturating_add(PHRED_OFFSET)`. The all-`0xFF` sentinel check at `:908` only catches records where *every* byte is `0xFF`; a **mixed** QUAL overflows `u8` — panic in debug, wrap to `' '` in release. One character, zero risk, symmetric with the write side, same defect family. (Flagged by Reviewer A as out-of-scope; folding in rather than filing, because leaving a known overflow two lines from a site being touched invites the same mistake.)
10. **`src/main.rs:237-260`** — add the O-3 `NOTE:` to the existing uBAM startup block, worded to cover both the input encoding and the output representation (§10 O-3).
11. **`CHANGELOG.md`** — under `### Unreleased` → `#### Fixes` (**the heading already exists**; append). Plus a distinct behaviour-change line per §5.1 item 4.
12. **Docs** — primary target `docs/src/content/docs/guide/quality.md` §"Phred encoding" (~`:29-31`), where `--phred64` is documented and where a user reaching for the flag will look. Optionally a one-line cross-reference in `outputs.md`, but **not** inside its "Rejected at CLI-validate time" list — the guard is deliberately *not* at validate time, so filing it there would state something false. Plain register per the docs convention.

---

## 6. Efficiency

No measurable impact. A field load replaces an immediate, hoisted out of the loop by the optimiser, on a path dominated by BGZF compression. No allocation change. The guard is two boolean tests once per run on already-computed values.

---

## 7. Integration

**Reads:** `cli.phred64`, `cli.phred_offset()`, `input_formats` / `any_bam`.
**Writes:** BAM `QUAL` bytes only. No header or report change.

**The guard preserves `clump_only`'s documented lossless invariant.** `src/clump_only.rs:721-726` states that every input record appears in the output with `R.id`, `R.seq`, `R.qual` **byte-identical**. Under `--phred64` + BAM input, the writer fix alone would zero `R.qual` (§5.1) — a direct breach of the mode's headline property. §3.1's guard is what keeps it intact. Record this in the PR body; a reviewer of #356 will look for it.

**FastQC on uBAM output is an unlisted beneficiary.** `--fastqc` with `--output-format ubam` runs fastqc-rust directly on the `.bam`, and fastqc-rust reads `QUAL` as raw Phred. Pre-fix, `--phred64 --output-format ubam --fastqc` produces per-base-quality plots centred on **Q71**; post-fix they are correct. Worth a validation row and a CHANGELOG clause — a user-visible artefact the unit test cannot reach.

**Byte-identity invariants — untouched.** `bam.rs` holds the codebase's only ASCII↔raw conversions (`:580` write, `:920` read); `fastq.rs` performs no offset arithmetic at all. The FASTQ writer is not modified, so the v0.6.11 `validation` matrix is out of scope by construction. State this explicitly in the PR body — it is the invariant most likely to draw scrutiny on a quality change.

**Backward compatibility — no golden-fixture churn.** The load-bearing check is a grep, and it is self-sufficient:

> `grep -rn "phred64" tests/ .github/workflows/` → **zero hits.** Nothing in the test suite or CI passes the flag, so `phred_offset()` returns 33 everywhere and an offset-conditional change is provably a no-op for every existing test and fixture.

The uBAM golden fixtures (`ubam_out_se_REFERENCE.bam`, `ubam_out_pe_REFERENCE.bam`) are generated without `--phred64` per the regen recipe in `test_files/README.md`, so they cannot move. Note the README lists "a bug fix in seq/qual handling" as a regen trigger — this qualifies textually but not materially; pre-empt the question in the PR body.

> **v1 also claimed "every fixture in `test_files/` is Phred+33," verified by minimum quality byte. That claim is false and is withdrawn.** `truncated.fq.gz` is Phred+64: min byte 66, max 99 → Q2..Q35 at offset 64 (a textbook Illumina 1.5 profile) versus Q33..Q66 at offset 33, and Q66 does not occur in Illumina data. Its `HWUSI-EAS611_0001:…#0/1` read names are Illumina GA/GAII vintage, and the trailing `B` runs are the Illumina 1.5 read-segment quality-control indicator (`B` = 66 = Q2 at offset 64). It is an orphan fixture — `grep -rn "truncated\.fq" src/ tests/ .github/` returns nothing — so there is no functional impact, but the heuristic was wrong and is removed rather than repaired. **Consequence for §9:** any encoding note added to `test_files/README.md` must not assert all fixtures are Phred+33; if it records encodings at all, `truncated.fq.gz` must be listed as the Phred+64 exception.

**Downstream:** consumers of uBAM output (`samtools fastq`, Bismark's uBAM reader and its cross-tool byte-identity CI) currently receive inflated `QUAL` for Phred+64 runs; post-fix they receive correct values. Phred+33 runs see no change.

---

## 8. Assumptions

**Fixed rules**
1. SAM/BAM `QUAL` stores true Phred `0..93`, not ASCII. Basis for the whole fix.
2. `BamReader` normalises to Phred+33 on read (`bam.rs:904-922`) and stays that way. The guard depends on it.
3. TrimGalore's FASTQ output is encoding pass-through, and that contract is preserved. Verified by absence of offset arithmetic in `fastq.rs`.
4. `--phred64` describes the **input** encoding, never a requested output encoding.

**Configurable / chosen**
5. `--phred64` + BAM input is a **hard error**, not warn-and-ignore. *Post-fix* rationale (stronger than v1's): warn-and-ignore would silently destroy **quality** rather than loudly destroy **reads**. Pre-fix it discards 100 % of reads; post-fix, without the guard, it writes an all-Q0 BAM with a zero exit status. Neither is acceptable, and the second is worse.
6. `input_phred_offset` is a **required** `create` parameter (§4).
7. Test call sites pass `33` literally rather than reusing a constant — keeps existing assertions verbatim and makes the Phred+33 assumption visible per site. Note `bam.rs:60` *does* define `PHRED_OFFSET = 33`, but it is read-side-only with different semantics (§4.1), so reusing it here would be actively misleading.

**~~Unverified~~ RESOLVED**
8. ~~The `clump_only.rs` / `specialty.rs` sites can reach a `phred_offset` without restructuring.~~ **Verified by both reviewers.** All four create the writer directly inside a function called from `main()`; `main.rs`'s three sites already take `cli: &Cli`. One added parameter each, **zero intermediate functions to thread through.**

---

## 9. Validation

Checks 1, 3, 4, 5, 6, 10, 11 are **integration tests**, not manual steps — the bug class is faulty *wiring*, and a unit test pins only arithmetic.

| # | What | Where | Expected | Fails pre-fix? |
|---|---|---|---|---|
1 | Bug 1 fixed, trim SE (`main.rs:1797`) | `tests/integration_ubam_out.rs` | raw qual == `[40; 34]` on the Phred+64 **FASTQ** fixture | **yes** (raw 71) |
2 | Phred+33 unchanged | same | raw 40 from `'I'` | no — **over-correction guard**, not a bug guard |
3 | Bug 1 fixed, specialty (`specialty.rs:127`) | same | `--hardtrim5 20 --phred64 --output-format ubam` on the Phred+64 **FASTQ** fixture → raw 40 | **yes** (raw 71) |
4 | Bug 1 fixed, clump-only (`clump_only.rs:766`) | `tests/integration_clump_only_ubam.rs` | `--clump_only --phred64 --output-format ubam` on the Phred+64 **FASTQ** fixture → raw 40 | **yes** |
5 | Guard: BAM in, FASTQ out | `tests/integration_ubam_out.rs` | `--phred64 -q 20 ubam_test.bam` → non-zero, stderr names `--phred64`, **no output dir created** (`ensure_output_dir` is at `main.rs:235`, after the guard) | **yes** |
6 | Guard: BAM in, **uBAM out** | same | `--phred64 --output-format ubam ubam_test.bam` → rejected | **yes** — the silent all-Q0 path |
7 | Guard: mixed input | same | `--phred64 --paired r1.fq.gz r2.bam` → rejected | yes |
8 | Guard: clump-only + BAM + uBAM out | `tests/integration_clump_only_ubam.rs` | rejected | **yes** — the other silent all-Q0 path |
9 | Guard: **early-returning mode** | `tests/integration_ubam_out.rs` | `--hardtrim5 20 --phred64 ubam_test.bam` → rejected | yes |
10 | Unit: offset arithmetic | `src/bam.rs` tests | `BamWriter` with offset 64, `'h'` in → stored raw **40**; offset-33 counterpart → 40 from `'I'` | **cannot compile** pre-fix (pins the arithmetic) |
11 | FastQC on uBAM | `tests/integration_ubam_out.rs` or manual | `--phred64 --output-format ubam --fastqc` → `fastqc_data.txt` mean quality ≈ 40, not ≈ 71 | **yes** |
12 | No regression | `cargo test`; `clippy --all-targets --release -- -D warnings`; `fmt --check` | all pass; **no test removed** (baseline **409** `#[test]`, not 407) | — |
13 | Golden fixtures unchanged | uBAM `assert_ubam_eq` references | pass unmodified — proves §7's no-op claim | no — by design |

**Check 9 rationale:** `main.rs:237-240` records code-review finding B-NIT-2 — the uBAM startup NOTEs "never fired on the hardtrim BAM paths (which return earlier)." This exact class of mis-placement has bitten this file before, so the guard's coverage of early-returning modes gets its own test.

**Unit-test read-back layer (check 10):** inspect at the raw BAM layer via `bam::io::Reader` + `Record::quality_scores()`, following the in-file precedent at `src/bam.rs:1522` (`bam_writer_aux_typed_int_and_float_round_trip`), which deliberately avoids the reader's `+33` normalisation. A `BamReader::open` round-trip also detects the bug but asserts the property indirectly.

**Fixture:** commit `test_files/phred64_test.fastq` — 4 reads × 34 bp, qual all `'h'`. Verified to run unaided (adapter auto-detect falls back to Illumina, 0 found, 4/4 written). No integration test in this repo generates inputs inline, so a committed fixture matches convention. Document it in `test_files/README.md` per that file's per-fixture provenance convention — subject to the §7 caveat about not asserting a false blanket encoding claim.

**CI:** optional. `samtools` is available in `validation-ubam` and (since #356) `validation`. The integration tests above cover the wiring more precisely than shell assertions. If a CI step is added, prefer `validation-ubam` and use the tempfile-split idiom rather than `unzip -l | grep -q`-style pipes, per the pipefail SIGPIPE hazard fixed in #357.

---

## 10. Questions or ambiguities

**Settled by the user before v1**
- Fix approach → **split fix**. (Rejected: reject-both; normalise-on-read — the latter would change FASTQ output from pass-through to converted, breaking Perl v0.6.11 parity on a path no CI covers.)
- PR shape → **one PR**. v2 note: this is now understood as a *correctness constraint* (§5.1), not a preference.

**Open (documented, not blocking)**
- **O-1.** `--phred33 --phred64` together silently yields 64 (no clap `conflicts_with`; `phred_offset()` tests only `phred64`). Pre-existing, not uBAM-specific; fixing it changes FASTQ-path behaviour. Interaction with the new guard is benign. Assumption: leave as-is.
- **O-2.** No extra rejection for `--clump_only` + **FASTQ**. Governing principle now stated explicitly in §3.1: reject where the input *format* makes the flag meaningless (BAM), not merely where a mode happens to ignore it (FASTQ). v1 held two contradictory principles here; resolved.
- **O-3. Adopted.** Emit a startup `NOTE:` for `--phred64` + uBAM output in the existing block at `main.rs:237-260`. Rationale is *not* stylistic consistency: it mitigates the all-Q0 silent-floor case in §3.2, which is otherwise indistinguishable from legitimate Q0 data downstream. Wording should cover both halves at once, e.g.:
  > `NOTE: input is Phred+64; output BAM QUAL stores true Phred scores (0–93) per the SAM spec, not ASCII.`

  This also resolves a reporting ambiguity Reviewer B found: the trimming report prints "Quality encoding type selected: ASCII+64" beside a BAM holding raw Phred. That correctly describes the *input* and is not a bug, but a user reading both could conclude the BAM is ASCII+64.

---

## 11. Self-Review

**Efficiency.** No complexity or allocation change (§6).

**Logic.** Traced the full quality lifecycle both directions. The two bugs are genuinely distinct — Bug 1 write-side (offset never applied), Bug 2 read-side (offset applied to already-normalised data) — so no single change fixes both. v2's key correction is that they are also **causally coupled**: fixing 1 without 2 opens a silent all-Q0 path.

**Round-trip integrity — improved, not compromised.** Post-fix, `--phred64 <fastq> --output-format ubam` writes true Phred; reading it back needs no flag (the reader normalises). There is no legal "Phred+64 BAM" — the spec forbids it — so refusing `--phred64` on BAM input creates **no** write/read asymmetry. The only BAMs the fix cannot re-read correctly are *pre-fix TrimGalore outputs*, which `--phred64` could never have fixed anyway (they need −31, not −64).

**Adjusted in v2 (from dual review):**
- **Inverted the implementation order** and rewrote the rationale — v1's "step 8 is independent" was wrong and would have introduced a worse, silent bug. Both reviewers proved it independently; I re-verified (byte-identical QUAL md5 today).
- **Withdrew the fixture-encoding heuristic.** `truncated.fq.gz` is Phred+64 (min 66 / max 99, Illumina 1.5 `B`-run indicator, GA/GAII read names). v1 also miscounted the ambiguous set — named six files, only three have min 73. Conclusion survives on the grep check alone.
- **Corrected §3.1's justification.** "The defect is read-side" is false for hardtrim and clump-only, where the flag is inert and the invocations work today. Replaced with the input-format argument, and reconciled the contradiction with O-2.
- **Promoted validation from manual to integration tests** (7 of 13), and added the two silent-all-Q0 combinations, the early-returning-mode case, and the FastQC row. Labelled which checks actually fail pre-fix — checks 2 and 13 do not, and are over-correction/no-op guards rather than bug guards.
- **Named the second falsehood** in the comment being rewritten ("FastqReader rejects sub-33 input"), and reframed `saturating_sub` as load-bearing rather than defensive.
- **Renamed the field `input_phred_offset`** after finding `PHRED_OFFSET` already exists in the same module with different semantics (§4.1).
- **Stopped citing `debug_assert_eq!`** as coverage — compiled out under `cargo test --release`, which is what CI runs.
- **Retargeted the docs** to `guide/quality.md`; `outputs.md`'s exclusion list is headed "Rejected at CLI-validate time", which the guard deliberately is not.
- **Folded in `bam.rs:920`'s `saturating_add`** rather than filing it — a real `u8` overflow on mixed-`0xFF` QUAL, two lines from a site being touched, same defect family.
- Marked Assumption 8 resolved; corrected counts (14 call sites, 409 tests) and line refs (`cli.rs:970-972`, `main.rs:237-260`); noted `specialty.rs` needs no clippy allow; noted `#### Fixes` already exists.
- Preserved the existing "PLAN §4 deviation" doc paragraph, which v1's replacement block silently dropped.

**Remaining risks**
- **Low.** The new Phred+64 fixture must be documented without repeating v1's false blanket claim (§7).
- **Informational.** Two currently-working invocations become hard errors (§5.1 item 4) — needs a CHANGELOG behaviour-change line, not just a fix line.
- **Informational.** Archived pre-fix `--phred64` uBAM output **cannot be repaired by re-running TrimGalore**: it stores `true + 31`, so reading it back yields Q71, and `--phred64` would subtract 64 rather than 31 — and is now rejected for BAM input anyway. Users must regenerate **from the original FASTQ**. State this explicitly in the CHANGELOG.

---

## 12. Revision history

**v2** (this document) — after dual `plan-reviewer`, both **APPROVE WITH REVISIONS**. Reports at `PLAN_review_reviewer-{A,B}.md`. All 3 critical and all 8 important items folded in, plus the worthwhile optional items. Every reviewer factual claim independently re-verified before adoption — including one that corrected Reviewer B (`truncated.fq.gz` orphan status: B's "zero references" is right for the fixture; a broader `grep truncated` matches unrelated error-message text).

**v1** — initial plan from investigation. Two material errors, both caught: the implementation ordering (would have introduced a silent all-Q0 corruption path) and the fixture-encoding claim (false).

---

## 13. Implementation notes

Implemented on `fix/phred64-ubam` from `dev` @ `2e3da11`. **Not yet pushed** — awaiting the dual `code-reviewer` + `plan-manager` verification gate.

### Order followed

§5.1's inverted order was honoured. The guard (step 1) was implemented and **verified in isolation before the writer was touched**: all four BAM-input routes (FASTQ-out, uBAM-out, `--clump_only`, `--hardtrim5`) rejected, and no output directory created. Only then were steps 2-12 applied. At no point did the tree contain the writer fix without the guard.

### Call-site count: 12 production edit points, not 7

The plan counted 7 `BamWriter::create` sites. Threading a parameter also required updating the **5 wrapper-function callers** in `main.rs` (`:306`, `:330`, `:491`, `:560`, `:588`) for `hardtrim{5,3}_to_bam` and `clump_only_*_to_bam`. The plan anticipated this in step 7/8 ("Update callers at…") but the §2 headline figure of 7 undercounts the edit surface. Test callers: 6 in `bam.rs` + 1 in `clump_only.rs` = 7. Total 19 edit points.

Every one was compiler-enumerated, which is exactly the property the required-parameter choice (§4) was made to obtain.

### Deviation: fixture composition changed after an empirical finding

The plan specified a fixture of "4 reads × 34 bp, qual all `'h'`". **Implemented as 24 × `'h'` (Q40) + 10 × `'B'` (Q2)** instead, after validation check 11 (FastQC on uBAM output) failed in an instructive way.

With the uniform-`'h'` fixture, post-fix FastQC reported per-base mean **9.0**, not the ≈40 the plan predicted:

| Fixture | FastQC `Encoding` | Reported mean |
|---|---|---|
| uniform `'h'` (all Q40) | `Illumina 1.5` — wrong | 9.0 |
| 24×`'h'` + 10×`'B'` | `Sanger / Illumina 1.9` — right | 40.0 / 2.0 |

Cause: `fastqc-rust` infers encoding from the **minimum** quality byte. A correctly-written all-Q40 BAM reads back as ASCII 73 throughout, leaving no byte below 64, so the heuristic guesses Phred+64 and subtracts 64 → Q9. The BAM was spec-correct; the *fixture* was degenerate.

This is the same min-byte-heuristic weakness §7 documents when withdrawing v1's fixture-encoding claim — encountered again from inside a dependency. Check 11's stated expectation ("mean ≈ 40, not ≈ 71") was therefore wrong as written: pre-fix it read ≈40 *by accident* (two errors cancelling, exactly as Reviewer A described in #358), and post-fix it reads 9 on a uniform fixture.

The mixed-quality fixture resolves it and is strictly better:
- Mirrors real Illumina 1.5 `B`-run tails (`truncated.fq.gz` has the same shape).
- Exercises **two** quality values, so the tests verify the offset is applied uniformly rather than pinning a single number.
- Restores FastQC's encoding detection, making check 11 meaningful.

Tests pass `-q 0 --length 0` so the low-quality tail survives to the writer; without it the default `-q 20` trims the `B`-run away and only one offset is exercised. Recorded in `test_files/README.md`.

### Deviation: `bam.rs:920` residual documented, not fully fixed

Step 9 applied `saturating_add` as planned, removing the `u8` overflow. A residual remains and is commented in place: a mixed-sentinel byte saturates to 255 rather than being mapped to `QUAL_MISSING_REPLACEMENT` per byte the way the all-`0xFF` branch does. Mapping it would be a semantic change to sentinel handling, beyond this fix's remit. Well-formed BAM does not mix sentinel and real scores.

### Iteration log

1. Guard added and verified standalone across all four BAM-input routes; no output created. Ordering constraint satisfied.
2. `input_phred_offset` threaded. `cargo check --all-targets` used iteratively to enumerate call sites; two `main.rs` sites were missed by a first regex pass (one ends `.with_context(…)` rather than `)?;`, one is a lib test) and were caught by the compiler.
3. Tests written; 418 passing. FastQC check 11 then reported mean 9.0 — investigated rather than accepted, root-caused to fixture degeneracy plus `fastqc-rust`'s min-byte heuristic (see deviation above).
4. Fixture rebuilt with mixed quality; assertions strengthened to a two-value expectation. FastQC detection correct; 418 still passing.

### Verification

| Gate | Result |
|---|---|
`cargo fmt --all -- --check` | clean |
`cargo clippy --all-targets --release -- -D warnings` | clean |
`cargo test` | **418 passing, 0 failed** (409 baseline + 9 new: 1 unit, 6 uBAM-out, 2 clump-only) |
Bug 1, trim SE | stored raw `[40×24, 2×10]`; pre-fix `[71, 33]` |
Bug 1, `--hardtrim5` | stored raw 40 across the kept 20 bp; pre-fix 71 |
Bug 1, `--clump_only` | stored raw `[40×24, 2×10]` |
Bug 2 guard | rejected on all four routes; no output dir created |
Phred+33 control | unchanged (`illumina_10K` → `BBBBBFFFFFFF`) |
FastQC on uBAM | `Sanger / Illumina 1.9`, base 1 mean 40.0, base 30 mean 2.0 |
Golden uBAM fixtures | `assert_ubam_eq` references pass unmodified — §7's no-op claim holds |

Diff: 9 files changed, +500/−14, plus the new `test_files/phred64_test.fastq`.

---

## 14. Post-review revisions (dual `code-reviewer` + `plan-manager`)

**Verdicts:** Reviewer A APPROVE WITH REVISIONS · Reviewer B APPROVE WITH REVISIONS · plan-manager **COMPLETE** (41 DONE / 1 PARTIAL / 0 MISSING). Reports: `CODE_review_reviewer-{A,B}.md`, `PLAN_manager_report.md`.

Both reviewers independently re-derived the §13 fixture deviation and cleared it. Reviewer B went further than §13 did and established that the FastQC misdetection also occurs on the **FASTQ output path** — untouched by this change — and that the lowest-quality-char rule is inherited FastQC 0.12.1 behaviour rather than a `fastqc-rust` invention. So the deviation is genuine fixture degeneracy, not a fix defect.

Every reviewer claim below was independently reproduced before being acted on.

### One real regression, introduced by step 9 and now fixed

`saturating_add` removed the `u8` overflow but replaced it with something worse: `255u8 as char` is U+00FF, **two bytes in UTF-8**, inside a `String` whose byte length the codebase assumes equals `seq.len()`. Reproduced:

```
saturating_add : chars=4 bytes=5  [73,73,195,191,73]   is_char_boundary(3)=false
wrapping (old) : chars=4 bytes=4  [73,73,32,73]        invariant intact
min(93)+33     : chars=4 bytes=4  [73,73,126,73]       invariant intact
```

`FastqRecord::truncate` / `clip_5prime` / `clip_3prime` index that `String` by byte offset, so the change traded a debug-only arithmetic panic for a release-mode char-boundary panic — or a silently wrong score count — in the exact input class it was added to harden. The in-code comment asserting the residual was benign was the untrue part.

**Resolved rather than documented.** Reviewer A's framing is right that mapping the sentinel per byte is not a semantic change but makes the mixed case *agree* with the all-`0xFF` branch two lines above. Sentinel → `QUAL_MISSING_REPLACEMENT`; anything else out of spec clamped to a new `MAX_PHRED_SCORE = 93`. Every output byte is single-byte ASCII, so the one-byte-per-base invariant now holds unconditionally. Pinned by `qual_ascii_conversion_is_single_byte_per_base`, which asserts byte length, ASCII-ness, and that byte-indexed slicing does not panic.

### Three false user-facing claims, all mine, all corrected

| Claim | Reality |
|---|---|
| CHANGELOG: `--fastqc` plots "were centred on Q71" | **False.** Pre-fix FastQC was *accidentally correct* — the writer's 31-point inflation and the heuristic's wrong offset cancelled exactly, giving mean 40.0. The `Encoding` label was wrong, not the plot. §13 had already recorded this; the CHANGELOG carried the withdrawn §7 prediction forward. |
| CHANGELOG: archived output "cannot be repaired by re-running" | **False, and harmful advice.** Verified recoverable in two commands with no access to the original FASTQ: `samtools fastq -n` reproduces the original Phred+64 ASCII, and re-running with `--phred64` gives byte-identical output (md5 `28f55edc…` both ways). Only the *single-invocation* repair is blocked. |
| "all-zero-quality BAM", in four places incl. a runtime `NOTE:` | **False.** The arithmetic is `q − 31` floored, not zeroed. Measured on `clock_10K_R1.fastq.gz`: **97.0% of bases non-zero** (Q7 × 604 444, Q3 × 48 659, Q6 × 24 729, Q0 × 21 565 = 3.0%). Worse, it *disarmed* the mitigation it justified — O-3's NOTE handed users a detection heuristic ("look for all-zero quality") that would not fire. The real artefact reads as plausible poor-quality data, which is harder to spot. |

Also corrected: the CHANGELOG's list of removed-working-invocations was presented as exhaustive at two but `-q 0` on the ordinary trim path is a third (`quality_trim_3prime` with `cutoff = 0` trims nothing), and it is the mode most users are in. Generalised to "any BAM-input run in which quality trimming was inert".

### Prose contradictions fixed

- `test_files/README.md` claimed the min-byte heuristic *misclassifies* `truncated.fq.gz`. It classifies it **correctly** (min 66 ≥ 64 → Phred+64). What fails is assuming Phred+33 without looking. The heuristic's real weakness runs the other way and `phred64_test.fastq` demonstrates it — in a file whose whole purpose is preventing this class of error.
- `guide/quality.md` said both "every read discarded" and "results are unchanged"; those cannot both hold. Now scoped by mode, matching how the CHANGELOG already scoped the same sentence.
- Added the `Phred+64 → uBAM → FASTQ` round-trip note: output is Phred+33 by design, since the BAM boundary discards the ASCII encoding. Correct behaviour, but otherwise filable as a bug.

### Coverage and maintainability

- **PE two-file writer now tested** (`phred64_ubam_out_pe_two_file_stores_true_phred`, + `phred64_test_R{1,2}.fastq`). Reviewer B verified all three reachable uncovered sites correct by hand, so this was coverage debt rather than a live bug; PE two-file FASTQ → uBAM is a realistic invocation.
- **Check 5's "no output dir created" clause now asserted** (`phred64_bam_input_rejected_before_output_dir_created`) — the plan-manager's sole PARTIAL. Uses a nested path `fresh_tmpdir` has not pre-created, so absence is attributable to the guard firing ahead of `ensure_output_dir`.
- **Guard→offset coupling made discoverable from both directions**: an `INVARIANT ESTABLISHED HERE` note at the guard naming every dependent site, doc bullets on the four wrapper functions pointing back, and an explicit paragraph in both `clump_only` `to_bam` doc comments recording that the `R.qual` half of their lossless invariant is enforced in `main.rs`, not locally. Reviewer B's suggested single local was not used: `main()`-scoped, it would reach only 5 of the 12 `cli.phred_offset()` sites, since the three `run_ubam_output_*` drivers take `cli` and call it themselves.
- `PHRED_OFFSET`'s doc corrected — 33 is fixed because the internal `FastqRecord` is Sanger, not because "the SAM spec fixes it" (that rule governs SAM text; BAM stores raw).

### Applied by Reviewer A during review

Two test defects, both re-verified here: `phred64_mixed_input_rejected` passed pre-fix (its `--paired` input was independently rejected by the pre-existing two-BAM check, so it could attribute nothing to the guard) — rewritten to a single-end shape that is otherwise a legal run; and the specialty test could pass vacuously on a header-only BAM. These also satisfy Reviewer B's M-2 and M-3.

### Undeclared micro-deviations, now declared

1. **Parameter named `input_phred_offset`, not `phred_offset`.** Deliberate per §4.1's naming-hazard analysis, but §13 did not say so.
2. **Check 2 (Phred+33 no-op guard) implemented at the unit layer** rather than as an integration test — it is the second arm of `bam_writer_subtracts_input_phred_offset`. No assertion lost; check 10's row already called for the same `'I'` → 40 assertion at that layer.

### Not adopted

- **A-L3** (value-domain check on the offset): `cli.phred_offset()` can only return 33 or 64, and A correctly notes a `debug_assert!` would be inert under `cargo test --release`. A `bail!` for an unreachable state adds a branch without adding safety.
- **A-L8** (`docs/.../reference/changelog.md` sync): that page is a manual copy already lagging two `Unreleased` entries; syncing it here would be inconsistent with practice and is release-time work. The H-2/H-3 corrections land in `CHANGELOG.md`, so the release-time sync will pick up the corrected text, not the false claims.
- **L-7** (hoisting the field load out of the map closure): LLVM hoists it; the repo's perf history (#248, #287) concerns per-read libm calls, not a register-resident byte.

### Follow-ups worth filing separately

- **FastQC encoding heuristic on BAM input.** For a BAM source the lowest-quality-char rule is never appropriate — BAM `QUAL` is raw Phred by spec, so the encoding is known. Any dataset whose minimum Phred ≥ 31 is reported 31 too low. Pre-existing and upstream-faithful, not introduced here, but this fix moves Phred+64 data from accidentally-correct to affected. An upstream issue against `ewels/FastQC-Rust` asking it to force Sanger for BAM input would close it.
- **Mixed-pair error message.** `--paired <fastq> <bam>` reports "two BAM files is not supported" when only one input is a BAM (surfaced by Reviewer A while rewriting F1).

### Verification after revisions

| Gate | Result |
|---|---|
`cargo fmt --all -- --check` | clean |
`cargo clippy --all-targets --release -- -D warnings` | clean |
`cargo test` | see §15 |
Recovery path (H-3) | `samtools fastq` → re-run → md5 identical to canonical |
`all-zero` claim (M-1) | 97.0% non-zero on `clock_10K_R1`, disproving it |
Multibyte invariant (H-1) | `min(93)+33` → 4 bytes for 4 bases; `saturating_add` → 5 |
