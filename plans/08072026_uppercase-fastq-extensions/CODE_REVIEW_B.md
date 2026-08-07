# Code Review B — case-insensitive FASTQ extensions (#384)

**Reviewer:** B (independent; Reviewer A reviews the same diff in parallel, no shared state)
**Target:** `fix/384-uppercase-extensions` off `dev` @ `533582a`, **uncommitted**
**Diff:** `src/io.rs` (+70/−3), `tests/integration_gzip_non_gz_extension.rs` (+111), `CHANGELOG.md` (+24)
**Plan:** `plans/08072026_uppercase-fastq-extensions/PLAN.md` (v2, §13 implementation notes)

> **Nothing was fixed.** The tree is shared with a concurrent reviewer, so every item below is a
> recommendation only. My single write is this file.

**Method.** I did not re-run the suite (it is green per §13 and proves little). Instead I probed
the built `target/release/trim_galore` (verified post-change: `SAMPLE.FASTQ.GZ` → `SAMPLE_trimmed.fq.gz`,
`1f8b`) with constructed invocations, and compiled a standalone `rustc` probe of the helper.
Scratch: `…/scratchpad/crev384b/`. No repo writes, no CWD-naming modes from the repo root.

---

## Summary

The change is correct and the helper is sound — I could not break it. No Critical or High findings.
Product behaviour matches the plan on every path I could reach, including the ones the plan
predicted would newly refuse. What I did find is **two Medium items about the change's edges rather
than its core**: one integration assertion that does not discriminate what the plan says it pins,
and one user-visible consequence class (the four specialty modes) missing from the CHANGELOG's
explicitly-enumerated list. Everything else is Low.

---

## 1. The helper — `strip_suffix_ignore_ascii_case` (`src/io.rs:521`)

**The char-boundary argument is airtight, and in fact stronger than the doc-comment claims.**

On a match, every byte of `name[idx..]` equals the corresponding suffix byte *after ASCII folding*.
Folding maps ASCII→ASCII and leaves ≥0x80 bytes untouched, so a non-ASCII name byte can only match
a non-ASCII suffix byte exactly. The suffix is itself a valid `&str`, so its first byte is never a
UTF-8 continuation byte; therefore the name byte at `idx` is never a continuation byte either, and
`idx` is a char boundary. This holds for *any* valid `&str` suffix, not just all-ASCII ones — the
all-ASCII framing in the doc-comment is a sufficient condition, not the necessary one. No change
needed; noted so a future maintainer adding a suffix does not think they are on thin ice.

`checked_sub` covers name-shorter-than-suffix; `name.as_bytes()[idx..]` cannot panic because
`idx ≤ len`; the `&str` slice is reached only on a match. Lifetimes tie the result to `name`. Clean.

**Folding is immune to non-ASCII lookalikes** — probed directly, all correctly `None`:
`x.ＦＱ` (fullwidth), `x.fı` (dotless i), `x.ß`, `x.G\u{212A}` (Kelvin sign) against `.fq`/`.gz`.
`straße.fq` → `straße` and `KELVIN\u{212A}.fq` → `KELVIN\u{212A}`: retained case and non-ASCII
bytes preserved intact.

**The naive form really does panic**, so the byte-wise choice is load-bearing rather than defensive:
`naive("😀", ".fq")` → `start byte index 1 is not a char boundary`. End-to-end confirmation that
this input is live: `😀.FASTQ.GZ` and `é.FQ.GZ` both trim to exit 0, producing `😀_trimmed.fq.gz`
and `é_trimmed.fq.gz`. The third unit test (`io.rs:1017`) therefore guards something real.

**Structure.** Two different mechanisms now express one idea: `OsStr::eq_ignore_ascii_case` in
`is_gzipped` (line 33) and the byte helper 490 lines away. That is the right call — the types differ,
and `is_gzipped` must keep working on non-UTF-8 `OsStr` — and the doc-comment cross-reference
carries the coupling. Naming mirrors `str::strip_suffix`. Placement beside its only caller is right.

## 2. Callers that must stay in step

I re-derived the caller set rather than trusting §2.3, and found no inconsistency beyond the
documented `.BGZ` asymmetry (confirmed live: `SAMPLE.FASTQ.BGZ` → `SAMPLE_trimmed.fq` plain, beside
`SAMPLE.FASTQ.BGZ_trimming_report.txt`).

- **`is_gzipped` → one decision producer.** `main.rs:345` is the sole `let gzip = …`; every namer
  (`io.rs`: SE/PE/unpaired/passthrough/clumped; `specialty.rs`: all four) takes `gzip: bool` as a
  parameter. Nothing recomputes it locally, so the fold cannot half-apply.
- **`clump_only.rs:277`/`:414`** now agree with the FASTQ-output path. Verified live (§3 below).
  `:1107` is a test helper. Note `clump_only.rs:781`/`:1019` set the same field from
  `input_is_compressed(fmt)` (content-based) on the uBAM-output path — a pre-existing
  name-vs-content split that the fold *narrows* for genuinely-gzipped uppercase files.
- **`demux.rs:118-121`** (case-sensitive `.gz`/`.fq` strips) is genuinely dead to this change: its
  argument is the trimmed output path, which this code always names lowercase — including under
  `--basename FOO`, where the stem is uppercase but the extension is still ours.
- **`specialty.rs:501-515`** — the `_R1`/`R2` strip. Newly *reachable* for uppercase input: before,
  the stem ended `.FASTQ` so the strip always missed. See L4/M2.
- **`main.rs:1645`, `:2219`** (interleaved-uBAM stems) use raw `file_stem()`. See L5.

## 3. The new tests

All three assert real, verified behaviour. I reproduced each with a constructed invocation:

| Test | Reproduced | Result |
|---|---|---|
| V5 `uppercase_extension_names_and_compresses_like_lowercase` | `-o out SAMPLE.FASTQ.GZ` | `SAMPLE_trimmed.fq.gz` (`1f8b`) + `SAMPLE.FASTQ.GZ_trimming_report.txt`, exit 0 |
| V4 `uppercase_and_lowercase_spellings_of_one_stem_now_collide` | both iterations | exit 1, `Output path collision …` **on stderr**, `out/` empty |
| V8 `clump_only_uppercase_gz_reports_gzip_and_ratio` | `--clump_only -o out SAMPLE.FASTQ.GZ` | `SAMPLE_clumped.fq.gz`; report `Input: … (gzip, 273332 bytes)` + `Compression ratio: 0.77x` |

**Discrimination spot-check (V2/`is_gzipped`), by inspection:** under `ext == "gz"`,
`is_gzipped("SAMPLE.FASTQ.GZ")` is `false`, so `io.rs:983` fails on its first line. Airtight — no
build needed. Same for V5: an unfolded stem names the file `SAMPLE.FASTQ_trimmed.fq*`, so
`trimmed.is_file()` fails at path resolution, exactly as §13 records.

**V4 on this machine's APFS:** the three `write_gz` calls create only **two** files —
`SAMPLE.FQ.GZ` truncates the existing `sample.fq.gz` in place, leaving one entry named
`sample.fq.gz` whose content is `"C"` (verified). The test still passes and still discriminates
on both platforms, because `preflight_output_collisions` is lexical and runs before any input is
opened; iteration 2's `SAMPLE.FQ.GZ` argument need not exist as a directory entry. See L3.

## 4. Interaction with #383/#385's collision keys

No surprise beyond the plan's §3.1 table. `path_identity_key` stays case-sensitive, so
`SAMPLE.FASTQ.GZ` + `SAMPLE.FQ.GZ` are two inputs, not a duplicate; the refusal comes from
`collision_key` on the outputs, which is the intended route and the message the test asserts.
Nothing in the diff can create a new *self*-overwrite (output ⊂ input) class: every namer appends
`_trimmed`/`_clumped`/`_val_*` to the stem, so a planned output can never fold onto its own input.

One class the §3.1 table does not cover, because it is out of a per-run pre-flight's reach: see L6.

---

## Findings by priority

### Critical — none

### High — none

### Medium

**M1 — V8's positive assertion cannot fail for the reason the plan gives it.**
`tests/integration_gzip_non_gz_extension.rs:454`, `assert!(report.contains("gzip"))`. The real
report (captured above) contains `Output: … (gzip level 1, …)`, so that substring is present
regardless of the `Input:` label — the assertion is satisfied by the output half. The complementary
`!report.contains("plain")` (`:457`) is likewise driven by the *output* label: under the §13
negative control (`is_gzipped` reverted) `gzip_output` is false, so the Output line reads `plain`
and the assertion fails for the output reason, not the input one. Net effect: **V8 pins the stem
and the output compression, but nothing in it pins C1's `Input: plain → gzip` shift**, which the
plan (§9 V8, §12 I6) names as its purpose. Recommend asserting the line, not the file:

```rust
let input_line = report.lines().find(|l| l.starts_with("Input:")).unwrap();
assert!(input_line.contains("(gzip,"), "{input_line}");
```

**M2 — the CHANGELOG's "Four visible consequences" omits a fifth: the specialty modes.**
`--hardtrim5/3`, `--clock` and `--implicon` all derive their stems from `strip_fastq_extensions`
(`specialty.rs:449`, `:471`, `:485`, `:501`) and their compression from the same `gzip` flag, so
uppercase input is renamed *and* newly compressed there too. Verified with the real binary
(run from scratch, `-o` given):

- `--hardtrim5 30 SAMPLE_R1.FASTQ.GZ` → `SAMPLE_R1.30bp_5prime.fq.gz`
  (was `SAMPLE_R1.FASTQ.30bp_5prime.fq`, plain).
- `--implicon SAMPLE_R1.FASTQ.GZ SAMPLE_R2.FASTQ.GZ` → `SAMPLE_8bp_UMI_R1.fastq.gz` /
  `…_R2.fastq.gz` (was `SAMPLE_R1.FASTQ_8bp_UMI_R1.fastq`) — note this name changes **twice**:
  the extension folds *and* the `_R1` strip at `specialty.rs:505` becomes newly reachable, which
  is a second, unmentioned improvement.

These four modes default to naming output into the **CWD** (the very next CHANGELOG bullet is about
that), so a renamed output there is the least recoverable of the set. Recommend either generalising
bullet 1 or adding a fifth. Also the shallower half of the same gap: no committed test covers a
specialty namer under uppercase input — the shared helper is unit-tested, which is defensible, but
the `_R1`-strip interaction is not covered anywhere.

### Low

**L1 — names that are *only* an extension yield an empty stem.** `.FQ.GZ` → `_trimmed.fq.gz`
(verified, exit 0; report `.FQ.GZ_trimming_report.txt`). Pre-existing for lowercase `.fq.gz`; the
fold makes a second spelling class reach it. Multi-input is caught, if oddly worded:
`Output path collision …: out/_trimmed.fq.gz and out/_trimmed.fq.gz`. No committed test pins it.

**L2 — stale comment adjacent to the new code.** `io.rs:975`: "Heuristic is extension-based
(matches FastqReader's gzip detection)". Since #374 the reader sniffs content, and `is_gzipped`'s
own doc-comment (lines 13-17) says so explicitly — the two now contradict each other, and the new
`test_is_gzipped_folds_case` was added directly beneath the stale line. Drive-by one-liner.

**L3 — V4's fixture set is filesystem-dependent, undocumented.** Only two of its three files exist
on a case-insensitive volume (§3). Harmless today; recommend one line in the test comment, because
adding any *content*-based assertion there later would pass on Linux and fail on macOS.
(Separately, `fresh_tmpdir` uses a fixed path per slug — pre-existing in this file, but three more
slugs now share the hazard of two concurrent `cargo test` runs.)

**L4 — `--implicon`'s `_R1` strip stays case-sensitive** (`specialty.rs:505-515`), so
`sample_r1.FQ.GZ` folds its extension but keeps the lowercase `_r1` marker. Pre-existing
asymmetry, but the fold is what makes the strip reachable at all for uppercase-extension input, so
this is the moment it becomes visible. Out of scope for #384; worth an issue.

**L5 — the plan's "cannot diverge" for `main.rs`'s `file_stem()` sites is slightly overstated.**
`main.rs:1645` / `:2219` (interleaved-uBAM → FASTQ / uBAM out) use raw `file_stem()`, so a
content-detected uBAM named `data.FASTQ.GZ` yields `data.FASTQ_val_1.fq.gz` there while the SE
dispatch yields `data_trimmed.*`. The divergence is **pre-existing** (lowercase `data.fastq.gz`
already split `data` from `data.fastq`); the fold only extends it to uppercase names. Contrived
input, no action needed — recorded because §2.4 asserts the opposite.

**L6 — cross-invocation aliasing on case-insensitive filesystems.** `-o out d1/sample.fq.gz` then
`-o out d2/SAMPLE.FQ.GZ` now target one file (`out/sample_trimmed.fq.gz` ≡ `out/SAMPLE_trimmed.fq.gz`
on APFS) where before the second wrote `SAMPLE.FQ_trimmed.fq`. The per-run pre-flight cannot see
across invocations. The class is pre-existing (`d1/sample.fq.gz` + `d2/SAMPLE.fq.gz` already
aliased); the fold adds mixed-*spelling* members. Inherent to the goal, not a defect.

**L7 — compression is run-wide, and the CHANGELOG bullet reads per-file.**
`gzip` comes from `input[0]` only, so `trim_galore A.FASTQ.GZ b.fastq` now gzips `b`'s output too
(it was plain before the fold). Pre-existing "first input decides" rule, new trigger. One clause in
bullet 2 would cover it.

**L8 — `.then(…)` must stay lazy.** `then_some(&name[..idx])` would evaluate the slice on the
non-matching path and reintroduce the `😀` panic. Clippy's `unnecessary_lazy_evaluations` does not
fire on indexing, so the risk is low and the doc-comment covers the reasoning; no action.

**L9 — `test_is_gzipped_folds_case` drops one plan V2 case:** lowercase `.bgz` → false is asserted
nowhere (only `.BGZ`). Substantively covered, since `.BGZ` false implies `.bgz` false under
folding. Cosmetic.

---

## Verification of the plan's own claims

- **Validation matrix unaffected** — re-confirmed by grep as §7 asked: no uppercase extension in
  `.github/workflows/ci.yml`, and no fixture in `test_files/` carries one (`PolyA.fastq.gz` and
  friends have uppercase *stems*, which the fold preserves — pinned by the `Sample.FastQ.Gz` case).
- **Report filename unchanged** — confirmed on every probe; the CHANGELOG's closing sentence and
  its refusal to claim MultiQC grouping are both accurate.
- **§13's six §5 steps** are all present in the diff as described; I found no undocumented deviation.
