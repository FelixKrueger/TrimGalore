# Plan — Match FASTQ extensions case-insensitively (#384)

**Issue:** [#384](https://github.com/FelixKrueger/TrimGalore/issues/384)
**Branch point:** `dev` @ `533582a`
**Revision:** v2 (2026-08-07) — see §12

---

## 1. Goal

Make an uppercase FASTQ extension behave exactly like its lowercase equivalent:
`SAMPLE.FASTQ.GZ` should produce `SAMPLE_trimmed.fq.gz` and
`SAMPLE.FASTQ.GZ_trimming_report.txt`, mirroring
`sample.fastq.gz` → `sample_trimmed.fq.gz`. One rule, no case-dependent behaviour to
document or remember.

Decision already recorded on the issue: fold case, rather than pin the current behaviour.
The reason both #382 reviewers gave for pinning — byte-identity with Perl v0.6.11 — does not
apply, because the current output is not Perl's either (§2.2).

---

## 2. Context

### 2.1 What is wrong today

Verified on `dev` @ `533582a`:

```
$ trim_galore -o out --dont_gzip SAMPLE.FASTQ.GZ
out/SAMPLE.FASTQ_trimmed.fq
out/SAMPLE.FASTQ.GZ_trimming_report.txt      ← disagree about the sample name
```

`strip_fastq_extensions` (`io.rs`) matches `[".gz", ".bgz", ".bgzf"]` and `[".fastq", ".fq"]`
with `strip_suffix`, which is case-sensitive. `.FASTQ.GZ` matches none, so the function falls
through to `Path::file_stem()`, which removes only the final component and leaves `SAMPLE.FASTQ`.

The file is *read* correctly — since #374 gzip-ness comes from the first three bytes, not the
name — which is why the run succeeds and the inconsistency is visible rather than fatal.

### 2.2 Perl v0.6.11 did not produce this mismatch

Its output naming (`trim_galore` lines 909-923) tests four suffixes case-sensitively and then
falls through to `else { s/$/_trimmed.fq/ }`, appending to the **full** filename; the report
(line 611) appends the same way. So Perl wrote `SAMPLE.FASTQ.GZ_trimmed.fq` beside
`SAMPLE.FASTQ.GZ_trimming_report.txt` — ugly but consistent.

Perl was case-sensitive, and its case-sensitivity still produced agreeing names, because its
fallback appended where ours strips. **The naming mismatch is ours**, and on naming, parity is
already lost.

Compression is the opposite case, and the plan must say so plainly: Perl's decision (line 925,
`if ($gzip or $filename =~ /\.gz$/)`, no `/i`) is case-sensitive, so Perl wrote `.GZ`-named
output **plain** — exactly what `dev` does today. Folding `is_gzipped` therefore *breaks* a
parity that currently holds, while folding the stem restores nothing Perl had. The fold is
chosen on internal consistency (§2.3), with the parity cost stated rather than leaned past.
The CI byte-identity matrix is unaffected only because no fixture is uppercase.

Restoring Perl's append-fallback would be both parity-correct and consistent, but is
unavailable: the fallback now also serves uBAM input, where `file_stem` is right
(`sample.bam` → `sample_trimmed.bam`, not `sample.bam_trimmed.bam`), and uBAM is a v2 feature
Perl never had. Separating the two populations requires deciding the case question anyway.

### 2.3 The adjacent function that makes this more than a one-liner

`is_gzipped` (`io.rs:29-31`) is also case-sensitive:

```rust
path.extension().is_some_and(|ext| ext == "gz")
```

It decides **output compression** — `gzip = !cli.dont_gzip && naming::is_gzipped(&cli.input[0])`
(`main.rs:345`) — and it has three more call sites, all in `clump_only.rs`: `:277` (SE) and
`:414` (PE) set `input_compressed`, which drives the `Input: … (gzip|plain)` label and gates the
`Compression ratio:` line in the `--clump_only` report; `:1107` is a test helper that
gzip-encodes by name. So folding also changes `--clump_only` **report content** for `.GZ`-named
input — `plain` → `gzip`, plus a ratio line. Benign in direction (the input genuinely is gzip)
but user-visible, so it goes in the CHANGELOG. Its doc-comment is explicit that it is
filename-based on purpose: it answers
"should the output be compressed", mirroring the input's name, not "is the input compressed",
which `format.rs` decides from content.

So folding only `strip_fastq_extensions` would trade the naming mismatch for a compression
mismatch:

| input | stem | output today | stem-fold only | both folded |
|---|---|---|---|---|
| `sample.fastq.gz` | `sample` | `sample_trimmed.fq.gz` | unchanged | unchanged |
| `SAMPLE.FASTQ.GZ` | `SAMPLE.FASTQ` | `SAMPLE.FASTQ_trimmed.fq` (plain) | `SAMPLE_trimmed.fq` (**still plain**) | `SAMPLE_trimmed.fq.gz` |

A silently uncompressed output is a worse surprise than a mismatched stem, particularly at
whole-run scale. Both functions fold, or neither. **This plan folds both.**

### 2.4 Extension handling elsewhere, deliberately untouched

- `demux.rs:118-121` strips `.gz` then `.fq` case-sensitively, but its input is the *trimmed
  output path*, which this code generates and always names with lowercase `.fq` / `.fq.gz`.
  Folding there would be dead code.
- `format.rs` decides input format from content (gzip magic, `BAM\1`), not names. Unaffected.
- `clump_only.rs` is **affected**, not cleared — v1 cited a test assertion (`:1469`) and missed
  the real consumers at `:277`/`:414` (§2.3). `clumped_output_name` also derives its stem via
  `strip_fastq_extensions` while `clumping_report_name` uses the full filename, so `--clump_only`
  carries the same stem/report mismatch #384 reports, and this change fixes it there too.
- Reviewer sweep (`grep -rn 'extension()\|file_stem()\|strip_suffix\|ends_with("' src/*.rs`)
  found no further filename-driven decision that could disagree with the fold: `main.rs:1646` /
  `:2219` and `--basename` read names but cannot diverge.

---

## 3. Behavior

1. `strip_fastq_extensions` matches the gzip-family suffixes and the FASTQ extensions
   case-insensitively, preserving the **rest** of the name's case: `SAMPLE.FASTQ.GZ` → `SAMPLE`,
   `MySample.FastQ.Gz` → `MySample`.
2. `is_gzipped` returns true for any ASCII-case variant of `gz`.
3. Everything not matching a FASTQ extension keeps `Path::file_stem()`, so `.bam` and unrelated
   extensions are untouched.
4. Output extensions this code *writes* stay lowercase (`_trimmed.fq.gz`), as now. Only
   *matching* becomes case-insensitive; nothing starts emitting uppercase.

### 3.1 Edge cases

| Input | Stem | Output | Note |
|---|---|---|---|
| `sample.fastq.gz` | `sample` | `sample_trimmed.fq.gz` | unchanged — the common path |
| `SAMPLE.FASTQ.GZ` | `SAMPLE` | `SAMPLE_trimmed.fq.gz` | **changed**: was `SAMPLE.FASTQ_trimmed.fq` |
| `Sample.FastQ.Gz` | `Sample` | `Sample_trimmed.fq.gz` | mixed case folds too |
| `SAMPLE.FQ` | `SAMPLE` | `SAMPLE_trimmed.fq` | plain in, plain out |
| `SAMPLE.FASTQ.BGZ` | `SAMPLE` | `SAMPLE_trimmed.fq` | `.BGZ` is not `.gz`, so still plain — consistent with `.bgz` today (#374/#381) |
| `sample.BAM` | `sample` | per uBAM naming | fallback keeps `file_stem`; unchanged |
| `sample.txt` | `sample` | `sample_trimmed.fq` | fallback; unchanged |
| `SAMPLE.FASTQ.GZ` + `sample.fastq.gz` in one run | — | rejected **already** | report names fold to equal keys today; the fold changes nothing (and on APFS the two names are one file) |
| `SAMPLE.FASTQ.GZ` + `sample.fq.gz` in one run | `SAMPLE` / `sample` | **newly rejected** | accepted today (report names differ); after the fold both trimmed outputs collide, on every filesystem |
| `SAMPLE.FASTQ.GZ` + `SAMPLE.FQ.GZ` in one run | `SAMPLE` / `SAMPLE` | **newly rejected** | same mechanism, not even case-dependent — a true duplicate output |

v1's claim that pure case-variants become "reachable where they were not before" was false —
they are already refused via the folded *report* keys, so a V4 built on them could never fail.
The genuinely new refusals are the mixed-spelling pairs above, which are integration-testable on
any filesystem; they are user-visible and belong in the CHANGELOG (a run that succeeded now
errors).

---

## 4. Signatures

No signature changes. Two bodies change:

```rust
pub fn is_gzipped(path: &Path) -> bool
pub fn strip_fastq_extensions(path: &Path) -> String
```

A small private helper keeps the matching in one place:

```rust
/// ASCII-case-insensitive suffix strip, preserving the case of what remains.
fn strip_suffix_ignore_ascii_case<'a>(name: &'a str, suffix: &str) -> Option<&'a str>
```

Implemented on **bytes**, slicing only after a successful match:

```rust
fn strip_suffix_ignore_ascii_case<'a>(name: &'a str, suffix: &str) -> Option<&'a str> {
    let idx = name.len().checked_sub(suffix.len())?;
    name.as_bytes()[idx..]
        .eq_ignore_ascii_case(suffix.as_bytes())
        .then(|| &name[..idx])
}
```

The naive `&str` form (`name[idx..].eq_ignore_ascii_case(...)`) panics when `idx` lands inside a
multi-byte character — verified with a probe: a file named `😀` works end-to-end today
(`😀_trimmed.fq`, exit 0) and would panic under the naive helper. Matching all-ASCII suffix bytes
guarantees `idx` is a char boundary, and `checked_sub` covers name-shorter-than-suffix. The
retained part keeps its original case. ASCII-only, matching `norm_path`'s convention (A1).

---

## 5. Implementation outline

**Step 1 — `src/io.rs`: add `strip_suffix_ignore_ascii_case`** beside `strip_fastq_extensions`.

**Step 2 — `strip_fastq_extensions`** uses it for both lists, leaving the `file_stem` fallback
and the two-sequential-strips structure (#382, @BenjaminDEMAILLE) intact.

**Step 3 — `is_gzipped`**: `ext.eq_ignore_ascii_case("gz")`. Extend its doc-comment to say
matching is case-insensitive and why it must agree with `strip_fastq_extensions` (§2.3).

**Step 4 — unit tests** in `io.rs`: §5 V1–V3 below.

**Step 5 — integration test** in `tests/integration_gzip_non_gz_extension.rs` (where #374/#382
put the related cases): one end-to-end uppercase run asserting the trimmed file and the report
share a stem, and that the output is gzipped.

**Step 6 — CHANGELOG** under `#### Changes`: the rename, the compression change, and the note
that this does not restore v0.6.11 parity either.

**Comment discipline:** one line per new comment; reasoning goes in the commit message.

---

## 6. Efficiency

`eq_ignore_ascii_case` on a suffix is O(len(suffix)) with no allocation, against
`strip_suffix`'s O(len(suffix)) memcmp. Five suffix probes per filename, once per input file.
Unmeasurable.

---

## 7. Integration

**Validation matrix — unaffected.** No fixture in `test_files/` and no `trim_galore` invocation
in `ci.yml` uses an uppercase extension, so no Perl byte-identity comparison changes. Worth
re-confirming by grep during implementation rather than trusting this sentence.

**Downstream:** a pipeline passing uppercase input gets differently-named *and* newly-compressed
output; some previously-successful multi-input runs become refusals (§3.1); and `--clump_only`
report content shifts for `.GZ` input (§2.3). All corrections, all visible: four bullets under
`#### Changes`, not two.

**MultiQC, stated honestly:** the fold makes uppercase output *internally* consistent; it does
not deliver downstream grouping. The report filename still carries `.FASTQ.GZ` (by design — full
input filename), so a case-sensitive extension-cleaner still derives a different sample name from
the report than from the trimmed file. Do not claim grouping is fixed in the CHANGELOG.

**#383's pre-flight:** gains reachable collisions it could not previously see (§3.1 last row).
That is the pre-flight working as designed.

---

## 8. Assumptions

- **A1 (fixed).** ASCII case folding only, consistent with `norm_path` and the rest of the crate.
  A Turkish-locale `İ` or similar is out of scope; FASTQ extensions are ASCII.
- **A2 (fixed).** Only *matching* becomes case-insensitive. Names this code writes stay
  lowercase, so no output filename gains uppercase.
- **A3 (fixed, weakened in v2).** Folding the stem without folding `is_gzipped` produces
  silently uncompressed output and is worse than folding neither — that direction is firm. The
  converse is softer: folding `is_gzipped` alone would be coherent but does not fix #384. And
  "both or neither" is not the whole option space — sourcing the compression decision from
  detected *content* (the follow-up #374/#381 deferred) would retire this coupling entirely;
  recorded in §10 as the direction, out of scope here.
- **A4 (fixed).** `.BGZ`/`.BGZF` fold as part of the gzip-family list, but — like lowercase
  `.bgz` — do not imply gzipped *output*, because `is_gzipped` tests only `gz`. That asymmetry
  is pre-existing and documented in `is_gzipped`'s doc-comment; this change does not widen it.
- **A5 (decision).** The renamed outputs are accepted as a behaviour change. The affected users
  are those currently receiving a mismatched pair.

---

## 9. Validation

**V1 — unit, `strip_fastq_extensions`.** `SAMPLE.FASTQ.GZ`, `Sample.FastQ.Gz`, `SAMPLE.FQ`,
`SAMPLE.FASTQ`, `SAMPLE.FQ.GZ`, `SAMPLE.FASTQ.BGZ`, `sample_R1.FQ.GZ` → stems `SAMPLE`,
`Sample`, `SAMPLE`, `SAMPLE`, `SAMPLE`, `SAMPLE`, `sample_R1`. The mixed-case case is the one
that proves the *retained* part keeps its case rather than being lowercased.

**V2 — unit, `is_gzipped`.** True for `.gz`, `.GZ`, `.Gz`; false for `.fastq`, `.bgz`, `.BGZ`,
and `sample.gz.fastq`.

**V3 — unit, regression guards.** Every existing `test_strip_fastq_extensions` and
`test_strip_fastq_extensions_bgz` case still passes unchanged; `sample.bam` / `sample.txt` keep
their current stems (the fallback must not start folding); and a non-ASCII name (`😀`, and one
shorter than the suffix) returns unchanged rather than panicking — the naive slice form panics
on exactly this input (§4), so this case is the guard against reintroducing it.

**V4 — integration, the new refusals.** `SAMPLE.FASTQ.GZ` + `sample.fq.gz` (accepted today,
verified) exits non-zero after the fold with the collision message; same for
`SAMPLE.FASTQ.GZ` + `SAMPLE.FQ.GZ`. Both run on any filesystem. v1's pair (pure case-variants)
is already rejected today via folded report keys and cannot discriminate — do not use it.

**V5 — integration, end-to-end.** Run **without** `--dont_gzip` (v1's §2.1 reproduction used it,
which hides the compression half). Assert exact literal filenames with no call to any function
under test in computing the expectations: `SAMPLE_trimmed.fq.gz` exists,
`SAMPLE.FASTQ.GZ_trimming_report.txt` exists, and the trimmed file's first two bytes are
`1F 8B`. Resolve which candidate path exists before the magic-byte check, so the naming half and
the compression half fail independently and legibly. ("Share a stem" is not assertable literally:
the report keeps the full input filename by design — see §7 note on MultiQC.)

**V6 — negative controls (manual one-offs, not committed guards).** Revert the stem fold alone →
V1 fails on its first case. Revert `is_gzipped` alone → V5 fails (via its path-resolution step —
the split assertion in V5 exists so the failure is legible). The *committed* guard for the
`is_gzipped` half is V2's `is_gzipped(".GZ") == true`; say so, so V6 is not mistaken for a
regression test.

**V8 — integration, `--clump_only` (the least-guarded consumer).** A `.GZ`-named input: the
report's `Input:` line says `gzip`, a `Compression ratio:` line is present, and the output is
`<stem>_clumped.fq.gz`. Pins C1's report-content change and the `--clump_only` stem fix.

**V7 — gates.** `cargo fmt --all -- --check`; `cargo clippy --all-targets --release -- -D
warnings`; `cargo test`; then rebuild and re-run the §2.1 reproduction.

---

## 10. Questions and ambiguities

**Open 1 — resolved.** Both reviewers independently answered **fold both**, with the same two
qualifications: (a) B found the compression half currently *has* Perl parity (§2.2), so the fold
knowingly breaks a parity the naming half never had — stated, not hidden; (b) both name
content-sourced compression as the follow-up that retires the coupling (out of scope, recorded
here as the direction).

**Open 2 — `#### Changes` or a minor version note?** Two visible changes for uppercase users.
The CHANGELOG entry is written either way; flagging in case a release note is wanted.

**Out of scope:** `demux.rs`'s case-sensitive strips (§2.4, would be dead code); restoring
Perl's append-fallback (§2.2, collides with uBAM naming); non-ASCII case folding (A1).

---

## 11. Self-Review

**Logic.** Traced both functions to their callers: `strip_fastq_extensions` feeds every output
namer in `io.rs` and all four in `specialty.rs`; `is_gzipped` feeds only `main.rs:295`'s `gzip`
flag, which is threaded to writers and namers alike. Folding both keeps them consistent, which
§2.3 shows is the actual requirement.

**Adjusted during review.** The first draft folded `strip_fastq_extensions` only — which was the
scope stated on the issue and in my own comment. Tracing `is_gzipped` showed that would leave
uppercase input producing silently uncompressed output, a worse outcome than the mismatch being
fixed. §2.3, A3 and V6 exist because of that.

**Edge cases.** §3.1, including the new collision the fold makes reachable, and the `.BGZ`
row where folding the stem does *not* imply folding compression.

**Remaining risks.** Open 1 is the real one. Beyond that, the change is small and its blast
radius is bounded by two functions whose callers are enumerated above.

---

## 12. Revision history

**v1 → v2**, after dual independent plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`); the
two load-bearing findings were re-verified against the binary and the v0.6.11 source before
adoption.

Adopted from **A**: C1 — v1's caller enumeration for `is_gzipped` missed all three
`clump_only.rs` sites and mis-cited the `main.rs` line, and §2.4 had "cleared" `clump_only.rs`
by checking a test assertion instead of the consumers (the exact §2.4-incompleteness failure the
reviewers were asked to hunt); C2 — v1's "newly reachable collision" was already rejected today
via folded report keys, so V4 could never fail (verified: exit 1 on dev, same inode on APFS);
C3 — the naive helper panics on a multi-byte boundary for filenames that work today (byte-wise
form adopted, non-ASCII case added to V3).

Adopted from **B**: I2 — Perl's *compression* decision is case-sensitive too, so `dev` is
parity-correct on compression and the fold breaks that parity knowingly (§2.2 rewritten; the
"parity is already lost" argument held for naming only); I1 — the genuinely new refusals are
mixed-spelling pairs, integration-testable everywhere (V4 rewritten); I3/I5 — V5 asserts literal
filenames without `--dont_gzip`, split so the two halves fail independently; I4 — MultiQC
grouping is *not* delivered and the CHANGELOG must not claim it; I6/V8 — `--clump_only` report
content pinned; V6 relabelled a manual one-off with V2 as the committed guard.

Adopted from **both**: A3 weakened from a biconditional to the one direction the argument
supports, with content-sourced compression recorded as the retiring follow-up; the CHANGELOG
grows to four bullets (rename, compression flip, new refusals, `--clump_only` shifts).

**Open 1 was put to both reviewers independently: both said fold both.**

---

## 13. Implementation notes

Implemented on branch `fix/384-uppercase-extensions` off `dev` @ `533582a`. All six §5 steps
done as specified in v2 — no deviations. Gates: `cargo fmt --check` clean,
`clippy -D warnings` clean, **538 tests pass** (533 prior + 3 unit + 3 integration, minus one
superseded count drift — 6 new tests total).

Diff: `src/io.rs` (+70/−3: the byte-wise helper, both folds, doc-comment, 3 unit tests),
`tests/integration_gzip_non_gz_extension.rs` (+111: V5 literal-filename/magic-byte,
V4 mixed-spelling refusals ×2, V8 `--clump_only` report), `CHANGELOG.md` (+24: four-bullet
entry under `#### Changes`, with the Perl compression-parity departure stated per B's I2 and
no MultiQC grouping claim per B's I4).

### Negative controls (V6), both run and both discriminating

| Control | Expected failures | Result |
|---|---|---|
| stem fold reverted alone | V1 unit + V5 integration | both FAILED; reverted → green |
| `is_gzipped` reverted alone | V2 unit + V5 + V8 | all three FAILED; reverted → green |

V5 failed at path resolution under control 2, as B's I6 predicted — legible because the naming
and compression assertions are split.

§2.1's reproduction re-run **without** `--dont_gzip` (B's I5): `SAMPLE_trimmed.fq.gz` with
`1f 8b` magic beside `SAMPLE.FASTQ.GZ_trimming_report.txt`, exit 0.

### Post-code-review round

Dual code review + coverage audit on the implementation. **0 Critical, 0 High** from both
reviewers; coverage INCOMPLETE by exactly one item. All findings applied:

- **Coverage / A-L6 (found by three passes independently):** V2's lowercase-`.bgz` false case
  was asserted nowhere — added, with the comment stating what it guards (the fold must not
  widen `is_gzipped` to the stem's gzip-family list).
- **B-M1 (verified against the real report):** V8's `contains("gzip")` was satisfied by the
  *Output* line (`gzip level 1`) regardless of the Input label. Now asserts the `Input:` line
  itself contains `(gzip,`.
- **A-M1 ≡ B-M2:** the CHANGELOG omitted the specialty modes — fifth bullet added, naming the
  newly-reachable `--implicon` `_R1` strip (same mechanism as #381's `.bgz` note).
- **A-M2:** paired uppercase integration test added, restoring the test file's stated
  dispatch-path matrix and pinning the `input[0]`-drives-run-wide-`gzip` coupling.
- **B-M2 (second half):** `uppercase_stem_reaches_specialty_output_names` unit test pins the
  implicon/hardtrim/clock namers, beside the #382 `.bgz` pin.
- **A-L1:** comment on the V4 fixture recording that two of its three names alias one file on
  APFS (both iterations still refuse on both filesystem kinds; content assertions forbidden).
- **A-L5:** V5's magic-byte check uses `starts_with`, so a short file fails the assert rather
  than panicking on the slice.
- **B-L2 / A-L3:** the stale "matches FastqReader's gzip detection" comment corrected.

One iteration during the batch: the specialty-test insertion anchored on the `fn` line and
landed between the existing test's `#[test]` attribute and its function — duplicated attribute
plus dead code, caught by clippy at the gate. The first repair then re-parented the #381
doc-comment onto the new test (the same defect #385's D6 introduced in `demux.rs`); restructured
so each test keeps its own doc-comment.

Final gates: **540 tests**, fmt clean, `clippy -D warnings` clean. Diff: 4 files, +269/−5.
