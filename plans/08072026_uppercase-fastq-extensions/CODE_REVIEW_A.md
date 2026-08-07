# Code Review A — #384 case-insensitive FASTQ extensions

**Target:** `fix/384-uppercase-extensions` off `dev` @ `533582a`, uncommitted.
`src/io.rs` (+70/−3), `tests/integration_gzip_non_gz_extension.rs` (+111), `CHANGELOG.md` (+24).
**Plan:** `plans/08072026_uppercase-fastq-extensions/PLAN.md` v2 (§13 implementation notes).
**Reviewer:** A, fresh context, independent of Reviewer B.

**Nothing was fixed.** The tree is shared with a concurrent reviewer; this file is my only write.
Every finding below is a recommendation.

**Method.** I did not re-run the suite (already green, proves little). Instead I drove the built
binary through 15 constructed invocations from
`/private/tmp/…/scratchpad/crev384a/` and hand-verified every assertion the three new
integration tests make, plus the cases the plan reasons about but does not test. Findings marked
*verified* were observed, not inferred.

**Verdict: no Critical, no High.** The change is correct and does what §3.1 says. Two Medium
(one CHANGELOG completeness, one test-matrix convention) and seven Low.

---

## What holds up

**The helper's boundary argument is airtight — and more general than its comment claims.**

```rust
fn strip_suffix_ignore_ascii_case<'a>(name: &'a str, suffix: &str) -> Option<&'a str> {
    let idx = name.len().checked_sub(suffix.len())?;
    name.as_bytes()[idx..].eq_ignore_ascii_case(suffix.as_bytes()).then(|| &name[..idx])
}
```

The comment argues from "a matched all-ASCII tail". The stronger statement also holds: for **any**
valid `&str` suffix, a byte-wise match forces `idx` to be a char boundary. `suffix`'s first byte is
either ASCII or a UTF-8 *leading* byte; a non-boundary position in a valid `&str` is a
*continuation* byte (0x80–0xBF), and `eq_ignore_ascii_case` folds only A–Z/a–z, i.e. ASCII→ASCII —
it can never make a continuation byte compare equal to a leading byte. So `&name[..idx]` cannot
panic, for any caller, present or future. `checked_sub` covers the short-name case. The slice is
only reached after a match, so the failing probes never index at all.

- **Non-ASCII lookalikes are immune** because the fold is byte-wise: no multi-byte sequence can
  fold to an ASCII letter, so `ß` (C3 9F) and the Kelvin sign (E2 84 AA) cannot match `z`/`K`.
  *Verified end-to-end:* `ß.FQ.GZ` → `ß_trimmed.fq.gz` (retained bytes preserved, gzipped);
  `😀.FQ.GZ` and a bare `😀` both exit 0.
- **The naive `&str` form really would panic**, so `test_…_non_ascii_and_short_names` guards
  something real: `"😀"` with `".gz"` gives `idx = 1`, a continuation byte. And `"é.fq"` exercises a
  *successful* strip whose retained part is multi-byte — the better of the two cases.
- **All four `is_gzipped` call sites fold together** (`main.rs:345`, `clump_only.rs:277`, `:414`,
  `:1107`). I re-ran the plan's §2.4 sweep over `src/*.rs`; the only other filename-driven
  decisions are `main.rs:1646`/`:2219` (interleaved-uBAM `file_stem()`, agrees with the fallback)
  and `specialty.rs:505-512` (`_R1`/`R2` mate token, case-sensitive — pre-existing, orthogonal,
  and the fold *improves* it: `S_R1.FASTQ.GZ` now reaches the `_R1` strip at all).
- **`demux.rs`'s dead-code claim (§2.4) checks out.** `demux_base_name` consumes the *generated*
  trimmed path, which this code always names lowercase `.fq`/`.fq.gz` — including under
  `--basename`, where the user string is a prefix, never the extension.
- **The Perl byte-identity matrix is untouched.** No fixture in `test_files/` and no
  `trim_galore` invocation in `ci.yml` carries an uppercase extension; no existing assertion pins
  the old uppercase behaviour.
- **All three new tests' assertions reproduce by hand.** V5: `out/SAMPLE_trimmed.fq.gz` with
  `1f 8b` beside `SAMPLE.FASTQ.GZ_trimming_report.txt`. V4: both pairs exit non-zero with
  `Output path collision …`, `out/` empty. V8: report reads `Input:  … (gzip, 110 bytes)`,
  `Compression ratio: 1.22x`, output `SAMPLE_clumped.fq.gz`.
- **The documented asymmetry stays put.** `SAMPLE.FASTQ.BGZ` → `SAMPLE_trimmed.fq`, plain
  (*verified*) — same as lowercase `.bgz`, exactly as A4 says.
- **Positive control for V4:** `SAMPLE.FASTQ.GZ` + `B.fq.gz` still succeeds, so the refusal is
  caused by stem equality, not by uppercase input per se.

---

## Medium

### M1 — the CHANGELOG enumerates "Four visible consequences" but four is not all of them

The entry lists SE naming, compression, new refusals, and `--clump_only`. It omits the four
specialty namers and `--passthrough`, which all derive their stem from `strip_fastq_extensions`
(`specialty.rs:449/471/485/501`, `io.rs:329`) and therefore also rename for uppercase input.
*Verified:*

| invocation on `P_R1.FASTQ.GZ` | before | after |
|---|---|---|
| `--hardtrim5 10 --paired` | `P_R1.FASTQ.10bp_5prime.fq` | `P_R1.10bp_5prime.fq.gz` |
| `--clock` | `P_R1.FASTQ.clock_UMI.R1.fq` | `P_R1.clock_UMI.R1.fq.gz` |

These three modes name output into the **current working directory** (the CHANGELOG says so 20
lines further down), which is where a positional `ls`/glob is most likely to be depended on — so
this is the sharper half of the rename, not the softer one. `--implicon` and `--passthrough` are
the same shape.

**Recommend:** a fifth bullet, or replace "Four visible consequences" with "Visible consequences"
and add one line covering the specialty and passthrough namers. Cheap, and it removes an
exhaustiveness claim that is currently false.

### M2 — the new tests skip the dispatch paths this file exists to cover

`tests/integration_gzip_non_gz_extension.rs`'s own module docstring states the selection rule and
why: tests are chosen *by dispatch path*, because in #374 "nine green checks … missed a
`--paired --cores 1` regression". Its table is SE / paired `--cores 1` / paired `--cores 2` /
`--clump_only`. #384 adds SE and `--clump_only` only — no paired, no specialty.

The risk here is genuinely lower than in #374 (that change altered reader construction *per path*;
this one changes two pure functions every path consumes), and I confirmed the paths by hand:
`--paired P_R1.FASTQ.GZ P_R2.FASTQ.GZ` → `P_R1_val_1.fq.gz` / `P_R2_val_2.fq.gz`, and a mixed-case
pair (`MIX_R1.fastq.gz` + `MIX_R2.FASTQ.GZ`) → `MIX_R1_val_1.fq.gz` / `MIX_R2_val_2.fq.gz`. Both
correct. So this is a coverage-convention gap, not a suspected bug.

**Recommend:** one `--paired` uppercase case (~8 lines, reusing `write_gz`/`sample_reads`), which
restores the file's stated matrix and pins the `input[0]`-drives-the-run-wide-`gzip` coupling that
paired mode adds. A specialty case is optional — note that `--hardtrim5/--clock` write into the
CWD, so such a test must set `current_dir`.

---

## Low

**L1 — the V4 fixture is silently lossy on a case-insensitive filesystem.** *Verified on this APFS
volume:* `sample.fq.gz` and `SAMPLE.FQ.GZ` are **one file**. The third `write_gz` reopens it, the
`"B"` body is overwritten by `"C"`, and the directory entry keeps the first-created spelling
(`ls` shows two files, not three). Consequences: the loop's second iteration feeds the same inode
under a different spelling, so on macOS it fires the *duplicate-output* branch
(`out/SAMPLE_trimmed.fq.gz and out/SAMPLE_trimmed.fq.gz`) rather than testing two distinct inputs;
on Linux all three files exist and it tests what the docstring says. **Both iterations still refuse
on both filesystems** (verified), so the test as written cannot be fooled — but the next person to
add a content assertion will be. One comment noting the APFS aliasing is enough.

**L2 — the newly-reachable refusals include the output-vs-input branch, which §3.1 and the
CHANGELOG do not mention.** *Verified:* `trim_galore -o . A.FASTQ.GZ a_trimmed.fq.gz` now fails
with "…would write output to `a_trimmed.fq.gz`, which is also one of its inputs"; before the fold
the planned output was `A.FASTQ_trimmed.fq` and the run was accepted. Correct on APFS (they alias);
a case-folded false positive on a case-sensitive filesystem — the trade-off `norm_path`'s
doc-comment already accepts. Worth half a clause in §3.1 since all three of its "newly rejected"
rows are output-vs-output.

**L3 — `is_gzipped`'s new "must stay in step with `strip_fastq_extensions`" comment overstates what
the code enforces.** The two use different matchers (`Path::extension()` vs byte-suffix strip) and
already disagree on dotfile names: `strip_fastq_extensions(".GZ")` treats `.GZ` as a gz suffix and
returns `""`, while `is_gzipped(".GZ")` is false because `extension()` is `None`. *Verified:* a file
named `.FQ.GZ` yields `out/_trimmed.fq.gz` (empty stem) where it used to yield `.FQ_trimmed.fq`.
Pathological, and identical in shape for lowercase `.gz` today — so not a regression — but the
comment now asserts a lockstep that holds only for the non-degenerate cases (and explicitly not for
`.bgz`). Recommend narrowing it to "the `gz` decision must agree with the `gz`-family strip".

**L4 — V8's substring assertions are matched against a report that embeds the full temp path.**
`report.contains("gzip")` / `!report.contains("plain")` would false-fail if `std::env::temp_dir()`
ever contained either word. The discrimination is carried by `!contains("plain")` and
`Compression ratio:` (`contains("gzip")` is satisfied by the `Output:` label too). Anchoring on the
line — e.g. `report.contains("(gzip, ")` — is both tighter and immune.

**L5 — `assert_eq!(&bytes[..2], …)` in V5 panics with a slice-index message, not the assert's, if
the output is <2 bytes.** Can't happen today (`clump_only.rs:343` and the gzip writer always emit a
member), but `bytes.starts_with(&[0x1f, 0x8b])` reads the same and fails legibly.

**L6 — the plan's V2 asks for a lowercase `.bgz` false case; nothing asserts it.** The pre-existing
`test_is_gzipped` (io.rs:965-977) never covered `.bgz`; the new test covers `.BGZ` only. Same code
path, so this is a one-line completeness item, not a hole.

**L7 — a reader primed by #381/#382 could misread the CHANGELOG on `.BGZ`.** Bullet 1 names
`.gz`/`.bgz`/`.bgzf` as the folded set; bullet 2 says `.GZ` output is now gzip-compressed. `.BGZ`
output is still plain (*verified*), matching lowercase `.bgz`. Half a sentence closes it.

---

## Efficiency

Nothing to raise. Five `eq_ignore_ascii_case` probes per input filename, once per file, no
allocation beyond the pre-existing `name: String`. `is_gzipped` swaps one `OsStr` compare for
another. Unmeasurable, as §6 says.

## Structure

The helper sits beside its only consumer, mirrors `str::strip_suffix` + `eq_ignore_ascii_case` in
name and signature, and is correctly private. Comments are one line each per house style, with the
non-obvious part (why byte-wise) in the doc-comment where it belongs. `.then(||…)` rather than
`then_some(…)` is not a style choice — `then_some` would evaluate the slice unconditionally and
reintroduce the panic. Test names are declarative and each carries the issue number.
