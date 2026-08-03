# CODE REVIEW — Reviewer B

**Target:** uncommitted working tree on `fix/phred64-ubam` (base `dev` @ `2e3da11`)
**Plan:** `plans/phred64-ubam/PLAN.md` (v2), implementation notes in §13
**Issue:** [#358](https://github.com/FelixKrueger/TrimGalore/issues/358)
**Scope:** 9 files, +517/−14, plus new `test_files/phred64_test.fastq`
**Verdict: APPROVE WITH REVISIONS**

No fixes were applied. Per the project's `CLAUDE.md`, source edits require an explicit
implementation trigger, so every item below is a recommendation.

---

## Summary

The core fix is correct and I could not break it. `input_phred_offset` is threaded to all
7 production `BamWriter::create` sites; the Bug 2 guard is sited ahead of every dispatch
branch and I found **no bypass** across 11 mode combinations; the ordering constraint in
§5.1 was honoured. `fmt` clean, `clippy --all-targets --release -D warnings` clean,
`cargo test --release` = **418 passed / 0 failed**, matching §13.

The most important thing I was asked to check — whether the fixture deviation papers over a
real defect — **it does not**. I reproduced the FastQC finding independently and traced it
to the correct cause (below). The output BAM bytes are provably right.

What does need work is almost entirely **user-facing accuracy**, and it is not cosmetic: two
CHANGELOG claims are demonstrably false (one of them tells users their archived data is
unrecoverable when it is recoverable in two commands), and the "all-zero-quality" framing
repeated across four strings is wrong for real modern data in a direction that *disarms* the
mitigation it justifies. Plus one code change: the `saturating_add` hardening trades a
debug-only panic for a release-mode invariant violation.

---

## The fixture deviation (PLAN §13) — independently verified, root cause confirmed

The implementer changed a test artefact in response to a failing check. I treated this as the
primary suspect and attacked it directly.

**Reproduced the reported symptom.** Uniform-`'h'` fixture, `--phred64 --output-format ubam
--fastqc -q 0 --length 0`:

```
raw QUAL bytes in BAM : I I I I …          (ASCII 73 → raw 40, CORRECT)
fastqc_data.txt       : Encoding  Illumina 1.5
                        base 1 mean  9.0
```

**The BAM is right; FastQC is wrong.** `samtools view` renders `QUAL` as `ASCII(raw+33)`, and
it renders `I` (73) → raw 40. The written bytes are spec-correct. FastQC then re-adds 33
(→ 73), finds no byte below 64, and applies its offset heuristic → subtracts 64 → 9.

**It is not BAM-specific, and it is not a `fastqc-rust` deviation.** The same input to the
**FASTQ** output path — a path this change does not touch — produces the identical
misdetection (`Encoding Illumina 1.5`, mean 9.0). And the lowest-quality-char rule is
inherited FastQC 0.12.1 behaviour, not something `fastqc-rust` invented.

**Conclusion:** genuine fixture degeneracy. Real data always carries some sub-Q31 base, and
the mixed fixture (min Q2) is strictly better — it exercises two offsets rather than one, and
mirrors real Illumina 1.5 `B`-run tails. The deviation is well-judged and correctly recorded.
Two wording nits follow from the above (L-2, L-7).

---

## High

### H-1 — `saturating_add` breaks the `qual.len() == seq.len()` byte invariant (release)

`src/bam.rs:962` — `.map(|&b| b.saturating_add(PHRED_OFFSET) as char)` collected into a
`String`.

`255u8 as char` is U+00FF, which is **2 bytes in UTF-8**. The codebase treats
`FastqRecord.qual` as one byte per base and slices it with byte indices derived from
`seq.len()`. Verified by execution (`rustc -O`, 5 bases with one `0xFF`):

```
seq bases = 5   qual String byte len = 6   char count = 5
bytes = [63, 63, 195, 191, 63, 63]

truncate(4)  -> OK but keeps [63,63,195,191]  = 3 bases of quality   (silent misalignment)
s[3..]       -> PANIC: "start byte index 3 is not a char boundary; it is inside 'ÿ'"
pre-change (wrapping) byte len = 5  bytes = [63, 63, 32, 63, 63]     (invariant intact)
```

`FastqRecord::truncate` (`src/fastq.rs:75`) calls `String::truncate`; `clip_5prime` /
`clip_3prime` (`:84`, `:97`) index the `String` directly. Both are on the hot trimming path.
Downstream consequences in release builds, where the neighbouring `debug_assert_eq!`s are
compiled out:

- `clip_5prime`/`clip_3prime` → **release panic** at a call site far from the cause, with an
  opaque char-boundary message.
- `truncate` → silently keeps the wrong number of quality bases.
- FASTQ output writes `qual.as_bytes()` → quality line longer than sequence line →
  structurally invalid FASTQ.
- `quality.rs` indexes `as_bytes()` in lockstep with `seq` → misaligned past the first `0xFF`.

For the sentinel byte specifically this is a **regression in release**: the previous `+`
wrapped to `0x20`, a single ASCII byte that kept the invariant (wrong value, intact
structure), and panicked loudly in debug. The change removes the debug panic and folds the
sentinel into a silent invariant violation.

Fairness note: the multi-byte class is **partly pre-existing** — any raw byte ≥ 95 already
produced a multi-byte char (`95+33 = 128` → U+0080, 2 bytes). The change widens it to
223–255 rather than creating it. But its stated purpose was to harden this exact input class,
and within that class it makes things worse.

The in-code comment asserts the residual is benign — "the point here is only to remove the
overflow" — which is the part that is not true. Fix stays inside the stated remit and is one
token:

```rust
.map(|&b| (b.min(93) + PHRED_OFFSET) as char)   // 93+33 = 126 = '~', always 1 ASCII byte
```

or, preserving sentinel semantics symmetrically with the all-`0xFF` branch above:

```rust
.map(|&b| if b == QUAL_MISSING_SENTINEL {
        QUAL_MISSING_REPLACEMENT as char
    } else { (b.min(93) + PHRED_OFFSET) as char })
```

Reachability is low (malformed BAM mixing sentinel and real scores), but so is the
reachability of the overflow the change was added to remove — they are the same input.

### H-2 — CHANGELOG's FastQC claim is false

> "`--fastqc` on uBAM output was affected too — per-base-quality plots were centred on Q71."

The opposite is true. Pre-fix, stored raw = `ASCII_in − 33`, so the bytes FastQC reads back
(`raw + 33`) are **exactly the original input's ASCII bytes**. For any legal Phred+64 input
those are all ≥ 64, so FastQC's heuristic picks Illumina 1.5 and subtracts 64 — recovering
the true score. Two errors cancelled.

Demonstrated on the project's own fixture. Case B is byte-identical to pre-fix output
(produced by running the same fixture *without* `--phred64`):

| | raw QUAL | FastQC Encoding | base 1 | base 31 |
|---|---|---|---|---|
| A: post-fix | `[40×24, 2×10]` | Sanger / Illumina 1.9 | **40.0** | **2.0** |
| B: pre-fix bytes | `[71×24, 33×10]` | Illumina 1.5 | **40.0** | **2.0** |

Identical numbers. PLAN §13 already records this correctly ("pre-fix it read ≈40 *by
accident*"); the CHANGELOG carries the pre-§13 claim forward. Delete the sentence, or replace
it with the accurate version: pre-fix FastQC numbers were coincidentally right, and post-fix
they are right for the right reason.

### H-3 — CHANGELOG's remediation advice is wrong; archived output IS recoverable

> **"Archived output cannot be repaired by re-running Trim Galore."** … "Regenerate from the
> original FASTQ."

Verified recoverable in two commands, with no access to the original FASTQ:

```
samtools fastq prefix.bam > p64.fq          # emits ASCII = true+64 — a genuine Phred+64 FASTQ
trim_galore --phred64 --output-format ubam p64.fq
```

Resulting QUAL is byte-identical to the canonical post-fix output (`[40×24, 2×10]`), which I
confirmed against case A above. This is exact for **any** Phred+64 input, because
`samtools fastq` reproduces the original input ASCII bytes verbatim.

The reasoning given in the entry is self-defeating: "reading it back adds 33 (yielding Q71),
and `--phred64` would subtract 64 rather than 31." Reading back yields ASCII 104; `--phred64`
declares offset 64; `104 − 64 = 40` — exactly right. The Q71 figure only appears if you
interpret the recovered ASCII at offset 33, which is precisely what `--phred64` tells the tool
not to do.

Only the **single-invocation** repair is blocked, by the new guard. Telling users to
regenerate from source data they may no longer have is a materially harmful inaccuracy in a
release note. Reword to: "cannot be repaired in one invocation (`--phred64` is now rejected
for BAM input), but the data is recoverable — round-trip through `samtools fastq`, which
yields genuine Phred+64 ASCII, then re-run with `--phred64`."

---

## Medium

### M-1 — "all-zero-quality BAM" is false for real data, and the truth is more dangerous

Four user-facing strings assert it:

- `src/main.rs:299-302` — O-3 startup NOTE: "running Phred+33 data with `--phred64` yields an
  all-zero-quality BAM"
- `src/main.rs:216-217` — guard error text: "silently zeroing quality"
- `docs/src/content/docs/guide/quality.md` caution: "**every** quality score is zero, because
  the subtraction underflows and floors"
- `CHANGELOG.md` behaviour-change entry: "silently write an all-zero-quality BAM"

Measured on `test_files/clock_10K_R1.fastq.gz` (Phred+33, realistic modern profile) with
`--phred64 --output-format ubam -q 0 --length 0`:

```
total qual bytes : 729806
raw 0 (floored)  :  21565  =  3.0%
non-zero         : 708241  = 97.0%
top values       : Q7 ×604444, Q3 ×48659, Q6 ×24729, Q0 ×21565, …
```

**97% of bases are non-zero.** Only bytes below the offset floor; every base at Phred ≥ 31
becomes Phred − 31. Uniform-Q40 Phred+33 input produces a uniform-**Q9** BAM, not Q0.

This matters beyond pedantry. O-3's whole rationale is that the silent-floor case is
"indistinguishable downstream from legitimate Q0 data", and the NOTE hands the user a
detection heuristic — look for all-zero quality — that **will not fire**. The real artefact is
a plausible uniformly-low-quality dataset, which is *harder* to spot than an all-`!` file.
Suggested wording: "collapses quality toward zero — every score is reduced by 31 and floored
at 0, so the output looks like uniformly poor-quality data rather than obviously empty."

(For `ubam_test.bam` specifically the claim happens to hold — it is uniform raw 30, and
`30+33 = 63 < 64` floors to 0 — which is likely why it survived testing.)

### M-2 — `phred64_mixed_input_rejected` is a worthless guard (verified passes pre-fix)

`tests/integration_ubam_out.rs` — the test asserts only `!output.status.success()`. Executed
the same invocation **without** `--phred64`:

```
$ trim_galore --paired test_files/phred64_test.fastq test_files/ubam_test.bam
EXIT=1
Error: --paired with two BAM files is not supported. …
```

So the pre-fix binary also exits non-zero and the assertion holds. The test cannot distinguish
the guard from an unrelated pre-existing rejection. Add the assertion its four sibling guard
tests already use:

```rust
assert!(String::from_utf8_lossy(&output.stderr).contains("--phred64"), …);
```

(Tangentially: that pre-existing error message is itself wrong for a mixed pair — only one of
the two inputs is a BAM. Pre-existing, out of scope, worth a follow-up.)

### M-3 — `phred64_ubam_out_specialty_stores_true_phred` passes vacuously on an empty BAM

```rust
for q in &quals(&out) {
    assert_eq!(q.len(), 20, …);
    assert!(q.iter().all(|&b| b == 40), …);
}
```

No emptiness check. A header-only BAM satisfies this test. Both siblings guard against it
(`assert!(!qs.is_empty())` in the SE test, `assert!(!tuples.is_empty())` in the clump-only
test). Add the same line. The `q.len() == 20` assertion is fine given the value check beside
it, so this is the only structural weakness in the specialty test.

### M-4 — CHANGELOG's list of removed working invocations is presented as exhaustive but is not

> "This does remove **two** invocations that previously worked, because `--clump_only` and
> `--hardtrim5/3` perform no quality arithmetic…"

`-q 0` also makes quality trimming inert on the ordinary trim path.
`quality_trim_3prime` (`src/quality.rs:22`) with `cutoff = 0` computes
`running_sum += 0 − q`, which is ≤ 0 for all bases, so `max_sum` never exceeds 0 and nothing
is trimmed. Therefore `trim_galore --phred64 -q 0 <input.bam>` also produced correct output
pre-fix and is now rejected — a third case, on the mode most users are in.

Generalise: "any BAM-input invocation in which quality trimming was inert — the specialty and
clump-only modes, and `-q 0` on the trim path."

### M-5 — 4 of 7 threaded writer sites have no test

| Site | Covered? |
|---|---|
`main.rs:1852` trim SE | ✅ `phred64_ubam_out_se_stores_true_phred` |
`specialty.rs:133` hardtrim5 | ✅ `phred64_ubam_out_specialty_stores_true_phred` |
`clump_only.rs:772` clump-only SE | ✅ `phred64_clump_only_ubam_out_stores_true_phred` |
`main.rs:1994` PE two-file FASTQ | ❌ |
`specialty.rs:183` hardtrim3 | ❌ |
`clump_only.rs:989` clump-only PE | ❌ |
`main.rs:2126` PE interleaved | ❌ (BAM-only input → unreachable with `--phred64`; correctly untestable) |

I exercised all three reachable uncovered sites by hand and every one is correct:

```
--paired --phred64 --output-format ubam r1.fastq r2.fastq  -> IIIIIIII…##########  (raw 40 / 2) ✓
--hardtrim3 20 --phred64 --output-format ubam              -> IIIIIIIIII##########             ✓
--clump_only --paired --phred64 --output-format ubam       -> IIIIIIII…##########              ✓
```

So this is coverage debt, not a live bug. PE two-file FASTQ + `--phred64 --output-format ubam`
is a realistic user invocation and the plan's own §9 framing ("the bug class is faulty
*wiring*, and a unit test pins only arithmetic") argues for pinning it.

Also uncovered: `--phred64` with the **default** `-q 20` on the uBAM path. Verified correct
(B-run trimmed, 24 bases at raw 40 retained), and it is the only case that exercises
phred64 quality-trimming *and* uBAM writing together — currently every test disables one or
the other via `-q 0`.

### M-6 — the guard→offset coupling is invisible at the 7 call sites

Every site passes `cli.phred_offset()` unconditionally. Correctness on BAM-input paths depends
entirely on a `bail!` at `main.rs:211`, up to ~1900 lines away. Nothing at
`clump_only.rs:989` or `specialty.rs:183` hints that passing 64 there would silently zero
quality.

A single local right after the guard, threaded to all sites, makes the dependency greppable
for one line of churn:

```rust
// Safe to use on BAM-input paths only because of the --phred64 guard above.
let input_phred_offset = cli.phred_offset();
```

This targets exactly the failure mode the required-parameter design (§4) was chosen to
prevent: "a site that never considered the offset."

---

## Low

- **L-1** `src/cli.rs:143` — `--phred64` help text unchanged ("Use Phred+64 quality encoding
  (Illumina 1.5). Default is Phred+33."), no mention of the hard rejection for BAM input. The
  file has an established convention here: `OutputFormat::UBam`'s doc enumerates its rejected
  flags, and `--clump_only`'s does the same.
- **L-2** `test_files/README.md` frames the FastQC misdetection as affecting "the round-tripped
  BAM" and as a `fastqc-rust` quirk. Neither is precise — I reproduced it identically on the
  FASTQ output path, and the lowest-char rule is upstream FastQC 0.12.1 behaviour. Say
  "FastQC's encoding heuristic (upstream behaviour, on either output format)".
- **L-3** `test_files/README.md` closing paragraph is backwards: it says `truncated.fq.gz` "is
  easy to misclassify because a naive minimum-quality-byte heuristic reads it as high-quality
  Phred+33 rather than moderate-quality Phred+64." Verified min ASCII 66 / max 99 — the
  min-byte heuristic classifies it **correctly** as Phred+64 (66 ≥ 64). What misreads it is an
  unconditional Phred+33 assumption. As written the paragraph contradicts the table row above
  it, in a file whose entire purpose is preventing this class of error.
- **L-4** `#### Behaviour changes` is a new third CHANGELOG heading. The Unreleased section's
  existing convention files rejections under `#### Changes` (e.g. "`--dont_gzip` +
  `--output-format ubam` now rejected at CLI validation"). Consider folding in.
- **L-5** PLAN §13 states "+500/−14"; actual `git diff --stat` is **+517/−14**.
- **L-6** `docs/guide/quality.md` now states "FASTQ output preserves the input encoding, so a
  Phred+64 run produces Phred+64 output." True for FASTQ→FASTQ, but a Phred+64 → uBAM →
  re-read → FASTQ round trip emits **Phred+33**, because the BAM boundary discards the ASCII
  encoding by design. Correct behaviour; one sentence would stop it being filed as a bug.
- **L-7** `--fastqc` on uBAM output under-reports quality by 31 for any dataset whose minimum
  Phred ≥ 31 (uniform-Q40 → mean 9.0, reproduced above). Pre-existing, upstream-faithful, and
  **not introduced here** — but for Phred+64 data specifically this fix moves the case from
  accidentally-correct to affected. For BAM sources the heuristic is never appropriate (BAM
  `QUAL` is raw Phred by spec), so an upstream issue against `ewels/FastQC-Rust` asking it to
  force Sanger for BAM input would close it cleanly. Out of scope for this PR; worth filing
  since this investigation surfaced it.

---

## Attacked and could not break

**Guard bypass — 11 combinations, no escape.** `bail!` fires for plain trim, `--rrbs`,
`--hardtrim5`, `--hardtrim3`, `--clock` (even input), `--implicon` (even input),
`--clump_only`, `--clumpify --cores 2`, `--demux`, `--paired` interleaved BAM,
`--output-format ubam`, and mixed FASTQ+BAM. `--clock`/`--implicon` with an odd file count
fail earlier in `cli.validate()` — correct precedence, not a bypass. `main.rs` is the only
`[[bin]]`, and nothing dispatches or returns between `cli.validate()` and the guard.

**`sanity_check_any` before the guard.** It *does* interpret quality (`BamReader::open` →
`next_record` → `bam_record_to_fastq`), but only via the fixed read-side `PHRED_OFFSET`; it
never consults `cli.phred_offset()`, so it cannot misinterpret. No output is produced before
the bail. Guard rejection creates **no output directory** — verified (`ensure_output_dir` is
after it). Plan check 5 asserted this; no test does, since `fresh_tmpdir` pre-creates the dir.

**O-3 NOTE reachability — B-NIT-2 not repeated.** Fires on trim SE, `--hardtrim5`, and
`--clump_only` uBAM paths (all verified by execution), and correctly does *not* fire for
FASTQ output. Sited before dispatch, so the early-returning specialty paths see it.

**Edge cases.** Read fully consumed by adapter trimming under `--length 0` → `*`/`*`, no
panic. Solexa/Illumina-1.0 bytes (ASCII 59–63) under `--phred64` → floor to raw 0, exit 0, as
§3.2 predicts. Sub-offset byte (ASCII 32 at offset 33) → floors to raw 0 — `saturating_sub` is
confirmed load-bearing, not decorative. Seq/qual length mismatch → **hard error in release**
from noodles ("sequence-quality scores length mismatch: expected 34, got 20"), so the
compiled-out `debug_assert_eq!` is not load-bearing, and nothing in the diff cites it as
coverage (checked).

**Round-trip / stranded workflows.** No write→read asymmetry introduced: post-fix output is
true Phred and re-reads correctly with no flag. Nothing legitimate is stranded — a genuinely
Phred+64 BAM is spec-illegal, and for the only real affected class (pre-fix TrimGalore output,
raw = true+31) `--phred64` was never the right correction, since it needs −31. Recovery path
exists and works (H-3).

**Fixture claims.** `truncated.fq.gz` Phred+64 verified: min ASCII 66 / max 99 → Q2..Q35 at
offset 64, `HWUSI-EAS611_0001:…#0/1` read names, trailing `B` runs. `test_files/README.md`
does **not** assert all fixtures are Phred+33 — it explicitly forbids it. Relative
`test_files/…` paths in `integration_ubam_out.rs` match that file's existing convention (26
prior uses; no `fixture()` helper there).

**Gates.** `cargo fmt --all -- --check` clean. `cargo clippy --all-targets --release
-- -D warnings` clean. `cargo test --release` = **418 passed, 0 failed** (359 lib + 59
integration), matching §13 exactly. All 9 new tests observed running and passing.

**Pre-fix failure analysis.** 8 of 9 new tests genuinely fail pre-fix (the unit test cannot
compile pre-fix, which is stronger). The one exception is M-2. The Phred+33 counterpart
assertion inside the unit test is correctly self-labelled as an over-correction guard rather
than a bug guard.

---

## Recommended order

1. **H-1** — one-token code change; the only source edit needed.
2. **H-2, H-3, M-1, M-4** — CHANGELOG + `main.rs` NOTE/error text + `docs/guide/quality.md`.
   These are what users read; two are currently false and one is actively misleading advice.
3. **M-2, M-3** — two added assertions; both are one line.
4. **M-5, M-6, L-1** — coverage and maintainability; cheap, and aligned with the plan's own
   stated design intent.
5. **L-2 … L-7** — wording and follow-up filings.

Nothing here questions the fix's correctness or the decision to ship it as one PR. The §5.1
coupling argument holds and the implementation respected it.
