# Code review — PR #382 `fix(io): strip .bgz/.bgzf so output and report share a stem (#381)`

**Reviewer A** · base `dev` (`c6528e5`) · head `ae55301` · author `BenjaminDEMAILLE`
Files: `src/io.rs`, `tests/integration_gzip_non_gz_extension.rs`, `CHANGELOG.md` (+165 / −35)

## Recommendation

**MERGE.** The rewrite is behaviour-preserving for every name that does not end in `.bgz`/`.bgzf` — proven exhaustively and analytically below, not sampled. The byte-identity invariant is untouched. All three CI gates pass locally. The recommended fixes are all documentation-accuracy and test-coverage items; none block.

## Summary

The change is small, correct, and the diagnosis in #381 is accurate. The rewrite is a genuine improvement on the fix proposed in the issue: the issue suggested extending the enumerated combined-suffix list to `[".fastq.gz", ".fq.gz", ".fastq.bgz", ".fq.bgz", ".fastq", ".fq"]`, which grows multiplicatively and would have silently left the bare-`.bgz` case to the fallback. Two sequential strips give a cleaner invariant, which is the strongest thing that can be said about this diff:

> **For every name `X`: `new(X + ".bgz") == new(X + ".bgzf") == old(X + ".gz")`, and `new(X + ".gz") == old(X + ".gz")`.**

That is, `.bgz`/`.bgzf` now behave *exactly* as `.gz` always did, and nothing else moves at all. I verified this over 2,352 generated assertions with 0 violations, and separately over 28,293 paths / 9,431 unique basenames with **0 differing basenames that do not end in `.bgz`/`.bgzf`**.

**Verification performed** (all from crate root, `pr382` checked out):

| Gate | Result |
|---|---|
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets --release -- -D warnings` | clean, exit 0 |
| `cargo test` | exit 0 — 383 lib tests + all integration suites, 0 failed |
| `cargo test strip_fastq_extensions` | 3 passed |
| `cargo test --test integration_gzip_non_gz_extension` | 7 passed |

---

## 1. The equivalence table (load-bearing evidence)

Old implementation = `dev:src/io.rs:435`; new = `src/io.rs:451`. I reimplemented both verbatim in a standalone harness and compared. `old` = `dev`, `new` = `pr382`.

| input basename | old (`dev`) | new (`pr382`) | | note |
|---|---|---|---|---|
| `sample.fastq.gz` | `sample` | `sample` | same | |
| `sample.fq.gz` | `sample` | `sample` | same | |
| `sample.fastq` | `sample` | `sample` | same | |
| `sample.fq` | `sample` | `sample` | same | |
| `sample.gz` | `sample` | `sample` | same | |
| `sample.txt.gz` | `sample` | `sample` | same | double-strip fallback |
| `sample.bam` | `sample` | `sample` | same | **uBAM path unaffected** |
| `sample.ubam` | `sample` | `sample` | same | |
| `sample` | `sample` | `sample` | same | |
| `sample.fastq.bz2` | `sample.fastq` | `sample.fastq` | same | pre-existing wart, preserved |
| `archive.fq.tar.gz` | `archive.fq` | `archive.fq` | same | pre-existing wart, preserved |
| `SAMPLE.FASTQ.GZ` | `SAMPLE.FASTQ` | `SAMPLE.FASTQ` | same | case-sensitivity limitation, preserved |
| `sample.fastq.gz.gz` | `sample.fastq` | `sample.fastq` | same | |
| `sample.fq.bgz.gz` | `sample.fq` | `sample.fq` | same | outer `.gz` wins in both |
| `.fq` | `""` | `""` | same | degenerate, pre-existing |
| `""` (empty) | `""` | `""` | same | `file_name()` → `None` in both |
| **`sample.fastq.bgz`** | `sample.fastq` | `sample` | **DIFF** | **the #381 fix** |
| **`sample.fq.bgz`** | `sample.fq` | `sample` | **DIFF** | the #381 fix |
| **`sample.fastq.bgzf`** | `sample.fastq` | `sample` | **DIFF** | the #381 fix |
| **`sample.fq.bgzf`** | `sample.fq` | `sample` | **DIFF** | the #381 fix |
| **`lib.2_R1.fastq.bgz`** | `lib.2_R1.fastq` | `lib.2_R1` | **DIFF** | the #381 fix, dotted stem preserved correctly |
| `sample.bgz` | `sample` | `sample` | same | bare form, pinned by the PR's test |
| `sample.bgzf` | `sample` | `sample` | same | bare form, pinned |
| **`data.tar.bgz`** | `data.tar` | `data` | **DIFF** | ⚠️ class (a) below — *not* mentioned in the PR |
| **`sample.bam.bgz`** | `sample.bam` | `sample` | **DIFF** | ⚠️ class (a) below |
| **`sample.FASTQ.bgz`** | `sample.FASTQ` | `sample` | **DIFF** | class (a); incidentally tidier |
| **`.bgz`** | `.bgz` | `""` | **DIFF** | ⚠️ class (b) below — degenerate dotfile |
| **`.fq.bgz`** | `.fq` | `""` | **DIFF** | class (b) |

**Analytic proof of the "same" rows** (so the result does not depend on my corpus choice). Let `N` be any name not ending in `.bgz`/`.bgzf`. Then the new gzip-strip can only match `.gz`, so:

- *Case A — `N` ends in `.gz`.* New sets `s = N[..-3]`. If `s` ends `.fastq`, new returns `s[..-6]`, which is exactly old's branch 1 (`.fastq.gz`). If `s` ends `.fq`, new returns `s[..-3]`, exactly old's branch 2. Otherwise new returns `file_stem(s)`; old cannot match branches 1–4 (it ends `.gz`, so not `.fastq`/`.fq`; and `s` ends with neither `.fastq` nor `.fq`, so not `.fastq.gz`/`.fq.gz`) and its fallback strips `.gz` to `s` and returns `file_stem(s)`. Identical.
- *Case B — `N` does not end in `.gz`.* New sets `s = N`; old cannot match branches 1–2, and its fallback leaves `N` untouched. Both then strip `.fastq`, else `.fq`, else `file_stem(N)`. Identical.

The branch orders are non-load-bearing in both implementations because no member of `{.fastq, .fq}` is a suffix of another, and likewise for `{.gz, .bgz, .bgzf}` — so the doc comment's disjointness claim at `src/io.rs:458-459` is **correct**.

### The two unadvertised difference classes

Both are confined to names ending `.bgz`/`.bgzf`, and both are *consistent with what `.gz` already did*, so neither is a defect. They are only a documentation problem.

**(a) Multi-component names lose one more component than "the inner `.fastq`".** `data.tar.bgz` → `data` (was `data.tar`); `sample.bam.bgz` → `sample` (was `sample.bam`). This is because new runs `file_stem` on the *`.bgz`-stripped* name, whereas old ran it on the full name. It is exactly `old("data.tar.gz") == "data"`, so it is the right answer — but it contradicts the claim in `CHANGELOG.md:46-48` and in the PR body that "unrelated extensions keep the stems they had". See Finding 1.

**(b) Degenerate dotfiles yield an empty stem.** A file literally named `.bgz` now gives `""` (was `.bgz`), so the output would be `_trimmed.fq`. This matches `old(".gz") == ""` and `old(".fq") == ""` — the degeneracy already existed for `.gz`/`.fq` and is not worth code. Informational only.

## 2. Blast radius

`strip_fastq_extensions` has **14 call sites** (`src/io.rs:90,116,143,179,180,210,211,245,283,315,316,343,369` and `src/specialty.rs:430,452,466,482`). Every one is a filename constructor. For a `.bgz`/`.bgzf` input the following user-visible names therefore move, and the PR mentions only the first:

`single_end_output_name` (`_trimmed`), `paired_end_output_names` (`_val_{1,2}`), `unpaired_output_names` (`_unpaired_{1,2}`), `passthrough_output_name` (`_passthrough`), `clumped_output_name` / `clumped_paired_output_names` (`_clumped`), `single_end_bam_output_name` / `paired_bam_output_name` / `clumped_bam_output_name` / `clumped_paired_bam_output_name` (`.bam`), `hardtrim_output_name` / `hardtrim_bam_output_name`, `clock_output_name`, `implicon_output_name`.

**uBAM is safe.** `sample.bam` → `sample` in both implementations (table above), so `single_end_bam_output_name`, `paired_bam_output_name` and `clumped_bam_output_name` are unchanged for their real inputs. The existing pins at `src/io.rs:806-816` still hold, and the comment at `src/io.rs:814` ("falls through to `file_stem()` for `.bam`") remains accurate.

**`--implicon` changes by more than one component**, which nothing in the PR predicts. `implicon_output_name` (`src/specialty.rs:482-500`) strips a trailing `_R1`/`R1` from the stem to avoid a doubled `R`. On `dev` the `.bgz` stem was `sample_R1.fastq`, which ends in neither, so the de-duplication silently failed. Verified empirically with the built binary:

```
# pr382, input sample_R{1,2}.fastq.bgz
sample_8bp_UMI_R1.fastq      sample_8bp_UMI_R2.fastq
# pr382, .fastq.gz control — identical stems, as the invariant predicts
sample_8bp_UMI_R1.fastq.gz   sample_8bp_UMI_R2.fastq.gz
# dev would have produced (derived from src/specialty.rs:482-500):
sample_R1.fastq_8bp_UMI_R1.fastq
```

This is the correct outcome — but it is a second-order repair, worth a CHANGELOG clause (Finding 2).

**`--demux` inherits the fix for free.** `demux.rs:129-141` derives its base name from the *already-trimmed* filename with its own ad-hoc `.gz`-then-`.fq` strip, not from `strip_fastq_extensions`, so it picks up the tidier stem automatically. No action.

**Output-collision exposure widens (pre-existing hole, not a regression).** Two SE inputs differing only in `.fastq` vs `.fq` now share one stem. Verified with the built binary:

```
$ trim_galore -o out coll.fastq.bgz coll.fq.bgz
out/coll_trimmed.fq                    ← one output file
out/coll.fastq.bgz_trimming_report.txt ← two reports
out/coll.fq.bgz_trimming_report.txt
```

No error, no warning. I confirmed the identical silent overwrite already happens on this branch for `coll2.fastq.gz` + `coll2.fq.gz`, so this is **pre-existing**: the SE FASTQ loop (`src/main.rs:783-791`) has no output-collision pre-flight, unlike paired (`src/main.rs:673`), paired-BAM (`src/main.rs:1799`), SE-BAM (`src/main.rs:609`) and specialty multi-pair (`src/main.rs:2410`). This PR widens the exposed name set by one class without creating the hole. Worth its own issue; **not** a reason to hold #382.

## 3. Byte-identity (both PR claims verified, not accepted)

- **"Nothing in `test_files/` is named `.bgz`"** — verified: `ls test_files/ | grep -Ei 'bgz|bgzf'` returns nothing.
- **"Nothing in the `validation` matrix is named `.bgz`"** — verified: `grep -rn bgz .github/workflows/` returns **no hits at all**. Every input in the validation job is `.fastq.gz` (`illumina_10K.fastq.gz`, `BS-seq_10K_R{1,2}.fastq.gz`, `clock_10K_R{1,2}.fastq.gz`, `SRR24766921_RRBS_R{1,2}.fastq.gz`, and the `${sample}_R{1,2}.fastq.gz` copies).

Independently of those claims, the md5 baseline cannot move: the equivalence proof shows `.fastq.gz` stems are unchanged, and the report *content* is stem-independent — `report.rs:192`, `:550` and `:663` all interpolate `config.input_filename` (the input basename), never the derived stem. The only report body that embeds an output path is `write_passthrough_stats` (`report.rs:480-481`), whose fixture is `BS-seq_10K_I1.fastq.gz` — stem unchanged. **Byte-identity is safe.**

## 4. Test sufficiency

**Are the tests capable of failing?** Yes — checked case by case rather than assumed.

- `test_strip_fastq_extensions_bgz` (`src/io.rs:505`): the first assertion alone fails on `dev` (`old("sample.fastq.bgz") == "sample.fastq"` vs expected `"sample"`), matching the failure the PR quotes. Five of its seven assertions fail on `dev`; the two bare-form assertions (`sample.bgz`, `sample.bgzf`) pass on both and are correctly described in-line as pins rather than regressions.
- `test_strip_fastq_extensions_leaves_other_names_alone` (`src/io.rs:534`): all four assertions pass on `dev` too. It cannot fail as a #381 regression test — it is a forward guard against the *next* rewrite. That is its stated purpose, so no objection.
- `bgz_output_and_report_agree_on_the_stem` (`tests/…:305`): asserts `sample_trimmed.fq` exists, which `dev` does not produce (it writes `sample.fastq_trimmed.fq`), so it fails on `dev` with the message the PR quotes.
- **The rewritten assertions in the author's own #374 tests are correct and not weakened.** `read_output_raw` (`tests/…:94-107`) `panic!`s when neither `{stem}.fq` nor `{stem}.fq.gz` exists, so changing `read_output(&dir, "pair_R1.fq_val_1")` → `read_output(&dir, "pair_R1_val_1")` is a real assertion that would fail on `dev`, not a relaxation. The direction is right: `pair_R1_val_1` is what `.fq.gz` has always produced. `read_output`'s non-empty check (`tests/…:84-89`) is untouched.

**What is left unpinned.** Three rows of the table above are not pinned by any test, and they are the rows most likely to be broken by a future rewrite of this function:

1. **A dotted stem** (`lib.2_R1.fastq.bgz` → `lib.2_R1`). This is the case that distinguishes suffix-stripping from `Path`-component semantics, and `sample.L001.fastq.bgz` is a realistic Illumina name.
2. **`archive.fq.tar.gz` → `archive.fq`** — the double-strip fallback, currently unpinned in either test.
3. **`SAMPLE.FASTQ.GZ` → `SAMPLE.FASTQ`** — see Finding 5.

Also unpinned: the `--clump_only` stem for `.bgz`. `clump_only_bgz_input_is_decompressed` (`tests/…:244`) asserts exit status only and never looks at the output name, which is why the PR did not need to touch it — but it means `clumped_output_name`'s changed `.bgz` behaviour rests on the unit test alone. Acceptable (same function), noted for completeness.

**Would a future rewrite be caught?** Mostly yes — the two unit tests plus the pre-existing `test_strip_fastq_extensions` pin 16 distinct inputs between them. Adding the three rows above would make the pin complete against the whole difference surface.

## 5. CHANGELOG placement

Correct. The new entry sits at `CHANGELOG.md:34`, inside `### Unreleased` → `#### Bug fixes` (opened at line 6, closed by `#### Changes` at line 51). That is the right section for a fix.

Good hygiene worth calling out: the PR also amends the #374 entry it supersedes, from "Two known limitations" down to "One known limitation remains" (`CHANGELOG.md:29-32`), rather than leaving a stale note claiming `.bgz` is not in the stripping list. I verified the surviving limitation is still true — `is_gzipped` (`src/io.rs:29-31`) tests `extension() == "gz"`, so `.bgz` input still yields plain `.fq` output.

---

## Fixes recommended

None applied — per instruction, no tracked file in this repository was modified. Exact edits given below.

### Finding 1 — CHANGELOG/PR claim "unrelated extensions keep the stems they had" is false — **Medium**

`CHANGELOG.md:46-48`. `data.tar.bgz` had stem `data.tar` and now has `data`; `sample.bam.bgz` had `sample.bam` and now has `sample`. Any `.bgz`-suffixed name with more than one preceding dot-component loses an extra component. The behaviour is right (it equals `.gz`), the sentence is not. Replacing it with the invariant I verified is both true and a stronger guarantee:

```diff
-  enumerated list of combined suffixes. Names with no inner `.fastq`/`.fq`
-  (`sample.bgz`), `.bam` inputs and unrelated extensions keep the stems they
-  had. Output *filenames* only: no fixture in `test_files/` or in the
+  enumerated list of combined suffixes. The effect is that `.bgz`/`.bgzf` now
+  behave exactly as `.gz` always has — `sample.bgz` → `sample`,
+  `data.tar.bgz` → `data`, just as `data.tar.gz` already did — and every name
+  not ending in `.bgz`/`.bgzf`, including `.bam` inputs and bare `.gz`, is
+  unchanged. Output *filenames* only: no fixture in `test_files/` or in the
   `validation` matrix is named `.bgz`, so the Perl-0.6.11 byte-identity
   comparison is untouched.
```

### Finding 2 — the specialty-mode and `--clump_only` name changes are unmentioned; `--implicon` moves by two components — **Medium**

Section 2 above. A maintainer reading "output filenames only, inner `.fastq` removed" would not predict `sample_R1.fastq_8bp_UMI_R1.fastq` → `sample_8bp_UMI_R1.fastq`. One clause after the sentence added by Finding 1:

```
  The same stem feeds `--clump_only`, `--hardtrim5/3`, `--clock` and
  `--implicon`, so their `.bgz` output names move too; `--implicon`'s `_R1`
  de-duplication now fires, which it could not when the stem ended in
  `.fastq`.
```

### Finding 3 — doc comments are four times longer than the convention allows, and the narrative is duplicated four times — **Low**

`CLAUDE.md` — "Default to a single line. Two is the maximum and needs a reason. State the fact, not the evidence… Detail that genuinely matters goes in the commit message, PR description, or a linked issue — not the source."

The #374/#381 narrative ("since #374 such a file is *read* correctly… the output and the report disagreeing about the sample name, which is what MultiQC groups on") appears **four times**: `src/io.rs:439-446`, `src/io.rs:498-503`, `tests/…:290-299`, and `CHANGELOG.md:34-42`. Only the last two earn it — the CHANGELOG is the user-facing record, and a test's docstring justifying *why the test exists* is the accepted exception. The function doc comment is the one to cut; paragraphs 1 and 3 are good API docs and should stay:

```diff
 /// Strip FASTQ extensions from a filename, returning just the base stem.
 ///
 /// Handles `.fastq` / `.fq`, each optionally followed by a gzip-family
 /// suffix: `.gz`, or the `.bgz` / `.bgzf` names `bgzip` output is often
 /// given. The two suffix groups are stripped in sequence rather than
 /// enumerated as a combined list, so adding one compression name does not
 /// multiply the cases.
 ///
-/// The `.bgz` forms matter because since
-/// [#374](https://github.com/FelixKrueger/TrimGalore/pull/374) such a file is
-/// *read* correctly (gzip-ness comes from the file's first three bytes, not
-/// its name). Before this stripped them too, `sample.fq.bgz` produced
-/// `sample.fq_trimmed.fq` alongside `sample.fq.bgz_trimming_report.txt`, the
-/// output and the report disagreeing about the sample name, which is what
-/// MultiQC groups on
-/// ([#381](https://github.com/FelixKrueger/TrimGalore/issues/381)).
+/// The `.bgz` forms have to be stripped here or the output and the report
+/// disagree about the sample name
+/// ([#381](https://github.com/FelixKrueger/TrimGalore/issues/381)).
 ///
 /// Anything else (a `.bam` input, an unrelated extension) falls through to
 /// `Path::file_stem`, which drops the final component only.
```

The two inline comments compress the same way — `src/io.rs:458-459` to `// One gzip-family suffix. No member is a suffix of another, so order is free.`, and `src/io.rs:465` ("Then the FASTQ extension itself, longest first.") to `// Then the FASTQ extension.`; see Finding 6 for why "longest first" is better dropped than kept.

### Finding 4 — add the three missing pins — **Low**

Section 4. Appended to `test_strip_fastq_extensions_leaves_other_names_alone` (`src/io.rs:534`), or as a third test:

```rust
// A dotted stem: the component boundaries must come from the suffixes, not
// from `Path`.
assert_eq!(
    strip_fastq_extensions(Path::new("lib.2_R1.fastq.bgz")),
    "lib.2_R1"
);
// The fallback drops one component only, so an inner `.fq` survives.
assert_eq!(
    strip_fastq_extensions(Path::new("archive.fq.tar.gz")),
    "archive.fq"
);
// Matching is case-sensitive; an uppercase name keeps its inner extension.
assert_eq!(
    strip_fastq_extensions(Path::new("SAMPLE.FASTQ.GZ")),
    "SAMPLE.FASTQ"
);
```

### Finding 5 — `.FASTQ.GZ` is the same bug as #381, still open — **Low, follow-up issue**

Matching is case-sensitive at `src/io.rs:461,466`, so `SAMPLE.FASTQ.GZ` → `SAMPLE.FASTQ` (both old and new). Since #374 made gzip-ness content-derived, such a file *reads fine* and then produces `SAMPLE.FASTQ_trimmed.fq` next to `SAMPLE.FASTQ.GZ_trimming_report.txt` — the exact #381 symptom, for the exact #381 reason (the extension list is a naming allowlist that the reader no longer consults). Not this PR's job, and case-folding the match is a behaviour change deserving its own discussion. Recommend a follow-up issue, plus the pin in Finding 4 so the current answer is at least recorded rather than accidental.

### Finding 6 — "longest first" implies an ordering constraint that does not exist — **Low**

`src/io.rs:465`. `.fastq` is not a suffix of `.fq`, so the order of that loop is free, exactly as the gzip list's is. A maintainer adding a third FASTQ spelling could rely on "longest first" as the rule that keeps things correct, when the actual rule is the one stated two comments earlier: *no member may be a suffix of another*. Either drop the rationale (Finding 3) or state the real one.

### Finding 7 — stale comment referring to the removed "strip list" — **Low**

`tests/integration_gzip_non_gz_extension.rs:287` still reads:

```rust
// `sample.fastq.gz` IS in the strip list, so this one keeps the tidy stem.
```

There is no longer a list of combined suffixes to be "in". Suggested: `` // `.fastq.gz`: both suffixes are stripped, so the stem is tidy. ``

### Finding 8 — "run twice" but the table has three rows — **Low**

`tests/…:302`: "The same bytes are run twice, once under each name". The table at `tests/…:306-310` has three entries (`.fastq.bgz`, `.fq.bgz`, and the `.fastq.gz` control). Suggested: "The same bytes are run under three names, and the assertion is the same for all three".

### Finding 9 — cross-reference the other filename-derived gzip decision — **Low**

`is_gzipped` (`src/io.rs:29-31`) keeps a narrower notion of the gzip family (`extension() == "gz"`) than `strip_fastq_extensions` now does (`.gz`, `.bgz`, `.bgzf`). The divergence is deliberate and is the CHANGELOG's surviving known limitation — but the two sites do not reference each other, so whoever fixes the compression limitation has to find both. One line on `is_gzipped` would close that:

```rust
/// Note: narrower than `strip_fastq_extensions`' gzip family on purpose —
/// output compression still mirrors the name, so `.bgz` in → plain `.fq` out.
```

### Finding 10 — SE FASTQ output-collision pre-flight is missing — **Informational, separate issue**

Section 2. Pre-existing (reproduced on this branch with `.gz` names), widened by one name class here. Every other dispatch path has a pre-flight; the SE FASTQ loop at `src/main.rs:783-791` does not, so `sample.fastq.gz` + `sample.fq.gz` in one invocation silently produce a single output with two reports. Given `CLAUDE.md`'s "Output-collision pre-flight — don't bypass it", this gap is worth an issue of its own. Not a reason to hold #382.

---

## Appendix — how the equivalence table was produced

Two standalone `rustc` harnesses (scratch only, nothing written into the repo) that reimplement `dev`'s and `pr382`'s function bodies verbatim:

1. **Differential sweep.** Names generated as `stem + token¹⁻³` over stems `{sample, lib.2_R1, "", x, sample_R1}` and tokens `{.fastq, .fq, .gz, .bgz, .bgzf, .bam, .txt, .tar, .bz2, .FASTQ, .GZ, .ubam}`, each under directory prefixes `{"", /data/, ./}`, plus the degenerate paths `""`, `.`, `..`, `/`, `dir/`, `a/b/`. Result: 28,293 paths, 9,431 unique basenames, 1,544 differing basenames — and **0 differing basenames that do not end in `.bgz` or `.bgzf`**.
2. **Invariant check.** For 784 base names `X` (including `_R1` and `.L001` tokens): `new(X + ".bgz") == old(X + ".gz")`, `new(X + ".bgzf") == old(X + ".gz")`, `new(X + ".gz") == old(X + ".gz")`. 2,352 assertions, **0 violations**.

Plus the analytic proof in section 1, which does not depend on the corpus.
