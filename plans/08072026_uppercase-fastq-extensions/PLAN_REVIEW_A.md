# Plan Review A — Match FASTQ extensions case-insensitively (#384)

**Plan:** `plans/08072026_uppercase-fastq-extensions/PLAN.md` (v1, 2026-08-07)
**Target:** `dev` @ `533582a`, clean tree
**Reviewer:** A (independent; no coordination with Reviewer B)
**Method:** read `src/io.rs`, `src/main.rs`, `src/clump_only.rs`, `src/specialty.rs`,
`src/demux.rs`, `src/cli.rs`, `tests/integration_gzip_non_gz_extension.rs`, `CHANGELOG.md`;
grepped the caller sets; ran `target/release/trim_galore` (build `080c1d0`, current with `dev`)
on constructed uppercase inputs; compiled a standalone `rustc` probe for the §4 helper.

Verdict: **the decision is right and the change is small, but three of the plan's checkable
claims are wrong**, and two of them sit underneath validation items — V4 as specified cannot
fail, and §4's helper as specified panics on input that works today.

---

## 1. Logic review

### C1 — §11's caller enumeration for `is_gzipped` is wrong (Critical)

§2.3 and §11 both state `is_gzipped` "feeds only `main.rs:295`'s `gzip` flag". Two errors:

- The flag is at **`main.rs:345`**, not 295.
- There are **three more call sites**, all in `src/clump_only.rs`:
  - `:277` — `let input_compressed = naming::is_gzipped(input);` (SE `--clump_only`)
  - `:414` — `naming::is_gzipped(input_r1) || naming::is_gzipped(input_r2)` (PE `--clump_only`)
  - `:1107` — test helper `write_synthetic_fastq`, which gzip-*encodes* based on the name

`input_compressed` lands in `ClumpOnlyStats` (`:305`, `:451`) and drives two things in the
`--clump_only` text report: the `Input: … (gzip|plain, N bytes)` label (`:596`) and the gate on
the `Compression ratio: N.NNx` line (`:640`, emitted only when input *and* output are
compressed). So folding changes **report content** on the `--clump_only` FASTQ path for a
`.GZ`-named input: `plain` → `gzip`, and a ratio line appears where none did.

Direction is benign — the input genuinely *is* gzip, so both are corrections, in the same
direction as the rest of the change. Nothing silently wrong results. But:

- It is a user-visible change the plan does not list, and §7/§5 Step 6 chose `#### Changes`
  precisely because visible changes get enumerated. This one is missing.
- §2.4 "clears" `clump_only.rs` by checking **the wrong line**: it cites `:1469`
  (`assert_eq!(out.extension(), Some("fq"))`, a test assertion on a lowercase-named input) and
  concludes "Unaffected". The exposure is at `:277`/`:414`, which §2.4 never looks at.
- `clumped_output_name` (`io.rs:362`) also calls `strip_fastq_extensions`, while
  `clumping_report_name` (`io.rs:464`) uses the full input filename — so `--clump_only` carries
  the *same* stem/report mismatch #384 reports, and the same fix resolves it. Also unlisted.

§11's closing sentence — "blast radius is bounded by two functions whose callers are enumerated
above" — is unsound as written, because the enumeration is incomplete.

**Fix:** correct the line number, add the three sites, state the report-label and ratio-line
change in §7, and mention the `--clump_only` stem fix in the CHANGELOG entry.

### C2 — §3.1's last row is false; V4 as specified cannot fail (Critical)

§3.1's last row and V4 both assert that `SAMPLE.FASTQ.GZ` + `sample.fastq.gz` in one run
becomes newly rejected — "this makes the collision reachable where it was not before", and V4
"pins §3.1's last row deliberately rather than leaving it to be discovered".

**That pair is already rejected on `dev` @ `533582a`.** Verified:

```
$ trim_galore -o out_coll SAMPLE.FASTQ.GZ sample.fastq.gz     # exit 1
Error: Output path collision (case-insensitive, for APFS/NTFS safety):
out_coll/SAMPLE.FASTQ.GZ_trimming_report.txt and out_coll/sample.fastq.gz_trimming_report.txt
would be written to the same file. …
```

The reason: `report_name` / `json_report_name` (`io.rs:480`, `:496`) build their filename from
the **full input filename**, not the stem, and `collision_key` folds case. So two inputs that
differ only in case already produce case-folded-equal *report* paths and are refused before the
trimmed names are ever compared. The fold changes nothing here. V4 as written would pass
identically before and after the change — a non-discriminating test resting on a false premise.

The fold *does* create newly-reachable collisions, but only for inputs that differ in
**extension spelling as well as case**, where the report names stay distinct. Verified accepted
today:

```
$ trim_galore -o out A.FASTQ.GZ a.fq.gz       # exit 0, no collision
out/A.FASTQ_trimmed.fq   out/a_trimmed.fq
out/A.FASTQ.GZ_trimming_report.txt   out/a.fq.gz_trimming_report.txt
```

After the fold both stems become `A`/`a` and both are gzipped → `A_trimmed.fq.gz` and
`a_trimmed.fq.gz` → `collision_key`-equal → refused. **Fix:** replace §3.1's last row and V4's
inputs with `A.FASTQ.GZ` + `a.fq.gz` (or `SAMPLE.FASTQ.GZ` + `sample.fq.gz`), and drop the
"reachable where it was not before" claim for pure case-variants.

Note also that on a case-**insensitive** filesystem (default APFS, i.e. the dev machine) the two
pure case-variants cannot coexist on disk at all — passing both spellings addresses one file
twice. The `A.FASTQ.GZ`/`a.fq.gz` pair works on both filesystem kinds, which is a second reason
to prefer it.

### C3 — §4's helper, as described, panics on filenames that work today (Critical)

§4 specifies "a length check plus `eq_ignore_ascii_case` on the tail, returning the untouched
prefix". Implemented the obvious way on `&str`:

```rust
let idx = name.len() - suffix.len();
if name[idx..].eq_ignore_ascii_case(suffix) { Some(&name[..idx]) } else { None }
```

`str` slicing panics when the index is not a UTF-8 char boundary, and
`name.len() - suffix.len()` can land inside a multi-byte character. Verified with a standalone
`rustc` probe — name `😀` (4 bytes), suffix `".gz"` (3 bytes), idx 1:

```
thread 'main' panicked: start byte index 1 is not a char boundary;
it is inside '😀' (bytes 0..4) of `😀`
```

Today's `str::strip_suffix` returns `None` safely for the same input, and the binary handles
such a file end-to-end — verified: a gzipped FASTQ named `😀` (no extension) gives
`😀_trimmed.fq` + both reports, exit 0. So the naive form converts a **working** input into a
**panic**. `strip_fastq_extensions` reads `file_name().to_string_lossy()`, so any non-ASCII or
invalid-UTF-8 filename is a candidate, and extensionless inputs are explicitly supported
(`io.rs:617` pins `sample` → `sample`).

**Fix:** compare on bytes and slice only after a successful match (a matched all-ASCII tail
guarantees the index is a boundary):

```rust
fn strip_suffix_ignore_ascii_case<'a>(name: &'a str, suffix: &str) -> Option<&'a str> {
    let idx = name.len().checked_sub(suffix.len())?;
    name.as_bytes()[idx..]
        .eq_ignore_ascii_case(suffix.as_bytes())
        .then(|| &name[..idx])
}
```

`checked_sub` also removes the `name.len() < suffix.len()` underflow case. Add a non-ASCII
regression case to V3.

### I1 — V5 will be vacuous if it follows the precedent in its named home (Important)

§5 Step 5 puts the integration test in `tests/integration_gzip_non_gz_extension.rs` "where
#374/#382 put the related cases". That file already contains the exact precedent —
`bgz_output_and_report_agree_on_the_stem` (`:304`) — and it passes **`--dont_gzip`** (`:317`),
deliberately, so the expected output name does not depend on `is_gzipped`. A copy-paste of that
shape for the uppercase case would assert only the stem and never observe compression, which
makes **V6's first negative control vacuous**: reverting `is_gzipped` would not fail it.

**Fix:** V5 must run *without* `--dont_gzip` and assert three things: the trimmed file is at
`SAMPLE_trimmed.fq.gz` (name — proving the retained case is preserved end-to-end), its first
two bytes are `1f 8b`, and **both** reports exist as `SAMPLE.FASTQ.GZ_trimming_report.txt` and
`…json` (the precedent checks both; V5 says "a report").

### I2 — "share a stem" is not literally true (Important)

§5 Step 5 and V5 say the trimmed file and the report must "share a stem". They never do:
`report_name` uses the full input filename, so even in the fully-correct lowercase case the pair
is `sample_trimmed.fq.gz` + `sample.fastq.gz_trimming_report.txt`. §1's concrete expected names
are right; the wording in §5/V5 is not, and an implementer who asserts
`strip_fastq_extensions(report) == stem` will write a test that cannot pass. State the invariant
the precedent states: trimmed == `{folded_stem}_trimmed.fq[.gz]`, report ==
`{input_filename}_trimming_report.*` — i.e. the report filename *begins with* the trimmed stem.

### I3 — Specialty-mode renames are neither tabulated nor validated (Important)

All four `specialty.rs` namers take the stem (`:449`, `:471`, `:485`, `:501`), so
`--hardtrim5/3`, `--clock` and `--implicon` outputs all rename for uppercase input. §11 knows
this; §3.1 tabulates SE trim only, and V1–V7 contain no specialty assertion. Verified today:

```
$ trim_galore --hardtrim5 30 -o out SAMPLE_R1.FASTQ.GZ
out/SAMPLE_R1.FASTQ.30bp_5prime.fq        → after the fold: SAMPLE_R1.30bp_5prime.fq.gz
```

`--implicon` changes **twice over**. Its `_R1` / `R1` / `_R2`… suffix strip (`:505-517`) runs on
the stem, so for uppercase input it is currently *unreachable* (the stem ends `.FASTQ`) and
becomes reachable after the fold: `SAMPLE_R1.FASTQ.GZ` goes from
`SAMPLE_R1.FASTQ_8bp_UMI_R1.fastq` to `SAMPLE_8bp_UMI_R1.fastq.gz` — stem, double-R avoidance,
and compression all change together. Add a §3.1 row and a unit assertion on
`implicon_output_name`; `--hardtrim5` and `--clock` are in the CI validation matrix, so a
reader will want the "no uppercase fixture, therefore unaffected" reasoning (§7) attached to
them explicitly.

### I4 — §1's "one rule" is overstated: `_R1` matching stays case-sensitive (Important)

Because that suffix strip is not an extension, it keeps case sensitivity: `SAMPLE_r1.FASTQ.GZ`
under `--implicon` retains `_r1` in the output while `SAMPLE_R1.FASTQ.GZ` does not. §1 promises
"One rule, no case-dependent behaviour to document or remember"; that is not delivered. Either
fold it too or name it in §10's out-of-scope list beside `demux.rs`. Recommend out-of-scope but
**stated** — an unstated residual case-sensitivity is what produced #384.

### I5 — A3's biconditional is stronger than the argument supports (Important)

§2.3 establishes one direction soundly: folding `strip_fastq_extensions` *requires* folding
`is_gzipped`, else uppercase input yields silently plain output. The converse is not
established. Folding `is_gzipped` alone is coherent — `SAMPLE.FASTQ.GZ` →
`SAMPLE.FASTQ_trimmed.fq.gz`, where name and framing agree — it simply does not fix #384,
which is the stem mismatch. §2.3's table omits that fourth column. Stating the implication
rather than "both or neither" makes the argument airtight and pre-empts exactly that
counter-proposal.

### I6 — V2/V3's case lists miss the cases the fold actually moves (Important)

- `sample.txt.GZ`: today → `sample.txt` (`file_stem` drops `.GZ`); after the fold → `sample`
  (the `.GZ` strips, no FASTQ ext, then `file_stem` drops `.txt`). §3.1 lists `sample.txt`
  ("unchanged", true) but not the changed uppercase-gz spelling; io.rs:616 already pins the
  lowercase form, so this is the fold making the two agree — desirable, but it is a stem change
  *outside* the FASTQ list and belongs in V3.
- V2's false-list has `.bgz` and `.BGZ` but not `.BGZF`.
- V1 has `SAMPLE.FASTQ.BGZ` but not `.BGZF` uppercase.
- Add the C3 non-ASCII case.

---

## 2. Claims verified as correct

| Claim | Verdict |
|---|---|
| §2.1 reproduction | **Correct.** `-o out SAMPLE.FASTQ.GZ` → `SAMPLE.FASTQ_trimmed.fq` + `SAMPLE.FASTQ.GZ_trimming_report.{txt,json}` |
| §2.3 `is_gzipped` body and doc-comment intent | Correct (`io.rs:29-31`; doc `:12-28`) |
| §7 no uppercase fixture in `test_files/`, none in `ci.yml` | **Correct.** My regex was capable of failing — it found the one uppercase literal that exists, `io.rs:756`'s `norm_path(Path::new("FOO.FQ.GZ"))` |
| §2.4 `demux.rs` sees only paths this code names | Correct. `demux_base_name` (`:112`) takes the trimmed output path; extensions there are always ours and lowercase, including under `--basename` (which replaces the stem only) |
| §2.4 `format.rs` unaffected | Correct — content-based; no filename-extension decision in it |
| §2.4 `clump_only.rs:1469` itself unaffected | Correct **about that line** — but the file is not unaffected, see C1 |
| §3.1 `.BGZ` row (`SAMPLE` / `SAMPLE_trimmed.fq`) | **Correct.** `.BGZ` folds in the stem list; `is_gzipped` tests only `gz`, so output stays plain |
| A4's citation | Correct — the `.bgz`/output asymmetry is already documented at `io.rs:23-28` and in the CHANGELOG's Unreleased §Bug fixes |
| `strip_fastq_extensions`'s "disjoint list" comment survives folding | Correct — `.BGZ`'s 3-byte tail is `BGZ`, which does not ASCII-fold-match `.gz`, so list order still does not matter |
| §5 Step 3 compiles | Correct — `OsStr::eq_ignore_ascii_case` is stable since 1.53, below the 1.88 floor |
| `#### Changes` exists under Unreleased | Correct |
| No new filename-based decision sites beyond the two | Correct — `cli.rs` uses only `path_identity_key` (case-**sensitive**, input identity, `:537-638`), deliberately unfolded and unaffected |

On A4's sharper question — does folding *change* anything for `.BGZ` inputs? **Yes**: the stem
moves from `SAMPLE.FASTQ` to `SAMPLE`, so the trimmed file is renamed. Only compression is
unchanged. §3.1's note ("still plain — consistent with `.bgz` today") reads as "nothing
changes"; say "renamed, still plain".

---

## 3. V6's negative controls — do they discriminate?

- **"Revert `is_gzipped` alone → V5's compression assertion fails."** Discriminates *only if*
  V5 omits `--dont_gzip` (I1). Even then, the failure surfaces as the missing-file assertion
  (`SAMPLE_trimmed.fq.gz` absent, because the plain name is used) *before* the magic-byte
  check. Reword to "confirm V5 fails". Note this control is partly redundant with V2, but it
  earns its place by proving the `main.rs:345` → writer wiring is live.
- **"Revert the fold in `strip_fastq_extensions` alone → V1 fails."** Discriminates cleanly —
  `SAMPLE.FASTQ.GZ` yields `SAMPLE.FASTQ`, not `SAMPLE`.
- **Missing control:** nothing fails for an implementation that lowercases the *whole* name
  (`name.to_ascii_lowercase()` then `strip_suffix`). V1's `Sample.FastQ.Gz` row catches it and
  the plan says so — but V5 must also assert the output *filename* is exactly
  `SAMPLE_trimmed.fq.gz`, or the end-to-end path has no such guard (I1).
- **V4:** does not discriminate at all as specified (C2).

---

## 4. Efficiency

§6 is right and the section is proportionate: five suffix probes per input file, no allocation,
once per file. `checked_sub` + a byte compare is no slower than `strip_suffix`. Nothing further.

---

## 5. Alternatives

- **Fold neither, close #384 won't-fix.** Coherent but leaves a real defect: the trimmed file
  and its report disagree about the sample name, and MultiQC takes the sample name from the
  report filename — the same reasoning the #381 CHANGELOG entry and the precedent test's doc
  comment (`tests/integration_gzip_non_gz_extension.rs:298`) already give. Reject.
- **Fold `is_gzipped` only.** Fixes compression, not the reported bug. See I5.
- **Converge output compression on `format::detect_input_format` instead of the filename.** The
  honest larger fix; it would subsume A4's `.BGZ` asymmetry and retire the "known limitation"
  paragraph already recorded in the Unreleased CHANGELOG. Blast radius is every run, including
  a plain file misnamed `.gz` — which is exactly why #374/#381 deferred it. Not for this plan,
  but §10 should link the two deferrals so they stay visibly connected; folding now does not
  foreclose it.

---

## 6. Action items

### Critical

1. **C1** — Fix §2.3/§11's `is_gzipped` caller list: `main.rs:345` (not 295) plus
   `clump_only.rs:277`, `:414`, `:1107`. Add the `--clump_only` report-label and
   `Compression ratio` change to §7 and the CHANGELOG. Re-check §2.4's `clump_only.rs` line
   (it cites `:1469`; the exposure is `:277`/`:414`).
2. **C2** — §3.1's last row and V4 are false: `SAMPLE.FASTQ.GZ` + `sample.fastq.gz` is already
   rejected today via the report paths. Replace with `A.FASTQ.GZ` + `a.fq.gz` (verified
   accepted today, collides after the fold) or drop the "newly reachable" claim.
3. **C3** — Specify the helper on bytes with `checked_sub`. The `&str`-slicing form panics on a
   filename like `😀`, which works today (verified both ways). Add a non-ASCII case to V3.

### Important

4. **I1** — V5 must omit `--dont_gzip` (the precedent test at `:304` uses it and would silence
   the compression half, making V6's first control vacuous); assert the exact output name, the
   `1f 8b` magic, and both report files.
5. **I2** — Replace "share a stem" with the actual invariant; `report_name` keeps the full
   input filename.
6. **I3** — Tabulate the specialty-mode renames in §3.1 and add a unit assertion on
   `implicon_output_name`, whose `_R1` strip becomes newly reachable.
7. **I4** — Name the residual case-sensitive `_R1`/`_R2` matching in §10's out-of-scope list, or
   soften §1's "one rule" claim.
8. **I5** — State A3 as an implication, not a biconditional; add the missing fourth column to
   §2.3's table.
9. **I6** — Extend V2/V3 with `sample.txt.GZ`, `.BGZF` uppercase, and the non-ASCII case.

### Optional

10. **O1** — Drop `--dont_gzip` from §2.1's repro; the output is plain anyway, which is the
    stronger demonstration and makes §2.3 self-evident.
11. **O2** — Cite MultiQC in §7 as the downstream stake, matching the #381 entry's wording.
12. **O3** — Link the content-based-compression deferral in §10.
13. **O4** — One line noting `--basename` bypasses the stem function entirely, so it is
    correctly unaffected.
14. **O5** — No edit needed to
    `tests/integration_gzip_non_gz_extension.rs:74-82` ("still looks at the input *filename*"
    stays true after folding); say so in §5 so an implementer does not have to work it out.

---

## 7. Answer to Open 1

**Fold both, as the plan does.** §2.3's argument is sound in the direction that matters:
folding the stem without folding compression makes uppercase input produce silently
uncompressed output, which is strictly worse than the mismatch being fixed. The middle ground
the plan missed (fold `is_gzipped` only) is coherent but does not fix #384, so it is not an
alternative to this plan — and the genuinely better long-term answer (drive compression from
detected content) is out of scope for the same reason #374/#381 deferred it. Folding neither
would leave a defect that is the rewrite's own, with a concrete downstream cost in MultiQC
sample naming. I would not fold neither.

**Open 2:** `#### Changes` is right. The audience is narrow (uppercase-extension users only)
and both effects are corrections toward the documented lowercase behaviour; a release note is
not warranted, but the entry should name the `--clump_only` and specialty-mode renames (C1, I3)
so the affected surface is complete.
