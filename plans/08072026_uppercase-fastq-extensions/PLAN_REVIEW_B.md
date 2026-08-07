# Plan Review B — Match FASTQ extensions case-insensitively (#384)

**Plan:** `plans/08072026_uppercase-fastq-extensions/PLAN.md` (v1, 2026-08-07)
**Target:** `dev` @ `533582a`, clean tree
**Reviewer:** B (independent; no coordination with Reviewer A)

**Verdict:** the approach is right and the plan is unusually careful, but its blast-radius
statement is wrong in one place (`is_gzipped` has four call sites, not one), and the parity
argument that licences folding `is_gzipped` does not survive reading the Perl source.
1 Critical, 6 Important, 6 Optional.

---

## What I verified (and that the checks could fail)

| Claim | Result |
|---|---|
| §2.1 reproduction | **Confirmed**, and stronger than stated — see I5 |
| `is_gzipped` feeds only `main.rs:295` (§2.3, §11) | **False.** 4 call sites; line is 345 → **C1** |
| `strip_fastq_extensions` callers (§11) | 12 in `io.rs`, 4 in `specialty.rs`; two namers bypass it → O2 |
| §2.4: `demux.rs` sees only self-named paths | **Confirmed** (`main.rs:97-99` → `demux_base_name`, `demux.rs:111`) |
| §2.4: `format.rs`, `clump_only.rs:1469` unaffected | **Confirmed** |
| §7: no uppercase extension in `test_files/` or `ci.yml` | **Confirmed**, and extended to `src/`+`tests/`+`docs/`: one hit, `io.rs:756` (`norm_path("FOO.FQ.GZ")`), which folding cannot affect. The grep found a hit, so it was capable of finding others. |
| §3.1: `SAMPLE.FASTQ.GZ` + `sample.fastq.gz` rejected after fold | **True by construction** (`collision_key` = lowercased absolute path; SE `planned` covers all inputs, `main.rs:825-834`) — but the scenario is unreachable on APFS → **I1** |
| §3.1 `.BGZ` row | **Accurate**; stem changes, compression does not → O1 |
| §2.2 Perl naming (`trim_galore` 909-923, 611) | **Confirmed verbatim** from `git show 0.6.11:trim_galore` — and line 925 adds something the plan missed → **I2** |

Reproduction, run in the scratchpad with `target/release/trim_galore` (built from `080c1d0`,
the pre-merge commit of #383; `git diff 080c1d0 HEAD -- src/io.rs` touches only
`collision_key`/`path_identity_key` and tests, leaving all three functions under review
byte-identical, so it is valid for these observations):

```
SAMPLE.FASTQ.GZ  → SAMPLE.FASTQ_trimmed.fq (ASCII text)  + SAMPLE.FASTQ.GZ_trimming_report.{txt,json}
SAMPLE.FASTQ.BGZ → SAMPLE.FASTQ_trimmed.fq (ASCII text)  + SAMPLE.FASTQ.BGZ_trimming_report.{txt,json}
ref.fastq.gz     → ref_trimmed.fq.gz                     + ref.fastq.gz_trimming_report.{txt,json}
```

---

## Critical

### C1 — `is_gzipped` has four call sites, and `--clump_only`'s report is one of them

§2.3 and §11 both assert `is_gzipped` "feeds only `main.rs:295`'s `gzip` flag". The call sites are:

- `src/main.rs:345` — the `gzip` flag (the line number in the plan is also stale)
- `src/clump_only.rs:277` — `input_compressed`, single-end `--clump_only`
- `src/clump_only.rs:414` — `input_compressed`, paired `--clump_only`
- `src/clump_only.rs:1107` — a test helper (`write_synthetic_fastq`), harmless

`input_compressed` is not cosmetic. In `write_clump_only_report` it selects the report's
input label (`clump_only.rs:596`: `"gzip"` vs `"plain"`) and gates the compression-ratio
line entirely (`clump_only.rs:640`: `input_compressed && output_compressed && output_bytes > 0`).
`gzip_output` on that path is the same folded `gzip` flag (`main.rs:526/534`, `550/558`). So for
a `.GZ`-named input under `--clump_only`, folding changes **three** things at once:

```
Input:  SAMPLE.FASTQ.GZ (plain, N bytes)   →   Input:  SAMPLE.FASTQ.GZ (gzip, N bytes)
Output: SAMPLE.FASTQ_clumped.fq (plain…)   →   Output: SAMPLE_clumped.fq.gz (gzip level 1…)
(no ratio line)                            →   Compression ratio: X.XXx
```

None of this is wrong — the label becomes *more* truthful, since a `.GZ` file really is gzip —
but it is an output-content change to a second mode that the plan does not enumerate, does not
put in the CHANGELOG (Step 6), and does not validate. Nothing catches it if it regresses: §7
confirmed no test anywhere uses an uppercase extension, so the whole clump path is dark here.

**Action:** correct §2.3/§11's enumeration; add the `--clump_only` report change to §3.1, §7
and Step 6; add a validation (see V8 below).

---

## Important

### I2 — the compression half *does* have a parity defence, contra §2.2

§2.2's load-bearing sentence is "parity cannot arbitrate this decision because parity is already
lost". Verified against `git show 0.6.11:trim_galore`: for the **name**, that is exactly right —
line 923's `else { s/$/_trimmed.fq/ }` appends to the full filename and line 611's report does the
same, so Perl agreed with itself and we do not.

But eighteen lines later, Perl decides compression:

```perl
if ($gzip or $filename =~ /\.gz$/){
```

No `/i`. Perl wrote `SAMPLE.FASTQ.GZ_trimmed.fq` **plain** — which is what `dev` does today.
So today's Rust is *parity-correct on compression* for `.GZ` input, and folding `is_gzipped`
breaks a parity that currently holds, while folding `strip_fastq_extensions` restores nothing
Perl had either. §2.2's "parity is already lost" is true of the naming half only; the plan
generalises it to license the compression half, which is the one place the argument is unearned.

This does not change my recommendation (see Open 1), but the plan should say so plainly rather
than leaning on a parity claim that inverts between the two functions — the CI validation matrix
is unaffected only because no fixture is uppercase, not because parity is irrelevant here.

### I1 — the new refusal §3.1 picks cannot happen on APFS; the one that can is missing

§3.1's last row and V4 use `SAMPLE.FASTQ.GZ` + `sample.fastq.gz`. Those two names are **one file**
on APFS — verified: the second `cp` prompted to overwrite the first. The row is therefore a
Linux/case-sensitive-FS-only scenario, correctly scoped as a *unit* test (the plan gets that
right) but it must never be promoted to an integration test, and it should be labelled as what
it is: a case-folding **false positive** on case-sensitive filesystems — two genuinely distinct
outputs refused. That is the trade-off `norm_path`'s own doc-comment (`io.rs:39-42`) already
accepts ("a loud early error rather than silent data loss"), so it is defensible; "the pre-flight
working as designed" undersells it.

Meanwhile the collision users will actually hit is absent from the plan. Verified on `dev` today:

```
$ trim_galore -o out_c1 SAMPLE.FASTQ.GZ SAMPLE.FQ.GZ     # succeeds, 6 output files
$ trim_galore -o out_c2 s2.fastq.gz     s2.fq.gz         # Error: Output path collision …
```

After the fold both uppercase inputs stem to `SAMPLE`, so the first command becomes a hard
error — on **every** filesystem, no case-sensitivity required, and not a false positive (the two
outputs really are one path). Same for `.FASTQ.GZ` + `.FASTQ.BGZ`. This is the row that belongs
in §3.1 and in the CHANGELOG, and it is integration-testable anywhere.

### I3 — V5's "share a stem" is either self-referential or false as written

`report_name` (`io.rs:480-492`) formats `{full input filename}_trimming_report.txt`. So the
trimmed file and the report never *literally* share a stem, in any case — verified above,
lowercase included (`ref_trimmed.fq.gz` vs `ref.fastq.gz_trimming_report.txt`). "Share a stem"
only holds under a stripping rule, and the only rule under which it holds is
`strip_fastq_extensions` itself — the function under test. A V5 implemented as
`strip_fastq_extensions(report_input_part) == trimmed_stem` passes for *any* self-consistent
fold, including a wrong one (e.g. one that lowercases the retained part: `sample` == `sample`).

**Action:** V5 asserts exact literal filenames — `SAMPLE_trimmed.fq.gz` exists,
`SAMPLE.FASTQ.GZ_trimming_report.txt` exists, and the first two bytes of the trimmed file are
`1F 8B` — with no call to any function under test in computing the expectations.

### I4 — the change does not deliver the MultiQC grouping §2.1 invokes

§2.1 and the #381 precedent frame the mismatch as "disagree about the sample name (which is what
MultiQC groups on)". After the fold, the report filename still carries `.FASTQ.GZ`, and
`report.rs` writes the input filename verbatim into the report body (`input_filename`,
`input_filenames`). A consumer stripping extensions case-sensitively — which is how MultiQC's
`fn_clean_exts` list works, worth confirming before writing this in the CHANGELOG — still
derives `SAMPLE.FASTQ.GZ` from the report and `SAMPLE` from the trimmed file. Two sample names,
as before.

The §1 goal (one rule, uppercase behaves as lowercase) *is* delivered, and that is a good enough
reason. But §2.1/§7 should not imply downstream grouping is fixed, because it is not — only the
lowercase path was ever groupable, and the fold makes uppercase internally consistent, not
downstream-clean.

### I5 — the §2.1 reproduction hides the compression symptom, so V7's re-run can't see it

§2.1 (and V7's "re-run the §2.1 reproduction") passes `--dont_gzip`, which forces plain output
and makes the compression half invisible. Verified without it: current output is still
`SAMPLE.FASTQ_trimmed.fq`, plain ASCII text — because `is_gzipped` returned false. Drop
`--dont_gzip` and the one command demonstrates both halves before and after. As written, V7's
final gate cannot observe the change C1/Open 1 are about.

### I6 — V6's negative controls: one discriminates, one is redundant and mis-aimed

Worked through both as specified:

- **"revert the fold in `strip_fastq_extensions`, confirm V1 fails"** — genuinely
  discriminating. `SAMPLE.FASTQ.GZ` → `SAMPLE.FASTQ` ≠ `SAMPLE`; V1 fails on its first case. ✅
- **"revert `is_gzipped`, confirm V5's compression assertion fails"** — the test does fail, but
  not where claimed. With the stem folded and `is_gzipped` reverted, the output is
  `SAMPLE_trimmed.fq`; V5 fails at *path resolution*, never reaching the magic-byte assertion.
  It discriminates, but not on the compression axis, and the failure message will point at a
  missing file. Fix by having V5 resolve whichever of the two candidate paths exists, then assert
  on the magic bytes separately — so the two halves fail independently and legibly.

Two further points. V6 is a **manual one-off**, not a committed guard; the permanent guard against
re-losing the `is_gzipped` half is V2, which already pins `is_gzipped(".GZ") == true`. Say that,
so V6 is not mistaken for a regression test. And no control covers C1's `--clump_only` report
path, which is currently the least-guarded consequence of the change.

**V8 (new):** `--clump_only` on a `.GZ`-named input — assert the report's `Input:` line says
`gzip`, that a `Compression ratio:` line is present, and that the output is
`<stem>_clumped.fq.gz`. This is the only proposed test that would catch a regression in C1's
newly-identified consumer.

---

## Optional

- **O1.** §3.1's `.BGZ` row is accurate on the Output column (verified: plain before *and* after)
  but is not marked as changed, though its stem moves `SAMPLE.FASTQ` → `SAMPLE`. Only the `.GZ`
  row carries a **changed** marker. One word.
- **O2.** §11's "every output namer" overstates: `main.rs:1646` and `main.rs:2219` (interleaved-uBAM
  paired FASTQ and BAM output) name from a raw `input.file_stem()` and bypass
  `strip_fastq_extensions` entirely. Both are BAM-input-only, so nothing to fold, but the
  enumeration should carry the exception so a later reader does not read it as a bug.
- **O3.** `--basename` substitutes for the stem (`io.rs:169`, `195`, `222`, …), bypassing
  `strip_fastq_extensions`, so `--basename SAMPLE.FASTQ` is unaffected by the fold. Consistent
  with A2; one line in §2.4 stops it being reported later as a residual inconsistency.
- **O4.** V1/V3 have no uppercase *non*-FASTQ cases. Two are worth pinning: `SAMPLE.TXT.GZ`
  (stem `SAMPLE.TXT` → `SAMPLE`, mirroring the existing lowercase assertion at `io.rs:616`) and
  bare `SAMPLE.GZ` (stem stays `SAMPLE`, but compression flips plain → gzip). The second is the
  only input class where compression changes while the name does not — exactly the case a reader
  of the CHANGELOG will not predict.
- **O5.** The disjointness comment at `io.rs:530-531` ("`.bgz` does not end in `.gz`") still holds
  under ASCII-insensitive matching, and `.bgzf` ends in neither, so `find_map` order stays
  irrelevant. Worth stating in the commit message rather than leaving the next reader to re-derive it.
- **O6.** Implementation nit for Step 2: the current chain is
  `find_map(|ext| name.strip_suffix(ext)).unwrap_or(&name)`, so the helper must return a
  `&'a str` borrowed from `name` (as §4's signature does) and the second stage must borrow from
  that same chain — not from a `String` temporary — or it will not borrow-check. §4 has this
  right; flagging only because it is the one place the implementation can trip.

---

## Assumptions

- **A1 (ASCII-only)** — sound, and consistent with `norm_path` (`io.rs:44`), which the collision
  pre-flight already relies on. Non-ASCII would additionally desynchronise the fold from
  `collision_key`, so ASCII is not just convenient, it is required for §3.1 to hold.
- **A2 (only matching folds)** — verified: every namer composes `format!` with lowercase literals.
  Nothing emits uppercase. Add O3 for completeness.
- **A3 (both or neither)** — sound in the direction the plan argues (see Open 1) but stated too
  strongly; there is a third option it does not name, and I2 weakens the parity ground it stands on.
- **A4 (`.BGZ` folds the stem but not compression)** — accurate. Verified that folding changes
  nothing about `.BGZ` compression: `extension()` is `BGZ`, and `eq_ignore_ascii_case("gz")` is
  false for it before and after. The asymmetry is genuinely pre-existing (`is_gzipped`'s own
  doc-comment, `io.rs:24-28`, already confesses it for lowercase `.bgz`).
- **A5 (renames accepted)** — the affected population is wider than "those currently receiving a
  mismatched pair": it also includes `--clump_only` users (C1) and multi-input runs that will now
  be refused (I1).

**Unstated assumption worth surfacing:** that `gzip` is derived from `cli.input[0]` alone
(`main.rs:345`). A run mixing `sample.fastq.gz` and `SAMPLE2.FASTQ.GZ` inherits the *first*
file's compression for both, before and after the fold. Pre-existing and documented in the
comment at `main.rs:338-344`, but the fold changes *which* runs are affected, so it is worth one
line rather than a surprise.

---

## Efficiency

Agreed and uncontroversial: `eq_ignore_ascii_case` over a ≤5-byte suffix, five probes, once per
input file, no allocation. Same order as `strip_suffix`'s memcmp. Nothing further.

---

## Alternatives

**The middle ground §2.3 misses.** `is_gzipped`'s doc-comment (`io.rs:24-28`) already names the
real fix and defers it: source the output-compression decision from `format.rs`'s *detected*
format instead of the filename. That resolves `.GZ`, `.BGZ`, `.BGZF` and a mis-named plain `.gz`
in one move, and it dissolves A3 — with compression sourced from content, folding
`strip_fastq_extensions` alone is safe and complete. The reason it is deferred is real (a plain
file misnamed `.gz` currently gets gzipped output and would stop), which is why I would not do it
under #384. But it should appear in §10 as the destination, so this fold reads as an interim step
rather than the end state.

**Restore Perl's append-fallback** — §2.2 rejects it on uBAM grounds and that reasoning holds
(`sample.bam_trimmed.bam` is indefensible). Agreed, no need to revisit.

**Fold neither, close won't-fix** — leaves a rule with no defence: the naming half is not Perl's
either (verified), so "we match Perl" would be false as stated. Rejecting this is right.

---

## Answers to the four questions put to me

1. **Open 1 — is A3 sound?** Directionally yes: folding the stem alone yields `SAMPLE_trimmed.fq`,
   a name that *looks* like the gzipped convention's sibling but is silently plain, which is worse
   than today's honest mismatch. But "both or neither" is false — the content-sourced compression
   decision above is a coherent third option, and I2 shows the compression half is the one place
   parity currently favours *not* folding. **Fold both**, as the plan proposes, and say in §10
   that the format-sourced decision is the follow-up that makes A3 unnecessary.
2. **A4 / `.BGZ`** — folding changes the `.BGZ` stem (`SAMPLE.FASTQ` → `SAMPLE`, verified) and
   changes nothing about its compression. §3.1's `.BGZ` row is what the code will do; it is only
   mislabelled as unchanged (O1). Refusing to widen the asymmetry is right.
3. **V6's controls** — one discriminates cleanly, one discriminates for the wrong reason and is
   redundant with V2; neither covers C1. See I6, including V8.
4. **Anything else reading the input filename** — §2.4's list is incomplete: `clump_only.rs:277`
   and `:414` (C1) are the material omissions. `main.rs:1646`/`:2219` and `--basename` read names
   but cannot disagree with the fold (O2, O3). `format.rs`, `fastq.rs`, `cli.rs`, `fastqc.rs`,
   `report.rs` and `demux.rs` are clean — swept with
   `grep -rn "extension()\|file_stem()\|strip_suffix\|ends_with(\"" src/*.rs`.

**Open 2** — `#### Changes` is right; no separate release note. But the entry needs three bullets,
not two: the rename, the compression flip, and the newly-refused multi-input runs (I1). Plus the
`--clump_only` report shift if C1 is accepted.

---

## Action items

**Critical**
1. Correct §2.3/§11: `is_gzipped` has four call sites (`main.rs:345`, `clump_only.rs:277`, `:414`,
   `:1107`), not one. Document the `--clump_only` report consequence in §3.1/§7/Step 6 and add
   V8 to cover it. (C1)

**Important**
2. Fix §2.2/§2.3's parity argument: Perl's `if ($gzip or $filename =~ /\.gz$/)` is case-sensitive,
   so today's compression behaviour *is* Perl-identical and folding departs from it. (I2)
3. Add the `SAMPLE.FASTQ.GZ` + `SAMPLE.FQ.GZ` collision to §3.1 and the CHANGELOG; relabel the
   existing last row as a case-sensitive-FS-only false positive and keep V4 a unit test. (I1)
4. Restate V5 as exact literal filenames plus a magic-byte check, computing no expectation from a
   function under test. (I3)
5. Drop the MultiQC-grouping implication, or state that it is not fixed for uppercase input. (I4)
6. Drop `--dont_gzip` from §2.1 / V7's reproduction. (I5)
7. Re-scope V6: split V5's path and compression assertions so each half-revert fails legibly; note
   V2 is the committed guard. (I6)

**Optional**
8. O1 — mark §3.1's `.BGZ` row as changed. 9. O2 — note the two uBAM namers that bypass the strip.
10. O3 — note `--basename` bypasses it too. 11. O4 — add `SAMPLE.TXT.GZ` and bare `SAMPLE.GZ` to
V1/V3. 12. O5 — record the disjointness argument in the commit message. 13. O6 — lifetime note for
Step 2.
