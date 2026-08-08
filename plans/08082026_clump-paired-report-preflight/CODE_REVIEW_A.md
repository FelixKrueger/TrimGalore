# Code Review A — #391: clumping reports join the `--clump_only` collision pre-flights

**Reviewer:** A (independent; Reviewer B reviews the same diff separately)
**Branch:** `fix/391-clump-report-preflight` @ `2704b52` (base `dev` @ `7ea2741`)
**Plan:** `plans/08082026_clump-paired-report-preflight/PLAN.md` (r2 + Implementation notes)
**Diff:** `src/main.rs`, `src/io.rs`, `tests/integration_output_collision.rs`, `CHANGELOG.md` (+464 / −65)

## Summary

**Verdict: approve.** No Critical or High findings. The change does what the plan says, in the
place the plan says, and the load-bearing claims hold under independent re-derivation: every
candidate set keys on exactly the inputs its writer keys on, the `--no_report_file` gate is
identical in polarity and position at all five arms, and the widened namer cannot alter a written
path. `--clock`/`--implicon` are behaviourally unchanged, not merely "probably unchanged" — a
2-element `Vec` extended into `planned` produces the identical candidate list, in the identical
order, that the two `push`es produced.

Three Medium findings, all additive rather than corrective: one CHANGELOG sentence missing (a real
new over-rejection class that #388's entry documented for its own change) and two test shapes that
#388's sibling test set has and this one does not. Eight Low items, most of them one-line
trivial fixes.

### Verification performed

| Check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets -- -D warnings` | clean (debug profile; CI's `--release` was verified by the implementer) |
| `cargo test` (full) | **557 passed, 0 failed** — matches the claimed count exactly |
| `cargo test --test integration_output_collision clump` | 9/9 pass |
| Live binary: #391 repro | rejected, message pins the report path |
| Live binary: same filename in two dirs, no `-o` | exit 0, correct file layout (see M3) |
| Live binary: fold-equal filenames, `-o` | rejected (see M1/M2) |
| `tempdir` tag uniqueness across the file | all 45 tags unique |

---

## Logic — re-derived and confirmed

Each item below was checked against the tree, not taken from the plan.

**1. Five candidate sets ↔ four writer gates.** Each arm passes exactly the report-keyed inputs
its writer uses:

| Arm (`src/main.rs`) | Candidate inputs | Writer | Writer's key |
| --- | --- | --- | --- |
| paired FASTQ namer (:558) | `&[r1, r2]` | `clump_only.rs:534-536` | `input_r1` **and** `input_r2` |
| SE FASTQ (:590) | `&cli.input` | `clump_only.rs:366` | `input`, per input |
| Shape B, N=1 uBAM (:633) | `&cli.input` | `clump_only.rs:1083` | `inputs[0]` |
| Shape A, per pair (:672) | `slice::from_ref(&chunk[0])` | `clump_only.rs:1083` | `inputs[0]` |
| SE BAM (:727) | `&cli.input` | `clump_only.rs:837` | `input`, per input |

Shape B's `&cli.input` is length-1 by its branch guard (`cli.input.len() == 1`), so it is exactly
`inputs[0]` — the whole-slice spelling is safe, not an accident. Shape A's
`slice::from_ref(&chunk[0])` is the one place a `chunk[1]` slip could hide, and
`clump_paired_bam_rejects_report_that_aliases_an_input` discriminates precisely that (keyed wrong,
pair 1 plans `y.fq_clumping_report.txt`, nothing collides, run proceeds).

**2. Gate parity.** `clump_report_candidates` returns `Vec::new()` iff `no_report_file`; all four
writers skip iff `!no_report_file` is false. Every call site passes the same `cli.no_report_file`
expression it passes the corresponding writer (`main.rs:559/577`, `591/607`, `634/649`, `673/706`,
`728/744`). Position is right at all five: candidates are appended *before* the
`preflight_output_collisions` call, which itself precedes every reader open.

**3. The widened namer cannot change a written byte.** `planned` in `run_specialty_paired`
(`main.rs:2543-2547`) is consumed only by `preflight_output_collisions`; the writers re-derive
their own paths from `naming::*`. Confirmed by reading the function end to end, not by trusting the
plan's three-way verification.

**4. `--clock` / `--implicon` unchanged.** `vec![a, b]` extended == `push(a); push(b)` — same
paths, same order, so even the *first-offender* path named in an error message is unchanged.
`clock_collision_gets_the_cwd_hint_not_the_false_advice` (`:642`) passes.

**5. `run_specialty_paired`'s doc-comment mode list is now complete.** Exactly three callers exist
(`:483` clock, `:498` implicon, `:547` clump paired); `--hardtrim5/3` use their own loops, so the
list "(`--clock`, `--implicon`, `--clump_only --paired`)" is exhaustive rather than merely longer.

**6. No fifth writer.** `grep clumping_report_name src/` returns only the `io.rs` definition, the
four `clump_only.rs` writer sites, two `clump_only.rs` tests, the two `io.rs` property tests, and
`main.rs:128` (the new helper). The trim-path `--clumpify` writes no clumping report.

**7. Hint constant has no string-continuation whitespace defect.** Every continuation line in
`main.rs` ends with a space before the `\` (verified by regex `[^ ]\\$` over the whole file —
zero hits). The rendered text was confirmed by running the repro:

> `Output path collision (case-insensitive, for APFS/NTFS safety): …/out/reads.fq_clumping_report.txt and …/out/reads.fq_clumping_report.txt would be written to the same file. Outputs and reports are named from the input filename alone, so inputs sharing a filename can collide when --output_dir (or a shared input directory) sends them to one place — rename one input, or pass --no_report_file if only the reports collide.`

The wording matches the plan's approved text verbatim (backticks dropped, correct for stderr). The
old doc comment on the constant was actually *wrong* — it described `CWD_OUTPUT_HINT`'s role — so
the rewrite fixes a stale comment as well as stale prose. Both pre-existing use sites (`:838` trim
paired FASTQ, `:1942` trim→uBAM) were re-read: their surrounding comments reference `#388` and the
`_val_` discriminator, not the hint's wording, so nothing went stale there.

**8. No new false-rejection class from report-vs-primary.** Report names always end
`_clumping_report.txt`; clump primaries end `_clumped.fq`, `_clumped.fq.gz`, or `_clumped.bam`. The
two families can never be equal, so the only new candidate-pair dimension is report-vs-report and
report-vs-input, both of which are genuine overwrites. The `--no_report_file` acceptance siblings
plus `clump_paired_accepts_distinct_filenames`'s exact-set assertion pin this from the other side.

**9. Property-test direction is sound.** The new 3×3 cross-product (primaries
`single_end_output_name` / `clumped_output_name` / `clumped_bam_output_name` vs secondaries
`report_name` / `json_report_name` / `clumping_report_name`) asserts nine instances of one true
implication: all three primaries key on the *stem*, all three secondaries on the *full filename*,
and both resolve directory identically (`output_dir` else the input's own parent). Distinct stems
imply distinct filenames, so "primaries differ ⇒ secondaries differ" holds in every cell.

---

## Efficiency

Nothing beyond the plan's O(paths) claim. Candidate growth is 2P → ≤4P on the paired FASTQ arm and
+1 per input/pair elsewhere; `preflight_output_collisions` remains two hash lookups per candidate
with `HashMap::with_capacity(planned.len())`. Allocation delta is one `Vec` per helper call plus one
per pair on the paired arm — on a code path about to read gigabytes. `slice::from_ref` in the Shape
A loop avoids a per-chunk temporary array. Making the helper return an iterator would remove one
`Vec` per call at the cost of complicating the gate; not worth it. **No findings.**

---

## Errors

No bugs found. Two non-issues checked and closed:

- `clump_report_candidates` returns `Vec` rather than `Result` — correct, since unlike
  `planned_secondary_outputs` (which reads the demux barcode file) it performs no I/O.
- `clumping_report_name`'s `file_name().unwrap_or_default()` yields an empty name for a path
  ending in `..`; pre-existing, and such a path fails as an unreadable input earlier.

The closure at `main.rs:553` newly captures `&cli` (for `cli.no_report_file`) alongside
`run_specialty_paired(&cli, …)`'s own immutable borrow. Both immutable, and the `run_pair` closure
already captured `cli`, so no new borrow surface. Compiles and passes.

---

## Findings and recommendations

### M1 — CHANGELOG omits the new case-only over-rejection (Medium)

Verified live on this build:

```
--clump_only --paired -o out a/Reads.fq b/reads.fq
→ Error: … out/Reads.fq_clumping_report.txt and out/reads.fq_clumping_report.txt
  would be written to the same file.
```

Both primaries (`Reads_clumped_1.fq` / `reads_clumped_2.fq`) differ even folded, so this run
succeeded before the patch and on a case-sensitive filesystem both reports genuinely could
coexist. It is the same loud-error-over-silent-loss trade #216 established — and #388's CHANGELOG
entry spells it out explicitly for its own change ("two genuinely distinct inputs whose filenames
differ only in case are refused even on a case-sensitive filesystem, where both reports could in
fact coexist"). #391's entry does not mention it, so a user hitting the refusal on ext4 has no
release note to point at.

**Recommendation:** add one sentence to the #391 bullet mirroring #388's, e.g. *"As with #388, two
inputs whose filenames differ only in case are now refused into a shared output directory even on a
case-sensitive filesystem, where both reports could coexist."*

### M2 — no fold-equal integration test on the clump arms (Medium)

`paired_rejects_fold_equal_filenames_into_shared_output_dir`
(`tests/integration_output_collision.rs:886`) is #388's dedicated fold-dimension test — six lines.
The clump set has no counterpart: `clump_paired_rejects_cross_pair_report_collision` uses
*exactly*-equal filenames (`a/x.fq a/y.fq b/y.fq b/x.fq`). The plan's Behavior 2b promised
"cross-pair **fold-equal** report names under `-o`" but the test it specced in step 8 is the
exact-equal shape, so the gap is inherited from the plan rather than introduced here. The
`io.rs` property test covers the fold dimension at the *namer* level; nothing exercises an actual
fold-equal clump-report rejection end to end.

**Recommendation:** add the six-line twin of `:886` (rejection asserted, so filesystem-independent);
M1's message shows both distinct filenames, so it can pin both. This is also the only clump test
that would exercise the `collision_key` fold path through the whole binary.

### M3 — no over-rejection guard for `output_dir = None` (Medium)

#388 has `paired_report_candidates_do_not_over_reject` (`:970`), whose third case exists precisely
to pin that report candidates honour `output_dir = None`. The clump set has no equivalent, and the
consequence is concrete: a helper that resolved report directories into `-o`-or-cwd instead of each
input's own parent would pass **all nine** new tests. Walking them —
`clump_paired_rejects_shared_mate_report_without_output_dir` and
`clump_se_rejects_report_that_aliases_an_input` both still reject under that mutation (paths
normalise to the same keys), the two `--no_report_file` cases plan no reports at all, and the rest
pass `-o`.

Verified the missing shape currently behaves correctly:

```
--clump_only --paired A/reads.fq B/reads.fq        (no -o)  → exit 0
A/: reads.fq  reads_clumped_1.fq  reads_clumped_2.fq  reads.fq_clumping_report.txt
B/: reads.fq  reads.fq_clumping_report.txt
```

Mitigating: `clumping_report_name`'s `None` rule is unit-pinned at `io.rs:958-959`, and the helper
delegates to it rather than re-deriving, so the mutation is not reachable without editing `io.rs`.
That is what keeps this Medium rather than High.

**Recommendation:** add `clump_report_candidates_do_not_over_reject` in the shape of `:970` —
same filename in two directories, no `-o`, assert exit 0 and both per-input reports beside their
own inputs. It closes the file's own stated convention ("every rejection case … paired with an
acceptance case on the same dispatch path") and doubles as executable documentation of the
primaries-follow-R1 / reports-follow-each-mate asymmetry the plan lists as a follow-up.

### L1 — `assert_dir_holds_only` no longer says which directory failed (Low, trivial fix)

The generalised message (`tests/…:86-89`) dropped "a rejected run must write nothing" without
gaining the directory name, and this diff is the first caller to invoke it three times in one test
(`clump_paired_rejects_shared_mate_report_without_output_dir` checks `p/`, `d/`, `q/`). On failure
you get two unlabelled file lists.

**Trivial fix:** `"{} must hold exactly the expected files", dir.display()`.

### L2 — three-line inline comment against the project's two-line max (Low, trivial fix)

`main.rs:630-632`. Suggested two-line form: `// #391 — defensive symmetry: with one input the
report can never alias it; / // the arm keeps its siblings' shape.` The other new inline comments
(`:550-551`, `:589`, `:670-671`, `:726`) are within budget. The 4-line `///` doc on
`clump_report_candidates` matches `planned_secondary_outputs`' 8-line doc directly above it, so it
reads as house style rather than excess.

### L3 — `io.rs:1239-1241` undersells and slightly misdescribes its own fixture (Low, trivial fix)

"primaries fold-collide (skipped pair)" reads as though the new `d/SAME.fastq.gz` input proves
nothing. Its actual role is stronger and worth recording: it *couples the skip predicate to the
assertion metric*. Under the old `a == b` skip, `d/SAME` vs `d/same` with `-o` would not be skipped
(the two `PathBuf`s differ), and the new `collision_key`-based secondary assertion would then
**fail**. The fixture is what makes the `collision_key` skip load-bearing rather than cosmetic.

**Trivial fix:** `// Fold-equal twin of d/same: the skip must use collision_key, not PathBuf
equality (#391).`

### L4 — `io.rs:1311-1313` doc invites a misreading (Low)

The doc now says "`_clumped_N` inverts it the same way, #391" three lines above an assertion that
`clumped_paired_bam_output_name` does **not** invert (it collapses three spellings onto one
primary). Both statements are true of different namers: the inverting one is
`clumped_paired_output_names` (FASTQ, `_clumped_1`/`_clumped_2`), which is deliberately covered by
*neither* property test because the property is false there — that is the whole of #391. Worth one
clause naming which namer inverts, so the next reader does not think the added `clumped_pe_bam`
block contradicts the doc.

### L5 — the helper's "unit-testable" justification is unrealised (Low, informational)

`clump_report_candidates` takes a bare leading `bool` (rather than `&Cli`, as
`planned_secondary_outputs` does) partly so it can be unit-tested in isolation. `src/main.rs` has
no `#[cfg(test)]` module and is the bin crate, so the helper is unreachable from both unit and
integration tests; the signature choice buys nothing today. Behaviour is fully covered by the nine
integration tests, and the shape matches its neighbour, so no change is needed — flagged only so
the claim is not relied on later. Moving it into the library (beside `naming`) would make it
testable and put it next to the namer it wraps; that is a refactor, not a fix.

### L6 — two idioms now coexist for one job (Low)

`main.rs` builds clump report candidates through a helper but trim report candidates inline at
`:825-832` and `:1931-1937` (both from #388). A `trim_report_candidates` twin would unify them and
remove the second inline copy. Follow-up material — doing it here would widen a bug-fix diff into
the trim paths.

### L7 — `docs/…/modes/clump-only.md` is wrong in exactly this area (Low, separate issue)

Both pre-existing, neither introduced here, and #388 set the precedent of not touching docs — but
they misdescribe the mechanism #391 is about:

- `:109` calls the report `<stem>_clumping_report.txt`. It is `<input-name>_clumping_report.txt`,
  as `:45` of the same page correctly says. Stem-vs-full-name is the exact confusion that produced
  #391.
- `:122` says "Multi-pair PE runs produce one report per pair" unscoped. True on the uBAM arm
  (`clump_only.rs:1083`), false on the paired FASTQ arm, which writes one report per *mate*
  (`clump_only.rs:534-536`) — as `clump_paired_accepts_distinct_filenames` now asserts.

### L8 — CHANGELOG placement (Low, editorial)

#391's bullet is first in `#### Bug fixes`, #388's is last, 116 lines apart, though they are one
defect family and #391's prose recapitulates #388's mechanism. The section is not in merge order
(#388 merged most recently), so there is no rule being broken — but the release notes would read
better with the two adjacent. Felix's call.

---

## Pre-existing residual, re-confirmed (no action)

The report-vs-report message prints the same path twice ("`…/out/reads.fq_clumping_report.txt` and
`…/out/reads.fq_clumping_report.txt`"), because `preflight_output_collisions` has no input
provenance. Confirmed live. Already recorded as a plan follow-up (B-A3); the hint's
`--no_report_file` clause is what makes the message actionable in the meantime.

---

## Test-convention audit

Conforms to the file's established conventions:

- Unique `tempdir` tags — all 45 in the file are distinct, and `remove_dir_all` targets exact
  paths, so the prefix-sharing pairs (`clump_rep` / `clump_rep_norep`, `clump_se_alias` /
  `clump_se_alias_ok`) cannot wipe each other.
- Rejection/acceptance pairing on the same dispatch path: paired FASTQ (4 rejections + 2
  acceptances), SE FASTQ (1 + 1), Shape A uBAM (1 + 0 — acceptance covered by the pre-existing
  multi-pair BAM tests). The two gaps are M2 and M3.
- Attributable read prefixes throughout (`CR1`, `CN1`, `D1`, `X1`, `MX`, `B1`, `SRC`, `REPORTY`,
  `PX`, `PR`), and acceptance cases assert content via `count_reads_from` rather than mere
  existence.
- Test names mirror the #388 family (`clump_paired_rejects_cross_pair_report_collision` beside
  `paired_rejects_cross_pair_report_collision`), which makes the two sets readable together.
- Doc comments explain *why each test would fail without the fix*, matching the module's tone.
  `clump_paired_bam_rejects_report_that_aliases_an_input`'s comment naming the exact surviving
  mutation is the strongest of the nine.

`assert_dir_holds_only`'s repurposing for acceptance-side exact-set assertions (documented as a
deviation) is the right call — it is the assertion that pins "candidate list == written set", the
invariant whose absence caused both #388 and #391.
