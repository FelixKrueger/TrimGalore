# Plan Review A — run the `--phred64` pre-flight whenever output is uBAM

**Plan:** `plans/11062026_ubam-phred64-preflight/PLAN.md`
**Reviewer:** A (independent; fresh context)
**Reviewed at:** `dev` @ `49152c1`, release binary `49152c1 — macos/aarch64 — built 2026-11-06T20:29:12Z`
**Verdict:** **Sound. No Critical findings.** The diagnosis is correct, the measurements reproduce
exactly, the placement is right, the mode set is complete, and the #362 regression it is most at risk
of causing does not occur. Three Important items concern the *plan's own reasoning and validation
plan* rather than the change: one prescribed mutation cannot fail the test the plan assigns it, the
§5.3 NOTE correction is declared required but has no validation entry, and §11 mislabels a
load-bearing guard as decorative.

---

## 0. What I verified independently

I re-derived or re-ran everything load-bearing rather than taking the plan's word for it. Scratch:
`/tmp/claude-501/revA-phred64-1794048757/`.

### 0.1 The §2 measurement table — reproduces

Fresh fixture, all-`I` (ASCII 73 = Q40) Phred+33, 4 × 20 bp reads. `samtools view` renders `QUAL`
as ASCII+33, so `*` (42) = Phred 9.

| Plan's claim | My result | Match |
|---|---|---|
| `--output-format ubam --phred64` (default `-q 20`) → exit 0, empty BAM | exit 0, `nrec=0`, reports written | ✅ |
| `… -q 5` → Q40 → Q9 | exit 0, `QUAL` = `***…` = Phred 9 | ✅ |
| `… -q 0` → Q40 → Q9 | exit 0, Phred 9 | ✅ |
| `--hardtrim5 20 --output-format ubam --phred64` → Q9 | exit 0, `p33_I.20bp_5prime.bam`, Phred 9 | ✅ |
| `--clump_only --output-format ubam --phred64` → Q9 | exit 0, `p33_I_clumped.bam`, Phred 9 | ✅ |
| `--hardtrim5 20 --phred64` (FASTQ out) → verbatim | `IIIIIIIIIIIIIIIIIIII` | ✅ |

One immaterial discrepancy: the plan records "2 records written" for the `-q 5` row; I got 4, because
my fixture has 4 reads. Fixture-size artefact, not an error in the finding.

### 0.2 The mixed-quality flattening — exact match

Fixture spanning Q0–Q40 through `--output-format ubam --phred64 -q 0 --length 0`:

```
in : 0  5 10 15 20 25 30 31 32 35 38 40
out: 0  0  0  0  0  0  0  0  1  4  7  9
```

Identical to the plan's table, and consistent with `max(0, q − 31)`. Eight distinct input levels
collapse onto Q0, i.e. seven levels destroyed — the plan's "7 of 12" is right. The "48 %" figure is
base-composition dependent (I measured 67 % on an evenly weighted fixture); it is a property of the
plan author's fixture, not a claim about the code, so it does not need correcting.

`saturating_sub` is confirmed load-bearing at `src/bam.rs:615` ff., with its own comment saying so.
The plan's §2 point 2 — lossy floor, not a reversible offset, so #358's recoverability does not carry
over — is correct.

### 0.3 The startup NOTE really is false on the run that corrupts

Captured verbatim from the stderr of `--hardtrim5 20 --output-format ubam --phred64`, the run that
then wrote Phred 9:

> `… Modes that do no quality arithmetic, such as --hardtrim5/3 and --clump_only writing FASTQ, pass the scores through untouched. …`

§5.3 is justified. (The sentence is at best ambiguous — "writing FASTQ" could be read as a trailing
modifier on both items — but it is printed unconditionally, including on uBAM runs, so on the exact
invocation that corrupts, the reassurance is either false or unparseable. Qualifying both modes fixes
it either way.)

### 0.4 Placement — "refuses before anything is written" holds for every mode

Verified line numbers at `49152c1`: target site `src/main.rs:386`; `ensure_output_dir`
`src/main.rs:546`; dispatch `--hardtrim5` ≈ 597, `--hardtrim3` ≈ 653, `--clump_only` 713. So the site
precedes directory creation and every dispatch branch.

Structural argument confirmed empirically using the **sibling guard 17 lines earlier** (`src/main.rs:369`,
`--phred64` + uBAM input), which occupies the same position in the control flow. With a deliberately
non-existent `--output_dir`:

```
=== ht5        exit=1 output_dir_created=NO ===
=== clumponly  exit=1 output_dir_created=NO ===
=== trim       exit=1 output_dir_created=NO ===
=== stray outputs in scratch root === (none)
```

A bail at that position creates no output directory and no files in all three dispatch modes,
including the two that return early from `main()`. The plan's claim holds.

### 0.5 Mode completeness — the output-format gate covers every BAM-writing site

Every `cli.phred_offset()` call that reaches a `BamWriter` lives in `src/main.rs` (622, 653, 854, 917,
958, 2346, 2507, 2634) and sits inside an `OutputFormat::UBam` match arm or a function only reached
from one. `src/specialty.rs:150,202` and `src/clump_only.rs:836,1086` receive the offset as a
parameter from those sites. So gating on `output_format == UBam` is not an approximation of the
affected set — it *is* the affected set. The plan's mode list (plain trim SE/PE, `--hardtrim5/3`,
`--clump_only`) is complete; `--clock`, `--implicon`, `--demux`, `--passthrough`, `--retain_unpaired`
and `--clumpify` are all rejected in `Cli::validate()` §3.4a before this point.

### 0.6 `--passthrough` + uBAM output — the plan is right

`src/cli.rs:630` rejects `--passthrough` under `--output-format ubam` unconditionally, in
`validate()`, before any of this. So the `.chain(cli.passthrough.iter())` arm genuinely contributes
nothing new in the uBAM case. §8's "*Verified in-tree, not assumed*" checks out.

### 0.7 #362 regression guard — genuine Phred+64 → uBAM keeps working

Four invocations on `test_files/phred64_test.fastq` (lowest char `B` = 66, so above the ASCII-64
floor and untouched by the new condition):

| Invocation | Exit | `QUAL` |
|---|---|---|
| `--output-format ubam --phred64 -q 0 --length 0` | 0 | `III…###` → Phred 40 / Phred 2 ✅ |
| `--hardtrim5 20 --output-format ubam --phred64` | 0 | Phred 40 ✅ |
| `--clump_only --output-format ubam --phred64` | 0 | Phred 40 / 2 ✅ |
| `--paired --output-format ubam --phred64` | 0 | 8 records, Phred 40 / 2 ✅ |

`'h'` (104) − 64 = 40 and `'B'` (66) − 64 = 2, both correct. None of these is at risk from the
widened condition. The deliberate #362 behaviour survives.

### 0.8 Prior art — accurately read

| Ref | Actual | Plan's reading |
|---|---|---|
| #358 | issue, CLOSED — "`--phred64` + `--output-format ubam` writes spec-non-conformant quality scores"; writer hardcoded `saturating_sub(33)`, so Q40 stored as raw **71** | "inverse bug, 31 points too **high**" ✅ |
| #362 | PR, MERGED — "honour `--phred64` on BAM output; reject it for BAM input" | option A of #358's two; blanket refusal would take it back ✅ |
| #532 | issue, CLOSED — "Validate the `--phred64` claim against the input before writing output" | ✅ |
| #536 | PR, MERGED — "check the `--phred64` claim against the input under `--fastqc`" | ✅ |

The plan quotes #536's rationale as "refusing on `--phred64` alone would reject `--hardtrim5/3` and
`--clump_only` runs that pass scores through untouched". #536's body reads: *"refusing there would
also reject `--hardtrim5/3` and `--clump_only` runs that pass the scores through untouched, and that
is a wider behaviour change than a point fix should make."* Faithful.

Worth adding to §2, because it strengthens the case: **#532's own body carries the same error as the
NOTE.** It states "With `--clump_only` and no `--fastqc` nothing reads the scores, so the output is
byte-identical to the input and the flag is genuinely inert" — true for FASTQ output, false for
`--clump_only --output-format ubam`, as §0.1 row 5 shows. The mistake is not isolated to one
`eprintln!`; it is the premise the whole #532 → #536 scoping rested on. That is the cleanest
one-line justification for this plan.

### 0.9 Non-findings I checked and cleared

- **FIFO / process-substitution input.** A second bounded read pass could in principle break
  `<(zcat x.fq.gz)`. It cannot: Trim Galore already rejects FIFOs and `/dev/stdin` with a dedicated
  message ("Write the stream to a file first"), and already re-opens every input several times
  (format detection, sanity check, adapter auto-detection). No new constraint.
- **`--clump_only` FASTQ output really is quality-byte-identical** under `--phred64` — checked on the
  mixed fixture, qual line in == qual line out. The plan's premise for leaving FASTQ alone holds.
- **Line references.** Mostly exact (`main.rs:386`, `:369`, `:546`, `:553-564`; `cli.rs:605`, `:630`;
  `bam.rs:534,550,615`). Minor drift: `main.rs:26` is the doc comment, the const is at **32**;
  `main.rs:343` is **340/345**; the invariant block starts at **420**, not 415. Immaterial.

---

## 1. Logic review

The core reasoning is correct and I found no gap in it.

**The gate axis is right, and provably so.** §0.5 settles the question the plan only argues: output
format is not a proxy for "where the declared offset feeds disk" — every such site is behind
`OutputFormat::UBam`. That is a stronger claim than the plan makes for itself, and it is worth
stating in §8, because it is what rules out the "some other path reads `phred_offset()`" risk that
§11's "Remaining risks" raises as open.

**Ordering is sound.** The walk in §11 (343 → 369 → 386 → 415 → 546 → 553 → dispatch) matches the
file, and nothing between `cli.validate()` and line 386 writes: `sanity_check_any` and
`detect_input_format` are read-only. §0.4 confirms the consequence empirically.

**One reasoning error — §11 mislabels the BAM skip.** §11 says:

> the BAM skip inside the loop is belt-and-braces for the `--passthrough` arm rather than
> load-bearing — it is retained unchanged.

The skip **is** load-bearing, on the `--fastqc` + FASTQ-output path. `--passthrough` with a uBAM
carrier is legal at line 386 (`any_bam` covers only `cli.input`, and the carrier's own format check
is at `src/main.rs:456`, *after* the scan). I built a uBAM carrier whose records carry raw Phred 2,
so `BamReader` synthesises ASCII 35, and ran
`--paired --phred64 --fastqc --passthrough lowq_carrier.bam <genuine p64 mates>`:

```
exit=1  output_dir_created=NO
Error: --passthrough is not supported with uBAM input in this release. …
```

That is the correct message, reached only because the skip let the loop fall through to line 456.
Without the skip, `lowest_quality_char` would read ASCII 35 off the synthesised quality and bail at
386 with *"contains quality character '#' (ASCII 35)"* — naming a character the BAM does not contain,
which is exactly the failure the code comment at `src/main.rs:379-382` exists to prevent.

The plan **retains the skip**, so there is no implementation defect here. The problem is that the
plan records a load-bearing guard as decorative, in the one document a later refactor would consult.
Correct §11, and — better — put the fact where it constrains the code: a short note on the skip
itself, per CLAUDE.md's "put an invariant on the thing it constrains".

**Edge-case table is accurate.** Every row checks out against the code, including the two the plan
calls unreachable. The Solexa row (ASCII 59–63) is worth one extra clause: genuine Solexa/Illumina
1.0 data with `--output-format ubam` is *currently accepted* and floored to Q0 by `saturating_sub`;
after the change it is refused. That is an improvement — the current output is already garbage — but
it is a second class of newly-refused run, and §7's "Behaviour change" list does not name it.

---

## 2. Assumptions

All three stated assumptions hold, one with a correction:

1. **`--passthrough` rejected with uBAM output** — ✅ verified, `src/cli.rs:630` (§0.6).
2. **10,000-read bound stays; near-moot on uBAM** — ✅ the reasoning is right (real Phred+33 data
   carries sub-ASCII-64 bytes in essentially every record). Accepting the residual is the right call.
3. **`--hardtrim5/3` and `--clump_only` stay off §3.4a** — ✅ correct, and §0.7 confirms the reason:
   those combinations work correctly on genuine Phred+64 today, so a blanket exclusion would be a
   regression against #362.

**Unstated assumption worth surfacing:** that the FASTQ-output trimming path is "loud" enough to need
no protection. I tested this, because it is the only thing holding the scope where it is. On
`--phred64` FASTQ output at default `-q 20` with mislabelled input: **exit 0**, an empty
`*_trimmed.fq`, and a report that says `Quality-trimmed: 80 bp (100.0%)` and `Sequences removed
because they became shorter than the length cutoff of 20 bp: 4 (100.0%)`. The report is unambiguous,
so the plan's characterisation is defensible — but the **exit code is 0**, which is precisely the
complaint #532 opened with ("a pipeline gating on exit status reads a successful trim as a failed
step", inverted). §2's phrase "largely self-defends" overstates it: the run destroys the library and
reports success. I am not asking to widen the scope — that is a bigger change than this plan should
make — but §8 should record the residual honestly rather than leaving the FASTQ path described as
self-defending.

---

## 3. Efficiency

Nothing to object to. One avoidable duplicate, pre-existing from #536 and inherited here:

- The loop calls `detect_input_format(input)?` for every path, but `input_formats` was already
  computed for `cli.input` at `src/main.rs:340`. The `cli.input` arm could index that vector and
  detect only the passthrough carrier. One 24-byte peek (plus at most one BGZF block) per input
  saved — negligible, and the rewrite costs the clean `.chain()` shape, so **not worth doing** as
  part of this change. Noting it only because §6 enumerates the cost without noticing the value is
  already in scope.

§6's note about the duplicate 10,000-read pass on `--fastqc`-plus-uBAM runs is correct and the
decision not to merge it is right.

---

## 4. Validation sufficiency

This is where the plan needs the most work. The six table entries are well chosen and would catch a
wrong implementation; the problems are in the mutation prescription and in what the table omits.

### 4.1 One prescribed mutation cannot fail the test it is assigned to — **Important**

§9 requires:

> flip `lowest < 64` to `lowest <= 64` → **test 4** must fail.

It cannot. Test 4's fixture is `fixture_copy(dir, "good_p64.fastq")`, which copies
`test_files/phred64_test.fastq`, whose distinct quality characters are `B` (66) and `h` (104). With a
minimum of 66, `lowest <= 64` is still false, so the mutant does not bail and test 4 passes. An
implementer running this loop would see the mutant survive and conclude the new tests are
insufficient — chasing a gap that is not there.

The boundary **is** covered, by a test the plan does not mention:
`tests/integration_fastqc_phred64.rs:508`, `preflight_accepts_ascii_64_which_is_phred64_q0`, which
builds a fixture whose lowest character is `'@'` (ASCII 64) and asserts success. It uses `--fastqc`
with FASTQ output, so it stays green after the change and still catches the `<=` mutant. Fix: point
the third mutation at that test instead, and add it to §9 row 5 as a second "stays green" guard.

While there: `good_p64.fastq` is a *destination filename*, not a fixture — `test_files/` contains no
such file. §9 row 4's "genuine Phred+64 fixture (`good_p64.fastq`, lowest `B`=66)" reads as if it
were one. Say "a copy of `test_files/phred64_test.fastq` (lowest `B` = 66) via `fixture_copy`".

### 4.2 §5.3 is declared a required correctness fix and has no validation entry — **Important**

§5.3 states the NOTE is *"wrong"* and the correction is *"required regardless of the rest of this
plan"*; §11 item 2 repeats that it is "a required correctness fix in its own right, not documentation
polish". The §9 table has **no row for it**. Nothing in the suite asserts anything about the NOTE's
text, so a mis-edit — or simply forgetting it, since it is the only item in §5 with no test attached —
ships silently. The NOTE is runtime stderr on a user-visible path, not prose; it is as testable as
the refusal message that row 6 covers.

Add a row: on a run that reaches the NOTE under uBAM output (e.g. the genuine-Phred+64 invocation of
row 4), assert the stderr does **not** claim that `--hardtrim5/3` passes scores through untouched.
That is one `assert!(!stderr.contains(...))` in the test row 4 already needs.

### 4.3 Mutation coverage after the change — mostly good, one note

Working through what each mutant breaks:

| Mutant | Caught by | Status |
|---|---|---|
| drop `\|\| matches!(… UBam)` | new rows 1–3 | ✅ as the plan says |
| `\|\|` → `&&` | new rows 1–3 **and** existing `preflight_refuses_mislabelled_input_before_writing_anything` | ✅ doubly |
| drop `cli.fastqc_requested()` arm | existing `preflight_refuses_mislabelled_input_before_writing_anything` (`--fastqc`, FASTQ out) | ✅ — worth stating, since the plan does not |
| drop the condition entirely (always scan) | existing `preflight_does_not_fire_without_fastqc…` (row 5) | ✅ — this is what #536's PR body calls "the scope test", and renaming it must not weaken it |
| `lowest < 64` → `<= 64` | `preflight_accepts_ascii_64_which_is_phred64_q0` (**not** row 4) | ⚠️ see §4.1 |

### 4.4 Smaller coverage gaps — **Optional**

- **`--hardtrim3` is untested.** Row 2 covers `--hardtrim5` only. `hardtrim3_to_bam` is a separate
  call site (`src/main.rs:653`). Both are downstream of 386, so one test is defensible — but say so,
  rather than writing "`--hardtrim5/3`" in the prose and testing one of them.
- **No paired uBAM-output case.** `--paired --output-format ubam` is the single-interleaved-BAM
  output shape; the loop is shared, so this is low risk, and `preflight_scans_every_input_not_just_the_first`
  covers the multi-input iteration on the `--fastqc` path. Acceptable to skip; worth one line in §9
  saying it is skipped deliberately.
- **An existing test already covers genuine-p64 → uBAM** under `--fastqc`:
  `bam_format_encoding_wins_over_the_phred64_flag` (`tests/integration_fastqc_phred64.rs:95`) runs
  `--phred64 --fastqc --output-format ubam` and asserts per-base means 40.0 and 2.0. After the change
  it exercises *both* arms of the widened condition and must stay green. Cite it in row 5 — it is a
  stronger #362 guard than row 4 alone, and it is free.

---

## 5. Documentation and convention issues

### 5.1 `docs/guide/outputs.md` is the second docs surface, and §5.7 misses it — **Important**

`docs/src/content/docs/guide/outputs.md:90-100` carries the canonical list, headed **"Rejected at
startup, before anything is written (v1 scope)"**, of what `--output-format ubam` refuses:
`--dont_gzip`, `--clock`/`--implicon`, `--demux`, `--passthrough`, `--retain_unpaired`, `--clumpify`,
`--rename`. That list already includes a *conditional* refusal in the same shape as this one
(`--rename` — "only when at least one input is FASTQ, and only when a clipping or hardtrim flag …
is set"), so there is precedent and a natural slot.

§5.7 touches only `guide/quality.md`. A reader who wants to know what uBAM output refuses will read
`outputs.md`, not the Phred-encoding section of the quality page, and will not find the new refusal.
Add a bullet, cross-referencing `/guide/quality/#phred-encoding` rather than restating the rationale.

I swept the docs tree: `quality.md` (lines 31, 33, 35, 37, 40, 57) is the only page that mentions
`--phred64`, and of those only 37 and 40 become stale, exactly as §5.7 says. Line 35 is about BAM
*input* and stays true. So the quality.md half of §5.7 is complete — it is only the uBAM-output page
that is missing.

### 5.2 The proposed source comment breaks the rule it cites — **Optional, but fix before implementing**

§5.2 says "Two lines maximum", then gives a three-line example — and the example is a *justification
of the change* ("on FASTQ output the flag only says how to read bytes that are copied verbatim"),
which CLAUDE.md puts in the commit message, not the source. It would also sit on top of the existing
five-line comment at `src/main.rs:381-385`, making an eight-line block above a two-line condition.

Per CLAUDE.md ("default to a single line… state the fact, not the evidence"), one line stating what
constrains the code is enough, e.g.:

```rust
// On uBAM output the declared offset is subtracted into QUAL and floors at 0.
```

Everything else belongs in the commit message. This matters more than usual here: the project's
memory records five rounds of maintainer pushback on exactly this.

### 5.3 §5.6's CHANGELOG instruction is self-contradictory — **Optional**

"one sentence, matching the file's format" — but the file's format is **not** one sentence. Every
`### Version 2.4.0` entry is a bold lede plus a 3–5 sentence paragraph; the #532 entry this plan
amends (`CHANGELOG.md:21-27`) runs seven lines. An implementer cannot satisfy both halves. Decide
which governs and say so. My read: match the surrounding entries' shape (bold lede + short
paragraph) but keep it to the shortest in the file, and keep the "point at the docs" rule — that
honours the intent of the CLAUDE.md rule without making one entry conspicuously unlike its neighbours.

Two smaller points on §5.6:
- The amendment target is correct: `CHANGELOG.md:27` does end "`--phred64` without `--fastqc` is
  unaffected", and that becomes false. ✅
- The file has an empty `### Unreleased` section at line 4 above `### Version 2.4.0` at line 6. The
  plan says to add under 2.4.0 without saying why. It is the right choice (v2.4.0 is genuinely
  unreleased — confirmed: no `v2.4.0` tag, latest is `v2.3.0`, `master` is 3 commits behind `dev`,
  all as the plan's header states), but one clause saying so saves the implementer a decision.

### 5.4 §10's release-date flag is real

Confirmed: `CHANGELOG.md:6` reads "Version 2.4.0 (Release on 3 November 2026)" and `CITATION.cff:18`
has `date-released: 2026-11-03`, while no `v2.4.0` tag exists. Correctly scoped out of this plan.

---

## 6. Alternatives

The plan considered four options and picked the right one. Two things worth adding to §8 rather than
reopening:

**The alternative the plan half-considered and should name explicitly.** Gating on *what will read
the scores* — `phred64 && (fastqc_requested() || ubam_out || quality > 0 || nextseq.is_some())` —
rather than on output format. §2's data is what argues against it: the FASTQ trim path's failure is
loud in the report (§2 above), so refusing it buys little and breaks `-q 0` FASTQ runs that are
genuinely inert. But the honest framing is "deferred because the FASTQ failure is visible, though it
still exits 0", not "the FASTQ path self-defends". Naming it as a considered-and-deferred option
makes the eventual #532 close-out easier.

**`--phred64` on §3.4a's exclusion list is correctly rejected.** This was #358's own option B, and
§0.7 shows why it is wrong now: four genuine-Phred+64 → uBAM invocations work correctly today. §8
states this; it is the single most important thing in the plan to keep, because it is the one
decision a reviewer unfamiliar with #362 would be tempted to overturn.

---

## 7. Action items

### Critical
None. The change as specified is correct, correctly sited, complete in its mode coverage, and does
not regress #362.

### Important
1. **Fix the third prescribed mutation (§9).** `lowest < 64` → `<= 64` cannot fail row 4 — that
   fixture's minimum is ASCII 66. Point it at the existing
   `preflight_accepts_ascii_64_which_is_phred64_q0` (`tests/integration_fastqc_phred64.rs:508`),
   which does catch it, and list that test under row 5 as a "stays green" guard. (§4.1)
2. **Add a validation row for §5.3.** The NOTE correction is declared required and has no check.
   One `assert!(!stderr.contains(…))` inside the row-4 test covers it. (§4.2)
3. **Correct §11's "belt-and-braces" claim.** The BAM skip inside the loop *is* load-bearing for the
   `--passthrough` arm on the `--fastqc` + FASTQ-output path — without it a uBAM carrier gets refused
   at 386 naming a character the file does not contain, pre-empting the correct message at
   `src/main.rs:456`. Demonstrated in §1. The code is right; the rationale is not, and §11 is where a
   future refactor will look. (§1)
4. **Extend §5.7 to `docs/guide/outputs.md:90-100`,** the canonical uBAM "rejected at startup" list,
   which already carries a conditional refusal (`--rename`) in the same shape. (§5.1)

### Optional
5. Rewrite §5.2's example comment as one line stating the fact, not three justifying the change; it
   currently exceeds the limit the same section sets. (§5.2)
6. Resolve §5.6's "one sentence / matching the file's format" contradiction, and add a clause on why
   the entry goes under `### Version 2.4.0` rather than the empty `### Unreleased`. (§5.3)
7. In §8, replace "the FASTQ trim path largely self-defends" with the measured residual: the run
   destroys the library, says so clearly in the report, and still **exits 0**. Name the
   gate-on-what-reads-the-scores alternative as considered-and-deferred. (§2, §6)
8. Name the Solexa (ASCII 59–63) case in §7's behaviour-change list — genuine Solexa + uBAM output is
   accepted-and-floored today and will now be refused. An improvement, but a second class of new
   refusal. (§1)
9. §9 row 4: `good_p64.fastq` is a destination name used by `fixture_copy`, not a file in
   `test_files/`. Say "a copy of `test_files/phred64_test.fastq` (lowest `B` = 66)". (§4.1)
10. Say in §9 that `--hardtrim3` and paired uBAM output are deliberately not given their own rows,
    rather than leaving the prose and the table disagreeing. (§4.4)
11. Fix the three drifted line references: the const is `src/main.rs:32` (26 is its doc comment),
    `input_formats`/`any_bam` are at 340/345 (not 343), and the invariant block starts at 420 (not
    415). (§0.9)
12. Add #532's own matching error to §2 — its body claims `--clump_only` under `--phred64` is
    "genuinely inert", which is false for `--clump_only --output-format ubam`. The mistake is in the
    premise the #532 → #536 scoping rested on, not just in one `eprintln!`, which is the tightest
    justification this plan has. (§0.8)

---

## 8. Note on what I did not run

I did not build a mutated binary. The two questions a build would settle were settled otherwise: the
"nothing is written" property by exercising the sibling guard at `src/main.rs:369`, which occupies
the same control-flow position (§0.4), and the `<=` mutation by reading the fixture's actual quality
alphabet (§4.1) — a build would have confirmed the mutant survives row 4, not changed the conclusion.
The working tree was not modified; all scratch output is under
`/tmp/claude-501/revA-phred64-1794048757/`.
