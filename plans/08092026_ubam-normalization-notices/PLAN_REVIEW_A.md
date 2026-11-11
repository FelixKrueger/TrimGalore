# Plan review A — fix reversed IUPAC warning + dropped-description notice (#406)

**Plan:** `plans/08092026_ubam-normalization-notices/PLAN.md`
**Base verified against:** `dev` @ `90dda11` (verified via `git log --oneline -1`)
**Reviewer:** A (independent; Reviewer B reviewed the same plan separately)
**Method:** static verification against `src/`, `tests/`, `docs/`, `.github/` **plus**
an empirical build. The tree was copied with `git archive HEAD` into a scratch
dir with a private `CARGO_TARGET_DIR` — the shared working tree was never
mutated and no branch was switched.

**Verdict in one line:** the plan's two chosen actions are both right and should
proceed, but **one of its load-bearing arguments (A-in-§Context, the
"nearly unreachable" shared-`OnceLock` claim) is factually wrong**, its
read-side "already correct" claim is wrong in one reachable configuration, and
the drafted description notice contains a clause that is false in exactly the
case where it fires. Three findings below are Critical.

---

## 1. Logic review

### 1.1 The two call sites and their directions — VERIFIED (claim 1)

`grep -rn emit_iupac_warning_once src/ tests/` returns exactly three hits:

| Line | Enclosing fn | Direction |
|---|---|---|
| `src/bam.rs:823` | `validate_and_normalize_seq_for_write` | FASTQ → uBAM (write) |
| `src/bam.rs:944` | `bam_record_to_fastq` | uBAM → FASTQ (read) |
| `src/bam.rs:995` | the definition (shared `static SEEN: OnceLock<()>` at `:996`) | — |

`validate_and_normalize_seq_for_write` has exactly one production caller,
`BamWriter::write_record` (`src/bam.rs:607`), so its input is always a
`FastqRecord` and its output always a BAM record. The reversal is real.

**Confirmed empirically, not just by reading.** A two-record FASTQ carrying
`R`/`Y` bases, run as `trim_galore --output-format ubam`:

```
Trimming: sampleF.fastq
Output:   out1/sampleF_trimmed.bam
WARNING: input uBAM contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V);
coerced to N for FASTQ output. …
```

There is no uBAM input and no FASTQ output in that run. Both halves are
backwards, exactly as the plan states. This also doubles as the plan's
validation 8 (pre-patch expected-fail control) — it is satisfied, and the recipe
is recorded in §4.1 below so the implementer does not have to re-derive it.

### 1.2 CRITICAL — the shared-`OnceLock` "nearly unreachable" argument is wrong (claim 2)

The plan (§Context, second bullet) argues:

> In practice the cross-direction suppression is nearly unreachable: on
> uBAM→uBAM the read side coerces first, so the write side sees only `ACGTN`
> and `iupac_seen` stays false.

The uBAM→uBAM half of that is correct and I verified it (see 1.3). **But it is
not the only way to reach both sides in one process.** Mixed FASTQ + uBAM
inputs in a single **single-end multi-input** run reach both, and such runs are
legal and explicitly supported:

- `any_bam` at `src/main.rs:240-242` is an `.any(…)` over per-input formats —
  mixed input is tolerated by construction, not merely unguarded.
- The design *depends* on mixed input being legal: `src/main.rs:296` tells the
  user to "convert **at least one** input to uBAM via `samtools import`" as the
  remedy for `--preserve-tags` with all-FASTQ inputs.
- The SE uBAM-output loop (`src/main.rs:1987-1994`) iterates `cli.input` and
  re-detects the format **per input** (`src/main.rs:2016`).
- The BAM/FASTQ mismatch rejection (`reject_bam_format_mismatch_in_pair`,
  `src/main.rs:347-362`) is gated on `cli.paired` — it does not apply to SE.

Empirically, with one IUPAC-bearing FASTQ and one IUPAC-bearing uBAM (distinct
stems, so the collision pre-flight does not fire):

```
# order A: FASTQ first
Trimming: sampleF.fastq
WARNING: input uBAM contains IUPAC degenerate bases …    ← wrong direction named
Trimming: sampleB.bam                                     ← correct message SUPPRESSED

# order B: BAM first
WARNING: input uBAM contains IUPAC degenerate bases …    ← correct
Trimming: sampleB.bam
Trimming: sampleF.fastq                                   ← write-side coercion SUPPRESSED
```

Order A is the worse outcome and it is the natural argument order: the user is
told the wrong direction *and* never told about the genuine uBAM coercion.

**Why this matters for the plan.** It does not change what the plan *does* —
splitting into two independent `OnceLock`s fixes this too. It changes three
things:

1. The severity. Splitting the static is not a tidiness improvement ("worth
   removing anyway", "two independent statics make each direction's guarantee
   local"); it repairs a reachable suppression bug. That belongs in the
   CHANGELOG bullet, which currently plans to mention only the reversed text.
2. Validation 3 is no longer hedged — see §4.1.
3. The plan should stop asserting unreachability, because a future reader will
   trust it.

### 1.3 CRITICAL — the read-side text is NOT "already correct" (claim 1, Behavior 2)

The plan states the read-side message is accurate and needs only a new guard
(Behavior 2, A3, outline item 3). That holds for uBAM → FASTQ. It does **not**
hold for uBAM → uBAM, which fires the same message:

```
$ trim_galore --output-format ubam -o out3 iupac.bam
WARNING: input uBAM contains IUPAC degenerate bases …; coerced to N for
FASTQ output. …
Output:   out3/iupac_trimmed.bam
```

The coercion is real (it happens into the internal FASTQ intermediate), but the
phrase **"for FASTQ output"** names an output that does not exist in this run —
the user's output is a BAM. So the plan would ship a fix for one misleading
direction claim while leaving a second one in place, in the *same* message, on a
path #400 already had to re-document.

The fix is small and makes the message unconditionally true: name the coercion,
not the destination. Suggested:

```
WARNING: input uBAM contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V);
coerced to N on read. This is legitimate in PacBio HiFi, ONT Dorado, and 10x
cellranger uBAMs.
```

"on read" is true for uBAM→FASTQ and uBAM→uBAM alike. Since A1 holds (nothing
pins the text — see §2.1), this costs nothing.

### 1.4 Description notice scope — VERIFIED (claim 3)

- `parse_name_and_data` (`src/bam.rs:707`) has exactly **one** production caller,
  `BamWriter::write_record` (`src/bam.rs:589`); every other hit is in the test
  module (`:1226`–`:1392`). Write-side only. ✅ A2 verified.
- All three uBAM output families route through `BamWriter::write_record`:
  trim (`trimmer.rs:678`, `:790`, `:791`), specialty hardtrim
  (`specialty.rs:170`, `:223`), clump-only (`clump_only.rs:694`, `:706`,
  `:707`). `BamWriter::create` call sites confirm the same three families
  (`main.rs:2025/2184/2316`, `specialty.rs:152/204`,
  `clump_only.rs:765/1003`). ✅ The notice is not clump-specific, as claimed.
- The drop is confirmed at runtime: `@readA 1:N:0:ACGTAC` → BAM `QNAME=readA`.

### 1.5 CRITICAL — `--rename` is reachable with uBAM output and its annotation is silently destroyed

The review brief asks whether `--rename` interacts. It does, worse than
"the notice might fire on IDs TrimGalore constructed".

`--rename` is **not** in the `--output-format ubam` exclusion list. `cli.rs:582-622`
(§3.4a) rejects `--dont_gzip`, `--clumpify`, `--passthrough`, `--clock`,
`--implicon`, `--demux`, `--retain_unpaired` — no `--rename`. It is rejected only
against `--clump_only` (`cli.rs:846-848`).

`FastqRecord::append_to_id` (`src/fastq.rs:173-181`) splices the
`:clip5:`/`:clip3:` suffix in **before the first TAB** — deliberately, per the
code-review C1 finding cited at `:171` — but has **no equivalent handling for a
SPACE**. So on a header with a description the suffix is appended after the
description, and `parse_name_and_data` then discards the whole tail.

Measured, with `--rename --clip_R1 3`:

| Input header | FASTQ output id | uBAM output QNAME |
|---|---|---|
| `@withspace 1:N:0:ACGTAC` | `@withspace 1:N:0:ACGTAC:clip5:ACG` | `withspace` — **`:clip5:ACG` gone** |
| `@nospace_read` | `@nospace_read:clip5:ACG` | `nospace_read:clip5:ACG` — survives |

`--rename` exists specifically so the clipped bases can be recovered later, so
this is a silent loss of the very thing the flag was asked to record — and it is
data-dependent (same flags, same command, loss or no loss depending on whether
the header has a space). It is strictly worse than losing an instrument ID.

This is arguably a separate bug, and I am **not** proposing the plan fix it —
that would be scope creep past Felix's decision. But it has two consequences
inside this plan's scope:

1. The drafted notice's wording is wrong for this case: it says *header
   descriptions* are not carried over, when what was also dropped is
   TrimGalore's own annotation. A user who sees the notice and does not care
   about instrument IDs will dismiss it and never learn their `--rename` record
   is gone.
2. It should be filed as a follow-up issue now, while the analysis is fresh, and
   referenced from the plan. The obvious candidates are teaching
   `append_to_id` to splice before the first SPACE too (matching its existing
   TAB behaviour), or adding `--rename` to §3.4a the way `--clump_only` already
   rejects it.

### 1.6 IMPORTANT — Behavior 6 misses two edge cases, one of which makes the notice's own text false

The plan's edge-case enumeration covers no-whitespace, TAB tail, and
TAB-then-space-inside-tag-value. I verified the third is handled correctly:
`@readY\tCB:Z:has space here` round-trips as `CB:Z:has space here` intact. ✅

Two cases are missing.

**(a) Trailing space / empty description → false-positive notice.**
`FastqReader` strips only `\n` and `\r` from the id line
(`src/fastq.rs:412`, `:446`: `trim_end_matches(['\n', '\r'])`), so a trailing
space **survives into `record.id`**. The codebase already contemplates this —
`src/fastq.rs:680` tests `append_to_id` on `id: "@read1  \n"`. For
`@trailspace ` the SPACE branch of `parse_name_and_data` is taken with
`rest == ""`, and confirmed output `QNAME=[trailspace]`.

An implementation that keys naively on "the SPACE branch was taken" — which is
literally what Behavior 6 and outline item 4 specify ("the notice keys strictly
on the SPACE branch being taken") — will announce a dropped description when
nothing was dropped. **Guard on `!rest.is_empty()`**, and add a unit test for
`@read1 ` alongside the existing `:1249` case.

**(b) SPACE before TAB → aux tags are silently lost, so the notice's third
clause is false.** The split takes the *first* whitespace, so a description
followed by a tag tail loses the tags too:

```
input header : @readX desc<TAB>CB:Z:ATCG
BAM record   : readX  4  *  0  255  *  *  0  0  ACGT…  IIII…
                                                   ↑ no CB tag at all
```

The drafted notice ends with **"Tab-separated TAG:TYPE:VALUE fields are
preserved."** In this input — a case where the notice *does* fire — that
sentence is false. This is the plan repeating the class of defect it exists to
fix: a message that is true for the common case and misleading for a reachable
one. Either qualify it ("…are preserved when they follow the read name
directly") or drop the clause.

Reachability: not from uBAM input (BAM QNAME cannot contain a space, so
`bam_record_to_fastq` always produces `@name\tTAG:…`), so this needs a
hand-built or third-party FASTQ. Rare, but the notice is a correctness
statement, and `write_record` parses the tail unconditionally — it does not
depend on `--preserve-tags`.

### 1.7 A4 — the near-universal-firing premise, and my judgement (claim 4)

**The two fixtures the plan names do carry descriptions**, as claimed, and the
dropped half is indeed the original instrument identifier:

```
BS-seq_10K_R1.fastq.gz        @SRR24827378.1 A00690:224:H7GC5DRX2:1:2101:1081:1016/1
SRR24766921_RRBS_R1.fastq.gz  @SRR24766921.1 A00686:91:HTHYKDMXX:1:1101:19696:1000/1
```

But "nearly every real FASTQ→uBAM run" overstates it. Surveying the first header
of all 20 `test_files/` FASTQ fixtures: **9 have no whitespace at all** and would
never fire the notice — `BS-seq_10K_I1` (`@SRR24827378.1`), `10K_150bp`,
`4_seqs_with_Ns`, `illumina_10K` (`@MS2_16089:1:1101:1740:15559#45/1`),
`illumina10K_with_polyA` (`@1`), `PolyA`, `PolyT`, `polyAT_R1/R2` (`@seq1`),
`truncated`. The classic non-SRA Illumina `#45/1` style has no space either.

**My judgement: fire it, once per run, and the plan's instinct is right.**
Reasons, in order of weight:

1. One extra line per process is not a noise problem, and the FASTQ→uBAM path
   *already* prints an unconditional informational line —
   `"NOTE: input is FASTQ; output uBAM records will have no aux fields."`
   (`src/main.rs:2021`). A second caveat on the same path is consistent with
   established behaviour rather than a new burden. Anyone who found one
   tolerable will find two tolerable.
2. The loss is irreversible and, per 1.5, is not limited to archival metadata.
3. Suppressing it behind a verbosity flag (the A4 fallback) is the worse
   failure: the users who most need it are exactly those who did not think to
   ask.

Two refinements I would apply, though:

- **Not `WARNING:`.** See §3.1.
- **Consider co-locating it with the `main.rs:2021` NOTE** rather than firing
  lazily from deep inside the writer. Trade-off: eager emission at startup
  groups the FASTQ→uBAM caveats into one readable block and cannot land
  mid-progress-output, but it over-discloses (it would fire for the 9
  whitespace-free fixtures too). Lazy emission is precise but interleaves.
  **I recommend keeping the plan's lazy approach** — precision matters here
  precisely because ~half of real headers have nothing to drop — but the
  implementer should know the existing NOTE is right there and is stylistically
  the closest sibling.

### 1.8 Minor — per-process vs per-file inconsistency

The adjacent `main.rs:2021` NOTE fires **per input file** (it sits inside
`run_ubam_output_single`, called from the SE loop). The new notice is specified
as **once per process**. In a 20-input SE run the two informational lines will
have visibly different cardinality. Once-per-process is defensible (that is what
the issue asked for) — just flagging that the neighbouring precedent differs, in
case Felix would rather they match.

---

## 2. Assumptions

### 2.1 A1 (nothing pins either message) — VERIFIED, with a caveat on method

Repo-wide grep for `IUPAC` / `coerced to N` / `degenerate` across `src/`,
`tests/`, `docs/`, `.github/`, `README.md`, `CHANGELOG.md` finds the message text
only in `src/bam.rs`. No test, no CI grep, no docs copy. **Rewording is free.** ✅

Caveat: three places describe the *behaviour* and would become stale if the
reworded message contradicted them — `src/cli.rs:187`, `CHANGELOG.md:163`, and
`docs/src/content/docs/modes/clump-only.md:52`. None quotes the message, so none
breaks; but the plan's greps were on exact phrases (`"IUPAC degenerate"`,
`"coerced to N"`) and would have missed all three, which say "IUPAC codes become
`N`" / "are coerced to `N`". The conclusion is right; the search was narrower
than the claim.

### 2.2 A5 (no test asserts stderr line counts) — CONCLUSION RIGHT, JUSTIFICATION WRONG

The plan says "no test asserts on total stderr line counts, so adding a line
cannot break an existing assertion." Two things do exist:

- **An exhaustive stderr assertion.** `tests/integration_no_args_help.rs:42`
  asserts `stderr.is_empty()`. It is on the **bare-invocation** path (no args →
  help), which never opens a BAM writer, so it cannot break.
- **A count-based CI grep.** `.github/workflows/ci.yml:457`:
  `test "$(grep -c 'Auto-detecting adapter type' /tmp/detect.log)" = "2"`.
  It counts a *specific substring*, not lines, so a new line cannot break it.

Every other stderr assertion I found is `contains` / `!contains` on a specific
string: `integration_non_restartable_input.rs` (`REJECTION` =
`"cannot be re-read from the start"`, plus negative checks for
`"valid UTF-8"`, `"doesn't seem to be in FastQ format"`,
`"Failed to decompress first block"`), `integration_ubam.rs:165/254`,
`integration_ubam_out.rs:329/462/485/507/532/714/737`. None of those substrings
appears in either drafted message.

So **A5 holds** — but restate it as "the one exhaustive stderr assertion is on
the bare-invocation path, and the one count assertion counts a specific
substring", because the current phrasing would let a future change to
`integration_no_args_help.rs` invalidate the assumption invisibly.

### 2.3 Unstated assumption worth surfacing

The plan assumes the notice's trigger and the *actual* drop are the same event.
Per §1.6(a) they are not (empty description) and per §1.6(b) the drop is larger
than the notice describes (tags). Both are implementation-level, both are cheap
to get right, and neither is currently written down.

---

## 3. Message wording — recommendations on both open questions

### 3.1 Prefix: use `NOTE:` (uppercase), not `Note:` and not `WARNING:`

This is settled by convention, not taste. Counting string literals in `src/`:

| Prefix | Occurrences |
|---|---|
| `WARNING:` | 10 |
| `NOTE:` | 5 |
| `Warning:` | 1 |
| `Note:` | **0** |

The plan's drafted `Note:` is used nowhere in the codebase. And the closest
sibling by meaning — the other informational line on the same FASTQ→uBAM path,
`src/main.rs:2021` — is a `NOTE:`. So: `NOTE:`, uppercase. It also gets the
register right, which was the substance of the plan's question: a format-mandated
transformation that fires on most runs is a note, not a warning.

### 3.2 Keep `1:N:0:INDEX` — and the plan's reason for doubting it is factually wrong

The plan hesitates because "the fixtures' dropped tails are instrument
identifiers, not index fields". That is true only of the two SRA fixtures it
looked at. **Five fixtures have dropped tails of exactly the `1:N:0:INDEX`
shape:**

```
clock_10K_R1.fastq.gz     @HWI-D00436:407:CCAETANXX:1:1101:1832:1891 1:N:0: CGATGTTT
clock_10K_R2.fastq.gz     @HWI-D00436:407:CCAETANXX:1:1101:1832:1891 3:N:0: CGATGTTT
demux_test.fastq.gz       @INST001:001:FLOWCELL1:1:1101:18269:1000 1:N:0:CTTNTA
nextera_100K.fastq.gz     @HWI-D00436:149:C6HE7ANXX:1:1101:1394:1945 1:N:0: TAAGGCGA_TATCCTCT
smallRNA_100K.fastq.gz    @HWI-M02293:93:000000000-ADUH0:1:1101:14115:1000 1:N:0:
```

Decisively: **the docs already use that exact example for this exact
behaviour.** `docs/src/content/docs/modes/clump-only.md:51`:

> Read ID (header line for FASTQ; name portion for BAM — a space-separated
> description such as `1:N:0:INDEX` is dropped on uBAM output)

Keeping `1:N:0:INDEX` aligns the stderr message with the documented wording.
Cutting it would gratuitously desynchronise them.

(Aside: `clock_10K_R1`'s description contains its own space — `1:N:0: CGATGTTT`.
Harmless, the first space wins and the whole tail goes, but it shows "multiple
spaces" is a real fixture shape, not hypothetical.)

### 3.3 Consolidated wording suggestions

Write side (replaces the reversed text):

```
WARNING: input FASTQ contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V);
coerced to N in the uBAM output. This is legitimate in PacBio HiFi, ONT
Dorado, and 10x cellranger data.
```

The plan's draft, unchanged — it is accurate and I have no improvement.

Read side (per §1.3, "for FASTQ output" → "on read"):

```
WARNING: input uBAM contains IUPAC degenerate bases (R/Y/M/K/S/W/B/D/H/V);
coerced to N on read. This is legitimate in PacBio HiFi, ONT Dorado, and 10x
cellranger uBAMs.
```

Description notice (`NOTE:` per §3.1; final clause qualified per §1.6(b)):

```
NOTE: FASTQ header descriptions are not carried into uBAM output — the BAM
spec forbids whitespace in QNAME, so text after the first space is dropped
(e.g. an Illumina "1:N:0:INDEX" field). Tab-separated TAG:TYPE:VALUE fields
are preserved when they follow the read name directly.
```

---

## 4. Validation sufficiency

### 4.1 Validation 3 IS constructible — stop hedging, and use this recipe

The plan writes validation 3 as "a run that would coerce on both sides, **if
constructible**". It is. Verified end to end:

```bash
# 1. FASTQ with IUPAC + a space-separated description
printf '@readA 1:N:0:ACGTAC\nACGTRYACGTACGTACGTACGTACGTACGTACGT\n+\n%s\n' \
  'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII' > sampleF.fastq

# 2. uBAM with IUPAC — BAM's 4-bit nybble encoding stores all 16 codes,
#    and samtools import does NOT coerce, so K/W survive
printf '@readC\nACGTKWACGTACGTACGTACGTACGTACGTACGT\n+\n%s\n' \
  'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII' > src.fastq
samtools import -0 src.fastq -o sampleB.bam

# 3. one mixed single-end run — BOTH directions must now warn
trim_galore --output-format ubam -o mixA sampleF.fastq sampleB.bam
trim_galore --output-format ubam -o mixB sampleB.bam sampleF.fastq   # both orders
```

Two harness traps I hit, worth passing on:

- **Distinct stems are mandatory.** My first attempt named both inputs
  `iupac.*`; they collapse to one `iupac_trimmed.bam` and the output-collision
  pre-flight aborts the run before any trimming. The test looked like it ran.
- **Do not grep for a term that appears in your filenames.** With files named
  `iupac.*`, `grep -ci IUPAC` matched the *collision error's* pathnames and
  reported a plausible-looking count from a run that never trimmed a record.
  Grep for `degenerate bases`.

Post-patch expectation: two distinct messages in both orders, each once.
Pre-patch (control): exactly one, and in order A it is the wrong direction.

### 4.2 Validation 6 — right check, but use the existing harness

Byte-comparing output is the right instinct (the change touches a write path).
The repo convention is not a manual byte compare, though: it is the Rust
comparator `assert_ubam_eq` in `tests/integration_ubam_out.rs:108`, which
compares `header_minus_pg` plus per-record tuples, against the checked-in
goldens `test_files/ubam_out_se_REFERENCE.bam` and `ubam_out_pe_REFERENCE.bam`
(`:154`, `:171`). Its doc comment at `:102-104` explicitly says it is *not* a
`samtools view -H | grep -v @PG` pipe-diff, for CI stability.

So validation 6 is simply:

```bash
cargo test --release --test integration_ubam_out
```

**Verified usable:** I ran it on the unpatched tree — `23 passed; 0 failed`, in
0.90s. So it is a working pre-patch baseline, not a hoped-for one, and the
implementer can diff against a known-green starting point.

That is stronger than an ad-hoc before/after diff (it pins against committed
goldens rather than against a build the implementer produced) and needs no new
fixtures. Recommend replacing the manual comparison with this.

### 4.3 Gaps — the silent-wrong-result modes are not covered

The validation table catches everything about *whether the messages fire*. It
catches nothing about the two silent-wrong-result modes in the blast radius:

| Uncovered | Why it matters | Suggested check |
|---|---|---|
| Trailing-space / empty description (§1.6a) | notice fires claiming a drop that did not happen | unit test `parse_name_and_data("@read1 ")` → no notice; assert the name is `read1` |
| SPACE-then-TAB tag loss (§1.6b) | the notice's own final clause is false here | unit test `@readX desc\tCB:Z:ATCG` → tags absent, documenting the behaviour the message must not contradict |
| `--rename` annotation loss (§1.5) | silent loss of user-requested data, out of plan scope but adjacent | follow-up issue, not a test in this plan |

Validation 5 ("does not fire when nothing is dropped") is where (a) belongs; as
written it names "fixture whose ids have no space (and one with a TAB tail)",
which misses the empty-description case entirely. Add it there.

One genuine strength worth calling out: **validation 8 (run the unpatched build
first) is the right instinct and it paid off** — it is how I found §1.3, because
seeing all three runs pre-patch is what exposed that uBAM→uBAM also emits "for
FASTQ output". Keep it.

---

## 5. Efficiency

Nil, and the plan's reasoning is sound. Both notices sit behind `OnceLock`, and
`parse_name_and_data` already computes the split, so no new scan and no new
allocation. The description notice is reached per record on space-bearing input
(most real Illumina data), but `OnceLock::get_or_init` on an initialised cell is
a single atomic load plus a branch — immeasurable against BGZF compression,
which dominates this path and is single-threaded here by design.

One pedantic correction: the plan says "one relaxed atomic load". `OnceLock::get`
issues an **Acquire** load, not Relaxed. On x86-64 and AArch64 an Acquire load is
a plain load instruction, so the performance claim is unaffected — but if that
sentence survives into a code comment it should not say "relaxed".

No behavioural coupling: nothing branches on whether a message fired, so the
statics cannot influence output. That is what makes §1.2 a message bug and not a
correctness bug.

---

## 6. Alternatives considered

1. **One emitter taking a direction enum, vs two emitters.** The plan offers
   both and picks two functions. Agreed — with a shared function the static must
   be per-variant, which in Rust means either a `match` returning `&'static
   OnceLock` or two statics anyway. Two functions is less machinery. Keep it.
2. **A single static with a bitmask** (bit 0 = read fired, bit 1 = write fired).
   Fewer statics, but reintroduces the coupling §1.2 is removing. Reject.
3. **Eager emission of the description notice at startup**, beside
   `main.rs:2021`. Discussed in §1.7 — better grouping, worse precision.
   Rejected in favour of the plan's lazy approach, but a legitimate option if
   Felix prefers the caveats in one block.
4. **Gate the notice behind a verbosity flag** (A4's fallback). Rejected: the
   users who need it are the ones who would not enable it, and once per process
   is already the cheapest possible disclosure.
5. **Reject `--rename` + `--output-format ubam`** the way `--clump_only` already
   does (`cli.rs:846`). Out of scope for #406 and it removes working
   invocations, but it is the more honest resolution of §1.5 than a notice.
   Belongs in the follow-up issue, for Felix to decide.

---

## 7. Other items from the brief

**CHANGELOG placement — verified and correct.** `### Unreleased` → `#### Bug
fixes` exists at `CHANGELOG.md:4-6` and already carries #391 and #381 entries. A
message that actively misinformed users belongs there, not under a docs or
feature heading. The plan's caution about the `#### Fixes` / `#### Bug fixes`
duplication hazard (#401) is well placed — the heading in the file today is
`#### Bug fixes`.

Two additions to the bullet: mention the cross-direction suppression fix (§1.2)
and the read-side "for FASTQ output" correction (§1.3), since both change
user-visible output beyond the reversed text.

**Implementation outline ordering — sound, one gap.** Steps 1→2→3 (split the
static, then retext each side) are correctly ordered: splitting first means each
text edit lands against its own guard and cannot be attributed to the wrong
direction. Step 5 fixes the doc comments at `bam.rs:796-800` and `:865` — both
confirmed to describe the coercion without naming a direction. It should also
cover **`bam.rs:697-706`**, `parse_name_and_data`'s own doc comment, which says
the space tail is "a descriptive FASTQ-only annotation (SPACE — discarded)"
without noting that a tag tail *after* a space is discarded with it (§1.6b) —
that is the comment that would have prevented the false clause in the drafted
message.

**Test plan (outline item 6) — agree, with one addition.** Not pinning message
prose is right and consistent with what both #400 reviewers concluded
independently. The existing mechanics tests do cover the branch distinction:
`parse_name_splits_on_space_and_discards_annotation` (`:1249`, and it uses the
real BS-seq header) and `parse_tab_separator_still_treated_as_tag_tail`
(`:1262`). The plan's "add a case only if the notice's trigger condition turns
out not to be exactly the SPACE branch" is the right trigger — and per §1.6(a)
**it has already turned out that way**, so the empty-description case is now a
required test, not a conditional one.

---

## 8. Action items

### Critical

1. **Retract the "nearly unreachable" claim and re-scope the static split as a
   bug fix** (§1.2). Mixed FASTQ + uBAM single-end input reaches both directions
   in one process; demonstrated in both argument orders. In FASTQ-first order the
   user gets the wrong direction *and* loses the correct message. Update
   §Context, A-in-context, and the CHANGELOG bullet.
2. **Fix the read-side text too — it is not "already correct"** (§1.3). On
   uBAM→uBAM it says "coerced to N for FASTQ output" when the output is a BAM.
   Change "for FASTQ output" → "on read". Update Behavior 2, A3, and outline
   item 3, which all currently assert the text needs no change.
3. **Drop or qualify the notice's final clause** (§1.6b). "Tab-separated
   TAG:TYPE:VALUE fields are preserved" is false when a description precedes the
   tab tail — verified: `@readX desc\tCB:Z:ATCG` produces a record with no `CB`
   tag. Shipping it would reproduce the exact defect this plan fixes.

### Important

4. **Guard the notice on a non-empty description** (§1.6a). A trailing space
   survives the reader (`fastq.rs:412` strips only `\n`/`\r`), so keying on "the
   SPACE branch was taken" as Behavior 6 specifies fires on `@read1 ` with
   nothing dropped. Key on `!rest.is_empty()` and add the unit test to
   validation 5.
5. **File a follow-up issue for the `--rename` loss** (§1.5) and reference it
   from the plan. `--rename` is not in cli.rs §3.4a's uBAM exclusions, and
   `append_to_id` splices before a TAB but not before a SPACE, so `:clip5:ACG`
   is silently destroyed on uBAM output for any space-bearing header — measured.
   Do not fix it here; do not let it go unrecorded.
6. **Use `NOTE:` (uppercase) for the description notice** (§3.1). `Note:` appears
   zero times in `src/`; `NOTE:` appears 5 times, including the sibling
   FASTQ→uBAM line at `main.rs:2021`. Resolves the plan's first open question.
7. **Keep `1:N:0:INDEX` in the message** (§3.2). Resolves the second open
   question, and note the plan's premise for doubting it is wrong: five fixtures
   have exactly that shape, and `docs/…/clump-only.md:51` already uses
   `1:N:0:INDEX` as the canonical example of this behaviour.
8. **Replace validation 6's manual byte compare with
   `cargo test --release --test integration_ubam_out`** (§4.2) — that is the
   repo's `@PG`-ignoring convention (`assert_ubam_eq`, `integration_ubam_out.rs:108`)
   and it compares against committed goldens rather than a self-produced build.
9. **Record validation 3's recipe and its two harness traps in the plan**
   (§4.1): distinct stems (or the collision pre-flight aborts the run) and do
   not grep for a term that appears in the filenames.

### Optional

10. Restate A5's justification (§2.2): an exhaustive `stderr.is_empty()`
    assertion does exist (`integration_no_args_help.rs:42`) but is on the
    bare-invocation path; CI's one count assertion (`ci.yml:457`) counts a
    substring. The conclusion holds; the stated reason does not.
11. Soften A4's "nearly every real run" (§1.7): 9 of 20 FASTQ fixtures have
    whitespace-free headers and will never fire it. The judgement to fire it is
    still right.
12. Extend outline item 5 to `parse_name_and_data`'s own doc comment
    (`bam.rs:697-706`), which omits the space-swallows-tag-tail behaviour (§7).
13. Consider documenting the drop outside `clump-only.md` (§1.4/§7): the drop
    applies to all three uBAM output families, but only the clump-only page
    mentions it. `guide/outputs.md` is the natural second home.
14. Decide whether the notice should be per-process or per-input-file (§1.8) —
    the adjacent `main.rs:2021` NOTE is per-file, so the two will differ in
    cardinality on multi-input runs.

---

## 9. Bottom line

Proceed. The scope Felix set is the right scope, both defects are real, and I
confirmed the headline one by running the unpatched binary. But the plan should
not be implemented exactly as written: it under-claims the shared static (a
reachable suppression bug, not a tidiness issue), it wrongly exempts the
read-side message from correction, and its drafted notice ships a sentence that
is false on the very input that triggers it. Items 1–3 change what gets written;
4–9 are cheap and make the difference between a fix and a complete fix.
