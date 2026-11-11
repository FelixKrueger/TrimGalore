# Plan Review B — #406 uBAM normalization notices

**Plan:** `plans/08092026_ubam-normalization-notices/PLAN.md`
**Reviewer:** B (independent; Reviewer A ran separately)
**Base verified against:** `dev` @ `90dda11`
**Method:** code reading plus eight live runs of the already-built `target/release/trim_galore`
(2.3.0, binary contains the current message string) against purpose-built fixtures in a
scratch dir. Every claim marked **verified** below was reproduced, not inferred.

Scope decision (direction fix + description notice, uppercasing stays silent) is taken as
given and not relitigated.

---

## 1. Logic review

### 1.1 The two call sites and their direction mapping — verified, with one correction

| Site | Function | Direction | Plan's verdict | Mine |
|---|---|---|---|---|
| `src/bam.rs:823` | `validate_and_normalize_seq_for_write` | FASTQ → uBAM (write) | reversed | **Confirmed reversed** |
| `src/bam.rs:944` | `bam_record_to_fastq` | uBAM → FASTQ (read) | "already correct" | **Only correct for uBAM→FASTQ; false for uBAM→uBAM** |

Both call sites confirmed at the stated lines, both routed into the single
`emit_iupac_warning_once` at `src/bam.rs:995-1003`.

Write side reproduced (FASTQ containing one `R`, `--output-format ubam`):

```
WARNING: input uBAM contains IUPAC degenerate bases (...); coerced to N for FASTQ output.
```

Input was FASTQ, output was BAM. Both halves backwards, exactly as the plan says.

The correction is on the read side. `bam_record_to_fastq` is called from all three
`BamReader` paths (`src/bam.rs:261`, `:380`, `:451`) — the threaded reader, the interleaved
de-interleaver, and the sync reader — **regardless of output format**. So on a
uBAM-in → uBAM-out run, which is a first-class supported path, the read-side message fires
and tells the user the coercion happened "for FASTQ output" when the output is BAM.
Reproduced:

```
$ trim_galore --output-format ubam -o outB iupac.bam
WARNING: input uBAM contains IUPAC degenerate bases (...); coerced to N for FASTQ output.
NOTE: uBAM-in → uBAM-out: missing-qual sentinel (BAM 0xFF) does NOT round-trip ...
```

The FASTQ intermediate is real, but it is an implementation detail the user never asked
for and never sees. This is the same defect class #406 was filed about — a message naming a
direction the run does not have — so Behavior 2 ("text unchanged in substance") would ship
half the fix. One-clause change: `coerced to N in the output` (or `in the trimmed output`).

One thing the plan does not check, and which the new write-side wording depends on: can the
write side ever fire with "input FASTQ" when the record actually came from a BAM? No.
`bam_record_to_fastq` is the only BAM→`FastqRecord` conversion in the crate and it coerces
before returning, so a BAM-derived record reaching `validate_and_normalize_seq_for_write`
carries only `ACGTN`. The new wording's precondition holds unconditionally, including in
mixed-input runs (each input is processed as its own file). Worth stating in the plan,
because the wording is only safe *because* of that.

### 1.2 The shared `OnceLock` — the plan's "nearly unreachable" argument is wrong

The plan argues cross-direction suppression is nearly unreachable because on uBAM→uBAM the
read side coerces first. That reasoning is sound for uBAM→uBAM, but it misses the case that
actually matters: **a single-end run with both a FASTQ and a BAM input.** Nothing rejects
that shape. `reject_bam_format_mismatch_in_pair` (`src/main.rs:347-362`) is gated on
`cli.paired && cli.input.len() > 1`; the single-end uBAM loop (`src/main.rs:1979-1998`)
iterates inputs and re-detects format per file with no cross-input agreement check.

Reproduced, both orderings:

```
$ trim_galore --output-format ubam -o out_mixed  iupac_fq.fastq iupac.bam   # 1 IUPAC message
$ trim_galore --output-format ubam -o out_rev    iupac.bam iupac_fq.fastq   # 1 IUPAC message
```

Both coercions demonstrably happened — `samtools view` shows `readA ACGTN…` (write-side
`R`→`N`) in one output and `bamread1 ACGTN…` (read-side `Y`→`N`) in the other — yet only one
of the two warnings was emitted, and which one depends on argument order. This is not
exotic: the message's own examples (PacBio HiFi, ONT Dorado) are exactly the data where a
user might hold one run as FASTQ and another as uBAM and trim both in one invocation.

The plan's **action** (two independent statics) is correct and unchanged. Only the
rationale needs replacing, and it matters: "nearly unreachable" invites a later reader to
undo the split as unnecessary. It also upgrades validation 3 from a hedge into a concrete
expected-fail control — see §4.

### 1.3 Scope of the description-drop notice — verified

- `parse_name_and_data` (`src/bam.rs:707`) has exactly one production caller,
  `BamWriter::write_record` (`src/bam.rs:589`). Every other reference is a unit test in the
  same file. Verified by repo-wide grep.
- `write_record` serves all three uBAM output families, so the notice is not clump-specific:
  trim (`src/trimmer.rs:678`, `:790`, `:791`), hardtrim (`src/specialty.rs:170`, `:223`),
  clump-only (`src/clump_only.rs:694`, `:706`, `:707`). Verified.
- The drop is structurally required (SAM QNAME forbids whitespace). Agreed.

**`--rename` interacts, and badly.** `--rename` is rejected with `--clump_only`
(`src/cli.rs:846`) but **not** with `--output-format ubam` — the §3.4a block
(`src/cli.rs:582-627`) does not list it, and `cli.rename` is passed straight into
`hardtrim5_to_bam` / `hardtrim3_to_bam` (`src/main.rs:443`, `:474`) and into
`TrimConfig.rename` (`src/main.rs:1083`). `FastqRecord::append_to_id`
(`src/fastq.rs:173-181`) inserts the suffix before the first **TAB** but has no notion of a
SPACE, so on a space-bearing header the `:clip5:<seq>` annotation lands *after* the
description — inside the part `parse_name_and_data` discards.

Reproduced:

```
$ trim_galore --hardtrim5 20 --rename --output-format ubam -o outR  iupac_fq.fastq
readA   ACGTNACGTACGTACGTACG          # header was "@readA A00690:224:…/1" → annotation GONE

$ trim_galore --hardtrim5 20 --rename --output-format ubam -o outR2 nospace.fastq
readF:clip5:ACGTACGTACGTACGTACGT      # header was "@readF" → annotation SURVIVES
```

So on headers shaped like this repo's own fixtures, `--rename --output-format ubam` silently
discards the annotation the user explicitly asked for, at exit 0. The planned notice *will*
fire on those runs (something was dropped), but its drafted text calls the loss a "FASTQ
header description", which under-reports what happened. Two consequences: the notice's
wording must cover "anything after the first space, including `--rename` annotations", and
the underlying interaction deserves its own issue (fixing it properly means either
rejecting `--rename` with `--output-format ubam` or making the uBAM path insert before the
first whitespace — the latter cannot be done in `append_to_id` globally without changing
FASTQ output, which is Perl-parity-sensitive).

To the question "could the notice fire on IDs TrimGalore itself constructed": yes, via
`--rename`, and that is precisely the case where firing is *most* warranted. No other
ID-mutating flag reaches this path — `--clock`, `--implicon`, `--demux`, `--passthrough` are
all rejected with `--output-format ubam` (`src/cli.rs:598-627`).

### 1.4 Trigger condition — Behavior 6 misses two cases

The splitting logic is `stripped.find(['\t', ' '])`, first-whitespace-wins
(`src/bam.rs:709-722`). The plan's three enumerated cases are right as far as they go, but:

**(a) Empty description after the space → false alarm.** `FastqReader::next_record` trims
only `\n` and `\r` (`src/fastq.rs:446`), so a header written as `@name ` (trailing space —
common padding from some writers) keeps the space in `record.id`, takes the SPACE branch
with `rest == ""`, and would fire a notice announcing a loss that did not occur. Gate on
`!rest.trim().is_empty()`. A trailing **TAB** is already safe (`src/bam.rs:730-732` skips
empty fields; tested at `src/bam.rs:1392`).

**(b) `@name desc<TAB>TAG:TYPE:VALUE` → the notice fires *and* real aux tags are dropped,
while the drafted message promises the opposite.** Because the SPACE is found first, `tail`
becomes `None` and the whole tab tail is discarded. Reproduced:

```
@readC\tCB:Z:AAAACCCC            → out: readC … CB:Z:AAAACCCC      (tag preserved)
@readD desc1:N:0:XX\tCB:Z:AAAACCCC → out: readD … (no CB)          (tag silently lost)
```

The drafted final sentence — "Tab-separated TAG:TYPE:VALUE fields are preserved" — is false
on exactly the record that triggered the notice. Reachability is narrow (TrimGalore's own
`BamReader` never produces a space before the tab, since SAM QNAME cannot contain one, so
this needs a hand-made or third-party FASTQ), but a blanket promise printed on nearly every
run should not be falsifiable at all. Reword to something conditional, e.g. "aux tags in a
tab-separated tail immediately after the name (as written by `samtools fastq -T`) are
preserved".

Multiple spaces (`@name  desc`) and no-whitespace (`@name`) behave as the plan says.

### 1.5 Message wording — recommendations on both open questions

**Prefix: use `NOTE:`, and this is decidable from the codebase rather than taste.** Grep of
`src/*.rs` string literals: `WARNING:` × 10, `NOTE:` × 5, `Note:` × 0 (one legacy
`Warning:`). More to the point, the three existing uBAM-path startup messages
(`src/main.rs:385`, `:396`, `:409`) plus `src/main.rs:2021` ("NOTE: input is FASTQ; output
uBAM records will have no aux fields.") are the same register as this one — format-mandated,
non-actionable, fires on ordinary runs — and they all use `NOTE:`. `WARNING:` in this
codebase is reserved for "you may have made a mistake" (`--preserve-tags` with no BAM input,
the IUPAC coercion). The description drop is the former. The plan's drafted `Note:` matches
nothing in the tree; `NOTE:` also groups visually with the two or three NOTEs already
printed on the same run.

**Naming `1:N:0:INDEX`:** defensible either way, and the plan's own worry is well-founded.
Of the fixtures reachable with uBAM output, the dropped tails are instrument identifiers
(`BS-seq_10K_R1` → `A00690:224:H7GC5DRX2:1:2101:1081:1016/1`; `SRR24766921_RRBS_R1` →
`A00686:91:HTHYKDMXX:1:1101:19696:1000/1`), not index fields. The fixtures that *do* carry
a literal `1:N:0:` tail — `clock_10K_R1` (`1:N:0: CGATGTTT`) and `demux_test`
(`1:N:0:CTTNTA`) — belong to modes rejected with `--output-format ubam`. Counterweight:
`docs/src/content/docs/modes/clump-only.md:51` and the `--clump_only` help text
(`src/cli.rs:186`) already use `1:N:0:INDEX` as *the* example for this exact behaviour, so
keeping it is consistent with what is already published.

My recommendation is to sidestep the choice: **echo the actual dropped text of the first
affected record**, truncated to ~60 chars. It is concrete, needs no representative example,
and is self-explanatory for whichever header shape the user actually has. It costs one
`&str` argument and one format inside the once-guard. If you would rather not, keep
`1:N:0:INDEX` for consistency with the two places it is already documented.

### 1.6 Implementation outline — ordering and one inaccuracy

Ordering is sound (bam.rs → comments → tests → CHANGELOG → fmt/clippy/test). Two notes:

- Validation 8 (capture the pre-patch reversed text) belongs *before* step 1, not after
  step 8. It is already done — the transcript in §1.1 is that capture, taken from the
  committed `dev` binary — so paste it into the PR rather than re-deriving it.
- Item 5 says "the three doc comments … (`bam.rs:796-800`, `:865`)": two locations named for
  three comments, and `:865` documents `bam_record_to_fastq`, i.e. the **read** side, not the
  write side. The write-side mentions are the inline comment at `src/bam.rs:605-606` and the
  doc comment at `:796-800`. Name all three precisely so none is missed.

### 1.7 CHANGELOG placement — correct as planned

`### Unreleased` → `#### Bug fixes` exists at `CHANGELOG.md:4-6` with two live bullets.
A stderr message that actively misinformed users is a bug fix, and a previously-silent loss
becoming disclosed is fixing an under-reporting defect, so one bullet covering both belongs
there. Do **not** open a second `#### Changes` heading under Unreleased — that is the #401
duplication hazard the plan already flags (`#### Fixes` at `:317` and `:529`, `#### Changes`
at `:143` are all historical sections). The bullet should say what is now disclosed
(everything after the first space, including `--rename` annotations) rather than only
"description".

---

## 2. Assumptions

| # | Plan's assumption | Verdict |
|---|---|---|
| A1 | Nothing pins either message | **Verified.** `"IUPAC degenerate"` / `"coerced to N"` appear only in `src/bam.rs`. Adjacent prose exists in `docs/.../clump-only.md:52` and `src/cli.rs:187` but neither quotes the message. |
| A2 | Description drop is write-side only, one production caller | **Verified** (§1.3). |
| A3 | Read-side text already accurate | **False for uBAM→uBAM** (§1.1). This is the one assumption that changes the work. |
| A4 | Near-universal firing is acceptable | **Endorsed** — see below. |
| A5 | No test asserts stderr line counts; adding a line cannot break a test | **Verified** with one boundary worth knowing. |

**A4 — my judgement: yes, fire it, once per run, at `NOTE:` level.** Verified that both
uBAM-reachable fixtures carry descriptions, and that the retained half is the SRA accession
while the dropped half is the original instrument identifier
(`@SRR24827378.1 A00690:224:H7GC5DRX2:1:2101:1081:1016/1`). Three reasons this is not noise:
(1) the loss is real and the more informative half is what goes; (2) once per process is
already the floor — one line on a run that prints a version banner, three or four other
NOTEs, a per-file Trimming/Output pair and a multi-line summary is not where the signal gets
lost; (3) on the trim and hardtrim uBAM paths this notice is currently the *only* disclosure
anywhere — `--output-format`'s own help (`src/cli.rs:254-258`) and
`docs/src/content/docs/guide/outputs.md`'s uBAM section say nothing about it; only the
`--clump_only` surface documents it. Suppressing it would leave two of three output paths
silent about a real loss. The `--rename` case (§1.3) makes the argument stronger, not weaker:
there the dropped text is something the user explicitly asked to be there.

**A5 — verified, with the boundary.** No test or CI step asserts stderr equality or line
counts. The negative assertions that exist (`tests/integration_adapter2.rs:331`, `:430`,
`:434`; `tests/integration_paired_format_guard.rs:94`, `:98`, `:108`, `:135`, `:139`, `:182`,
`:274`, `:300`; `tests/integration_non_restartable_input.rs:142`, `:169`, `:240`, `:321`;
`tests/integration_output_collision.rs:234`, `:604`, `:634`, `:670`;
`tests/integration_clump_only_ubam.rs:442`) all match on substrings our text does not
contain. The only count-style checks are `tests/integration_adapter2.rs:581`
("Scanning for poly-G content"), `.github/workflows/ci.yml:457`
(`grep -c 'Auto-detecting adapter type' == 2`) and `ci.yml:916` (`grep -c "^@PG"` on a
header dump) — none touch our lines. **The one boundary:**
`tests/integration_no_args_help.rs:42` asserts `stderr.is_empty()`. It exercises the no-args
help path, which never reaches BAM code, so it is safe — but it is the single place where a
stray stderr line would fail, so keep the emitters behind their branches and do not hoist
either to startup.

One consequence the plan does not mention: ~16 unit tests in `src/bam.rs` call
`parse_name_and_data`, and at least one (`:1249`) takes the SPACE branch, so `cargo test`
output will gain the notice once. Harmless, but expect it.

---

## 3. Efficiency

Agreed: nil. Both notices sit behind an already-computed branch, add no scan and no
allocation, and after the first hit cost one atomic load plus a branch per record. Two
corrections of detail, neither affecting the conclusion:

- `OnceLock::get_or_init` is an **Acquire** load, not relaxed as the plan states.
- "No new allocation" stops being true if you adopt the echo-the-dropped-text suggestion
  (§1.5) — but the format runs once, inside the guard, so the per-record cost is unchanged.

The write-side call is already guarded by `if iupac_seen` per record, so the common case
(no IUPAC anywhere) does not touch the static at all. The description notice, by contrast,
will call its emitter once per space-bearing record, i.e. essentially every record on real
data — still a single acquire load, which is invisible next to BAM encoding, but it is the
one of the two that actually executes per record.

---

## 4. Validation sufficiency

**Validation 3 is constructible — resolved.** Two independent ways, and the plan should stop
hedging:

1. *End-to-end, and it doubles as an expected-fail control:*
   `trim_galore --output-format ubam <iupac.fastq> <iupac.bam>` (single-end, mixed formats,
   accepted — §1.2). Pre-patch it prints **one** IUPAC warning though both sides coerced;
   post-patch it must print **two**. I built the BAM with `samtools import -0`, and noodles
   encodes IUPAC faithfully (`=ACMGRSVTWYHKDBN`,
   `noodles-bam-0.89.0/src/record/codec/encoder/sequence.rs:92`), so an in-repo fixture can
   equally be written from a test with the `noodles` dep already in `Cargo.toml` — no
   external samtools needed for CI.
2. *Structurally,* two separate `static`s make the property local. If you also want it
   asserted: it must be asserted from the **integration** binary (subprocess), never a unit
   test — `OnceLock` is process-global and `validate_seq_coerces_iupac_to_n`
   (`src/bam.rs:1746`) already consumes the write-side guard in the lib test binary, so any
   in-process "fires exactly once" assertion is order-dependent and will flake.

**Validation 6 is the right check, and the repo already has the harness — use it by name.**
`assert_ubam_eq` (`tests/integration_ubam_out.rs:108`) compares header-minus-`@PG` plus
per-record `(name, flags, seq, qual, sorted aux)` tuples against the committed goldens
`test_files/ubam_out_se_REFERENCE.bam` and `ubam_out_pe_REFERENCE.bam`
(`:154`, `:171`). Those run under plain `cargo test`, so validation 6 is a subset of
validation 7 — no manual before/after byte diff needed. CI adds a second net via the
`--clump_only` uBAM record-multiset diffs (`.github/workflows/ci.yml:860-890`).

**Validation 5 — name the fixtures; they already exist.**
`test_files/BS-seq_10K_I1.fastq.gz` has ids with no whitespace at all (`@SRR24827378.1`) →
must not notify. `test_files/ubam_test_with_tags.bam` with
`--preserve-tags CB,UB --output-format ubam` exercises the TAB branch end-to-end (same shape
as the CI step at `ci.yml:1117-1134`) → must not notify, and the `CB`/`UB` tags must still
land in the output.

**Gaps — the silent-wrong-result modes the current table would not catch:**

- **Read-side text on uBAM→uBAM** (§1.1). Validation 2 says "uBAM fixture → FASTQ output",
  which is the one direction where the existing text is fine. Add: uBAM fixture with IUPAC
  → `--output-format ubam`, assert the message does not say "for FASTQ output".
- **Empty description** (§1.4a). A header ending in a bare space must not notify.
- **`@name desc<TAB>TAG`** (§1.4b). Notice fires, tag is dropped — assert the message does
  not claim tags are preserved.
- **`--rename --output-format ubam` on a space-bearing header** (§1.3). Notice must fire;
  document that the annotation is lost.
- **Multi-input runs.** Once-per-process means a five-file run where only file 4 has
  descriptions still notifies (at the moment it happens), but a run where files 1 and 4 both
  have them notifies once. That is intended; state it so it is not read later as a bug.

---

## 5. Alternatives

- **Emit from `write_record`, not from `parse_name_and_data`.** Have the parser return
  `description_dropped: bool` (or a small struct) and let `BamWriter::write_record` decide to
  emit. Keeps stderr out of a pure parse function, and — the real win — makes the trigger
  condition unit-testable without touching global state, which is what the plan's item 6 is
  currently dodging around. Cost: ~16 destructuring sites in `src/bam.rs`'s own test module.
  Worth it if you adopt the two new edge cases in §1.4, because those are much cheaper to
  assert on a returned `bool` than through subprocess stderr.
- **One emitter with a direction enum** instead of two functions: the plan already picks two
  functions and gives the right reason. Agreed, no change.
- **Per-writer rather than per-process notice.** Would give one notice per output file in a
  multi-input run. Arguably more useful, but it is a scope change from "once per run" and I
  would not take it now.
- **Gate the notice behind a verbosity flag** (the plan's own fallback if A4 is rejected):
  I would not. There is no verbosity flag today, inventing one for a single line is worse
  than the line, and per §2 the trim and hardtrim uBAM paths would then have no disclosure at
  all.

---

## 6. Action items

### Critical

1. **Fix the read-side text too.** `bam.rs:995-1003`'s "coerced to N **for FASTQ output**" is
   false on uBAM-in → uBAM-out, a supported path (reproduced, §1.1). Behavior 2 and A3 assert
   the opposite. Make it output-neutral ("coerced to N in the output"). Without this, #406
   ships with the same defect it was filed for, just on the other side.

### Important

2. **Drop the "Tab-separated TAG:TYPE:VALUE fields are preserved" promise, or make it
   conditional.** Verified false when the header is `@name desc<TAB>TAG:TYPE:VALUE` — the tag
   is silently discarded on the very record that triggers the notice (§1.4b).
3. **Widen the notice's wording to cover `--rename`, and file the interaction separately.**
   `--rename --output-format ubam` silently discards the `:clip5:`/`:clip3:` annotation
   whenever the header has a space — reproduced on this repo's own fixture shape (§1.3).
   `--rename` is not in the §3.4a rejection list. Fixing it (reject, or insert before the
   first whitespace on the uBAM path only) is out of #406's scope; the notice not misdescribing
   it is not.
4. **Gate the notice on a non-empty description.** `@name ` (trailing space) takes the SPACE
   branch because `FastqReader` trims only `\n`/`\r` (`src/fastq.rs:446`) — a notice there
   announces a loss that did not happen (§1.4a).
5. **Use `NOTE:`, not `Note:`.** Zero `Note:` in the tree; the four sibling uBAM-path
   messages in the same register all use `NOTE:` (§1.5).
6. **Replace the A4 "nearly unreachable" rationale.** Cross-direction suppression is
   reachable and reproduced via a mixed FASTQ+BAM single-end run (§1.2). Keep the split;
   correct the reason, and promote validation 3 to a real expected-fail control.
7. **Add a docs/`--help` line, or say explicitly why not.** The plan has no docs step, yet the
   description drop / uppercasing / IUPAC coercion are documented **only** on the
   `--clump_only` surface (`src/cli.rs:181-193`, `docs/.../clump-only.md:51-52`).
   `--output-format`'s help (`src/cli.rs:254-258`) says nothing, and
   `docs/.../guide/outputs.md`'s uBAM section claims "The output records themselves round-trip
   losslessly" — which is not true for FASTQ input (description, case, IUPAC) nor, strictly,
   for uBAM input with IUPAC. One sentence in each place.
8. **Name the existing byte-identity harness in validation 6** (`assert_ubam_eq` +
   `test_files/ubam_out_{se,pe}_REFERENCE.bam`) instead of a manual diff; it already runs under
   `cargo test` (§4).

### Optional

9. Echo the actual dropped text (truncated ~60 chars) instead of inventing an example —
   resolves the `1:N:0:INDEX` question without guessing (§1.5). If you keep an example,
   `1:N:0:INDEX` at least matches what is already published.
10. Consider returning `description_dropped: bool` from `parse_name_and_data` and emitting
    from `write_record` — testability plus no I/O in a parser (§5). ~16 test call sites.
11. Never assert once-ness in a unit test; `OnceLock` is process-global and
    `validate_seq_coerces_iupac_to_n` already consumes the write-side guard in that binary
    (§4).
12. Fix the two detail slips: "relaxed atomic load" → Acquire (§3); item 5's "three doc
    comments (`:796-800`, `:865`)" names two locations, one of which is the read side — the
    write-side mentions are `src/bam.rs:605-606` and `:796-800` (§1.6).
13. Move validation 8 to step 0 — already captured in §1.1, paste it into the PR (§1.6).
14. Name the negative fixtures in validation 5: `test_files/BS-seq_10K_I1.fastq.gz`
    (no whitespace) and `test_files/ubam_test_with_tags.bam` + `--preserve-tags CB,UB`
    (TAB branch) (§4).
15. State in Behavior that once-per-process is per **run**, not per input file, so a
    later reader does not file it as a bug (§4).
