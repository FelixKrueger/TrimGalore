# PLAN — #415: reject whitespace in a uBAM QNAME instead of silently dropping tags

**Issue:** [#415](https://github.com/FelixKrueger/TrimGalore/issues/415)
**Branch base:** `dev` @ `1fc00f6` (line numbers re-derived at that tree — the issue's
`bam.rs:709` predates #417, #420 and #422)
**Type:** behaviour change (a currently-succeeding input class starts failing), read-side only
**Decisions taken before drafting (Felix, 2026-08-11):** hard error always; reject all ASCII
whitespace.
**Decisions taken at manual review (Felix, 2026-08-11):** O1 **upheld** — the check stays
read-side on every output format, for the reason in §2.3. `\v` (VT, `0x0B`) is **accepted**;
the reject set for QNAME is `u8::is_ascii_whitespace`, and its principle is splitter-protection
rather than spec-enforcement (§2.7).
**Decisions taken after dual plan review (Felix, 2026-08-11):** the partial-output residue
(A7, now known false) is **tested and documented here, cured in a separate issue** covering all
five per-record bails. The `Z`-tag-value corruption found by Reviewer A is **folded in** — see
§2.8 and the second predicate in §4.

---

## 1. Goal

A BAM whose QNAME contains whitespace violates the SAM spec, but samtools builds one
happily and we accept it — then silently drop the read-name tail **and** every aux tag the
user explicitly asked to preserve. Turn that into a refusal at the point of conversion.

This also makes #408's assumption A5 enforced rather than inherited. #408 gated `--rename`
on input format, arguing that uBAM input is safe *by SAM spec* because a BAM read name
cannot carry a description. That reasoning is correct about the spec and unenforced in the
code; this is the one input shape where the admitted arm is not actually safe.

---

## 2. Context

### 2.1 The mechanism, confirmed by running it

`bam_record_to_fastq` (`src/bam.rs:871`) builds the FASTQ id as `@{QNAME}` and then appends
`\tTAG:TYPE:VALUE` per preserved tag (`:898-922`). No whitespace validation on QNAME.

`parse_name_and_data` (`src/bam.rs:707`) splits that id on the **first** of `\t` or `' '`
(`:709`). A space inside the QNAME therefore wins the race against the tab that introduces
the tag tail, so the whole tail — tags included — is classified as a discardable FASTQ
description.

Reproduced at `1fc00f6` with a samtools-built fixture (`QNAME = "name with space"`, one
`CB:Z:AAACCC` tag), and independently re-reproduced byte-for-byte by both plan reviewers:

| Run | Result |
|---|---|
| `--output-format ubam --preserve-tags CB` | **exit 0**; output QNAME `name`; `CB` tag **gone**; notice fires echoing `"with space\tCB:Z:AAACCC"` |
| `--preserve-tags CB` (FASTQ out) | **exit 0**; header `@name with space\tCB:Z:AAACCC` — name and tag text both physically present. A tab-aware consumer recovers the tag but reconstructs a spec-invalid QNAME; a whitespace-splitting one folds `with space\tCB:Z:AAACCC` into a comment. Nothing we were given is destroyed on this path — see §2.3 |
| `--clump_only --output-format ubam` | **exit 0**, tag dropped silently — and `--clump_only` advertises byte-identical preservation, which is the exact ground on which `cli.rs:853-855` refuses `--rename`. A documented-lossless mode silently breaking its own contract (Reviewer B) |

### 2.2 Four symptoms, one root cause

The issue names the first. Probing found the rest. Ordered by severity — symptom 2 is what
carries §2.3, and symptom 4 carries nothing:

1. **Aux tags dropped on the uBAM-output path** — the reported defect.
2. **LF or CR in a QNAME corrupts the output FASTQ structurally.** `FastqRecord::write_to`
   (`src/fastq.rs:78-96`) writes the id verbatim and then `\n`, so a newline inside the name
   emits a five-line record and desynchronises **every record after it** in the file. Silent,
   unbounded, and it lands on the FASTQ-output path — the one where the reported defect does
   not reach. **Demonstrated end-to-end** by both reviewers with a byte-level fixture: three
   records in, **13 lines out, exit 0** (§5 step 3). Not in the issue.
3. **The same corruption via a preserved `Z` tag value.** `append_tag_type_and_value`
   (`src/bam.rs:1058-1061`) writes `Z:{s}` verbatim into the same id with no validation, so an
   LF inside a preserved tag value breaks framing identically — same function, same id, same
   `write_to`. A tab inside one splits the tag tail into a bogus field and surfaces as
   `tag name must be 2 chars` (`bam.rs:748`), which says nothing about the input. Found by
   Reviewer A; **in scope** per §2.8.
4. **`fastq::read_id_prefix` truncates the paired sync key at the first ASCII whitespace**
   (`src/fastq.rs:208`). Two templates agreeing up to their first whitespace would pass the
   R1/R2 sync check that exists to catch desynchronised mates. **Stated with its weakness:**
   on the likely real input — a converter that put the whole Illumina header line into QNAME —
   truncating at the first space is precisely what `read_id_prefix` is *for*, the prefix is
   the template name, and the check still works. The failing case is contrived, so this
   symptom is not load-bearing anywhere in the plan.

### 2.3 The fix is read-side, so it also fires on the FASTQ-output path

**Settled at review (Felix, 2026-08-11): upheld.** `bam_record_to_fastq` is the BAM→FASTQ
conversion point, so validating there covers every output format. The scope is wider than the
issue implies, and the class of runs that stop working is not exotic: a uBAM whose QNAMEs
carry whitespace, converted to FASTQ. Its likeliest real origin is a converter that put the
whole Illumina header line into QNAME, and those files pass through correctly today.

**Narrowing by output format is not available.** The obvious narrower fix — refuse only when
`--output-format ubam` — was ruled out, not preferred-against: symptom 2 lands on the
FASTQ-output path. An LF in a QNAME desynchronises the whole output file, which is a larger
loss than the lost tag the issue reports. A gate on output format would close a lost-tag hole
and leave a lost-file hole open. So the real choice was between refusing on every path and
splitting the rule by character class (reject LF/CR/FF always, space/tab only on uBAM output).

**Why every path, and not the class split:** input validity must not be discovered by changing
an output flag. Under the class split, the same BAM is accepted until the day the user adds
`--output-format ubam`, and the malformation was there the whole time. Space and tab on the
FASTQ path are also not harm-free so much as harm-deferred — we hand the next tool a name
that samtools' `-T` convention round-trips back into a spec-invalid QNAME, and a bwa `-C`-style
comment passthrough turns `with space\tCB:Z:AAACCC` into invalid SAM tags. The class split
also costs two messages and a rule that reads "whitespace is illegal, except space and tab,
unless you are writing BAM".

**The cost, priced and accepted:** runs that lose nothing today will fail, and the remediation
is renaming reads at whatever scale the file is. Accepted on the judgement that header-in-QNAME
uBAMs are a rare conversion mistake rather than a common shape. **If that judgement is wrong,
the class split is the fallback** — it is a strictly smaller change than this plan, not a
rewrite of it.

**The priced cost is narrower than v2 thought.** Both reviewers measured that `\v`, `@` and
non-ASCII QNAMEs *already* fail on the uBAM-output path today, at the noodles writer (§2.5). So
the only class where this change causes a genuinely **new** failure is space / tab / CR / FF on
the **FASTQ-output** path. That is exactly the class this section prices, and nothing more.

**Three arguments deliberately retired at review, so no reviewer leans on them:**

| Retired argument | Why it does not hold up |
|---|---|
| "The header is malformed for any whitespace-splitting consumer" | `@name\tCB:Z:AAACCC` is *already* comment-territory for bwa/bowtie2 — tab is whitespace to them. A tab-aware consumer parsing `@name with space\tCB:Z:AAACCC` recovers the tag intact. The real harm narrows to `-C`-style comment passthrough emitting invalid SAM, which is downstream and conditional |
| "Symptom 4 applies regardless of output format" | True but near-empty: on the likely input, prefix truncation is correct behaviour (§2.2) |
| "Gating on output format re-conditionalises the invariant #408 wanted closed" | Output-format-conditional refusal is the **house style**: `cli.rs::Cli::validate` §3.4a already rejects `--demux`, `--clock`, `--clumpify` and `--retain_unpaired` specifically for uBAM output. #408 wanted its `--rename` refusal off an unenforced assumption; it never asked for an unconditional QNAME rule |

### 2.4 The chokepoint is single; the entry points are four

**`bam_record_to_fastq` has exactly three production callers.** v2 mislabelled two of them
(inherited from the session handoff), and §5 step 4's test design was wrong as a result.
Corrected against the code:

| Call site | What it actually is | Error wrapper | Reached by |
|---|---|---|---|
| `src/bam.rs:261` | **paired interleaved de-interleaver** (`FREAD1`/`FREAD2` split, `tx_r1`/`tx_r2`, `MAX_SLACK`) | `:268`, `"{path}: BAM record {n}: {msg}"` | `--paired` on one interleaved BAM |
| `src/bam.rs:380` | **threaded single-stream reader** (`open_threaded_inner`, batches into one `tx`) | `:394` | SE with `--cores ≥ 2` |
| `src/bam.rs:451` | sync reader (Direct) | — | default SE run, and the sanity check |

`parse_name_and_data` has exactly one caller, `BamWriter::write_record` (`:589`).

**Four entry points feed those three call sites, and three of them run before the trimming
loop.** This is what makes the validation design in §9 non-obvious:

1. `main.rs:243` — `sanity_check_any(&cli.input[0])`, before all dispatch: `BamReader::open` +
   one `next_record()` (`main.rs:38-46`) → record 1 through the `:451` reader.
2. `adapter.rs:115` — adapter auto-detection, up to `MAX_SCAN_READS = 1_000_000`
   (`adapter.rs:84`), skipped when an adapter or preset is given.
3. `adapter.rs:268` — `detect_poly_g`, also ≤ 1 M records, skipped by `--no_poly_g` / `--poly_g`.
4. Only then the trimming loop.

Consequence, measured by both reviewers: a bad **record 1** dies at the sanity check with no
writer open; a bad record 2 under `-a` alone dies in the poly-G scan; a bad record 2 under
`-a --no_poly_g` reaches the loop — and only that configuration exposes the partial output in
§8's A7. On a real file over 1 M reads the pre-scans stop shielding anything, so the
`-a --no_poly_g` arm is the faithful one, not the exotic one.

### 2.5 Our own *record encoder* cannot emit this — and its rule is far stricter than whitespace

> **Corrected in round 3.** This section originally said "we cannot produce one", which is only
> true of the record encoder. `bgzf::Writer` has no validation, so Rust can write these bytes —
> and now does, in `examples/mk_bam_fixture.rs`. The looser reading sent v3 to a Python
> generator and survived five reviews.

`noodles`' record encoder rejects the write outright. The actual predicate
(`noodles-bam/src/record/codec/encoder/name.rs`) is:

```rust
(MIN_LENGTH..=MAX_LENGTH).contains(&buf.len())        // 1..=254
    && buf != MISSING                                 // != b"*"
    && buf.iter().all(|&b| b.is_ascii_graphic() && b != b'@')
```

`is_ascii_graphic()` is `0x21..=0x7E`, so the write side rejects **everything except printable
ASCII minus `@`**: all whitespace, VT, every non-ASCII byte, `@`, and names over 254 bytes. Two
consequences:

- Trim Galore can never *produce* such a BAM, so the defect is strictly read-side, and neither
  fixture can be built with our own code (§5 step 3).
- **The three classes §2.7 accepts are not "harmless" on the uBAM-output path** — they fail
  today, exit 1, with `Error: failed to write BAM record / Caused by: invalid input parameter`
  and a partial BAM left behind. The message never mentions the read name. Both reviewers
  measured this; the edge-case table in §3 now says so. Pre-existing and out of scope here,
  but worth its own issue (§5 step 8).

### 2.6 Fixing the primary defect makes the secondary defect unreachable

The issue's second half: the #406 notice says "**FASTQ** header text after the first space"
even on a uBAM-input run (`emit_description_dropped_once`, `src/bam.rs:1029-1043`).

After this change, a space can no longer reach `parse_name_and_data` from BAM input. The id
built by `bam_record_to_fastq` is `@{QNAME}` plus tab-prefixed tags; with whitespace QNAMEs
rejected, the only remaining space would have to sit inside a `Z` tag value, and the tab
that introduces that tag always precedes it — `find(['\t', ' '])` returns the earliest
match, so the tab wins and the tail parses as tags. The notice therefore fires only on
FASTQ input, where "FASTQ header text" is accurate.

Both reviewers checked this argument and both found it structurally sound rather than merely
probable: `str::find` over a `char` array returns the earliest match of *either* needle, and
`append_to_id` splices the `:clip5:` suffix *before* the first tab (`fastq.rs:176-177`), so
nothing can reorder them.

**Decision: resolve that half by deleting the possibility, not by rewording.** A test pins
that a whitespace-QNAME BAM now errors, and a comment records why the wording is
unconditionally correct. #406 split its IUPAC messages per direction because both
directions were reachable; that is not true here, and adding words to defend an unreachable
state is how a message becomes untrustworthy.

### 2.7 The QNAME reject set comes from our splitters, not from the SAM spec

v1 framed the check as spec-enforcement and then listed exceptions to the spec (`@`,
non-ASCII). That is backwards, and it produced a self-contradiction: §6 proposed
`u8::is_ascii_whitespace` while the edge-case table demanded `\v` be rejected. **Rust's
`is_ascii_whitespace` follows WHATWG Infra and excludes VT** — verified with `rustc`, and
independently re-verified by both reviewers:

```
space 0x20  is_ascii_whitespace=true      VT 0x0B  is_ascii_whitespace=false
  tab 0x09  is_ascii_whitespace=true      FF 0x0C  is_ascii_whitespace=true
   LF 0x0A  is_ascii_whitespace=true      CR 0x0D  is_ascii_whitespace=true
```

An implementer following v1's §6 would have written the §3 test and watched it fail.

The rule that resolves it: **reject the bytes that our own code treats as a boundary.** v2
claimed three such sites and that their union was *exactly* `u8::is_ascii_whitespace`. Both
reviewers falsified that claim — there are five, and the union is **Unicode** whitespace:

| Site | Boundary set |
|---|---|
| `parse_name_and_data`, `src/bam.rs:709` — `find(['\t', ' '])` | space, tab |
| `fastq::read_id_prefix`, `src/fastq.rs:210` — `split_ascii_whitespace` | space, `\t`, `\n`, `\r`, `\f` |
| FASTQ line framing, `src/fastq.rs:78-96` — id written verbatim, then `\n` | `\n` (and `\r`, which survives as a trailing byte on the id line) |
| **`FastqRecord::append_to_id`, `src/fastq.rs:175` — `self.id.trim_end()`** | **Unicode `White_Space`** — a strict superset: includes VT, NBSP (U+00A0), U+2028, U+3000 |
| `demux.rs:223`/`:240` — `record.id.replace(' ', "_")` | space only, so it does not extend the union — but it means demux silently *renames* such reads |

**So VT's acceptance is a judgement, not a derivation.** The fourth row strips a **trailing**
`\v` from the output read name — but only under `--rename` together with a clip flag that
actually clipped (`trimmer.rs:261-274`), and `--rename` stays admitted for BAM-in → FASTQ-out
after #408. Both reviewers demonstrated it end-to-end (`@readtrail:clip5:ACG`, the VT gone at
exit 0). That is the same "silently repairing a name changes read identity" harm the edge-case
table refuses elsewhere.

Accepting VT anyway is still right: rejecting it read-side would over-reject on exactly the
grounds this section gives, and trailing-VT-plus-`--rename` is vanishingly narrow. What changed
is the **status** of the acceptance — knowingly scoped collateral, not a proven-safe derivation.
V5 pins VT-acceptance, so this distinction has to be written down or the next reader will find
`trim_end` and conclude the plan was wrong rather than deliberate.

Spec-invalidity remains context in the error message, not the reason for the refusal.

### 2.8 The `Z`-tag-value predicate (Reviewer A's C5, folded in)

Symptom 3 corrupts the output through the same id and the same `write_to` as symptom 2, so
stopping at QNAME would ship a half-applied rule: `LF` in a read name refused, `LF` in a
preserved tag value still silently corrupting the file. A user cannot tell the two apart.

**The predicate is different, and deliberately so.** A space inside a `Z` value is legal
(`Z` is `[ !-~]*`) and common, so the QNAME predicate must not be reused. What breaks us:

| Byte in a `Z` value | Effect | Verdict |
|---|---|---|
| `\n`, `\r`, `\f` | breaks FASTQ line framing exactly as symptom 2 | reject |
| `\t` | splits the tag tail into a bogus field → `tag name must be 2 chars` (`bam.rs:748`) | reject |
| space | legal and meaningful; the introducing tab precedes it, so `parse_name_and_data` still parses the tail correctly | **accept** |

Only reachable with `--preserve-tags`, and only from a byte-level-constructed BAM — but that is
the same reasoning that failed for QNAME, and both reviewers built such a BAM in an afternoon.

---

## 3. Behavior

1. In `bam_record_to_fastq`, after the existing empty-name check (`src/bam.rs:894-896`) and
   before the id is built, scan the QNAME bytes with `reject_whitespace_in_qname`.
2. If any boundary byte is found, return an error naming the record, the offending name with
   whitespace escaped so a tab or newline is visible, and a remediation.
3. In the tag-appending loop (`:900-922`), scan each `Z` value with `reject_framing_in_tag_value`
   before it is appended (§2.8). Space stays legal.
4. Otherwise proceed exactly as today. No other behaviour changes.
5. Every BAM input path inherits both checks — single-end and paired, FASTQ and uBAM output,
   `--clump_only` and the specialty modes included — because all of them read through the three
   callers in §2.4.
6. Nothing on the FASTQ-input side is touched: a space in a FASTQ header remains a
   legitimate description and continues to be dropped with the existing notice.

### Edge cases

| Input | Behaviour |
|---|---|
| QNAME with a space | error |
| QNAME with a tab | error |
| QNAME with `\n`, `\r`, `\f` | error |
| QNAME with `\v` (VT, `0x0B`) | **accepted** — §2.7. Not in `u8::is_ascii_whitespace`. Known collateral: a *trailing* VT is stripped by `append_to_id`'s `trim_end` under `--rename` + a clip flag. Accepted knowingly |
| QNAME with non-ASCII UTF-8 | **accepted by this check.** On uBAM output the noodles writer still rejects it today (`invalid input parameter`, partial BAM) — pre-existing, §2.5 |
| QNAME containing `@` | **accepted by this check.** Same writer caveat as above |
| QNAME longer than 254 bytes | unchanged — accepted here, rejected by the writer on uBAM output with the same opaque message. Pre-existing, §5 step 8 |
| Leading/trailing whitespace in QNAME | error — no trimming; silently repairing a name changes read identity |
| Empty QNAME | unchanged, already an error (`:894`) |
| QNAME that is invalid UTF-8 | unchanged, already an error at `:897`. The new check runs **before** it, so a name that is both invalid UTF-8 *and* whitespace-bearing reports the whitespace — the more actionable of the two |
| A `Z` tag value containing a space | **accepted** — legal, and the introducing tab precedes it (§2.8) |
| A tag value containing `\n`, `\r`, `\f` or `\t` | error — §2.8. Applies to **every** aux type that emits text, not just `Z`: the scan runs over whatever the emitter appended |
| A `Z` tag value that is both invalid UTF-8 **and** framing-bearing | reports the **UTF-8** error, the opposite precedence to the QNAME row above. The tag scan runs after the emitter, which needs a `&str`; the QNAME scan runs on raw bytes before the conversion. Deliberate asymmetry, recorded so it is not read as an inconsistency |

---

## 4. Signature

```rust
/// Reject a BAM read name holding any byte our own splitters treat as a boundary.
///
/// The set is `u8::is_ascii_whitespace`: `parse_name_and_data` splits on space/tab,
/// `fastq::read_id_prefix` on ASCII whitespace, and a FASTQ id occupies one line.
fn reject_whitespace_in_qname(name: &[u8]) -> Result<()>;

/// Reject bytes that would break the FASTQ id this `Z` value is appended to.
///
/// Space is legal in a `Z` value; `\t` would forge a tag field and `\n`/`\r`/`\f` would
/// break the single-line id.
fn reject_framing_in_tag_value(tag: &str, value: &str) -> Result<()>;
```

Both kept as free functions rather than inlined so they can be unit-tested per character class
without constructing a BAM record.

The docstrings state the current rule and its reasons in two lines each, per the comment
convention: no spec citation (that belongs in the user-facing error, §5 step 1), no account of
the defect being fixed, and each invariant sits on the function that enforces it.

---

## 5. Implementation outline

1. **`src/bam.rs`** — add both predicates; call `reject_whitespace_in_qname` from
   `bam_record_to_fastq` immediately after the empty-name check, and
   `reject_framing_in_tag_value` inside the tag loop. QNAME error text:

   > `BAM read name contains whitespace: "name with space". Trim Galore reads a FASTQ header
   > as ending at the first whitespace, so any --preserve-tags aux tags after it are lost on
   > uBAM output, and a newline would break the FASTQ record structure. The SAM specification
   > forbids whitespace in QNAME. Rename the reads before trimming — the uBAM section of the
   > docs gives a samtools recipe for space-only names; a tab or newline in a name cannot be
   > expressed in SAM text at all.`

   Three review corrections are baked into that wording. v1 asserted aux-tag loss
   unconditionally, which is false on the FASTQ path; v2 fixed that but still said "a FASTQ
   header ends the read name at the first whitespace", which a user looking at their own
   FASTQ output (where the whole name is visibly present) reads as contradicted — so the rule
   is now attributed to **Trim Galore's reader**. And the tag clause is qualified with
   `--preserve-tags`, since without it there are no tags to lose.

   **The remediation claim is now specific, because Reviewer B falsified the general one.** v2
   said no clean samtools one-liner exists. This works, and was verified running (I reproduced
   it as well — tag intact, then round-trips through Trim Galore):

   ```bash
   samtools view -h in.bam \
     | awk 'BEGIN{FS=OFS="\t"} /^@/{print; next} {gsub(/ /,"_",$1); print}' \
     | samtools view -b -o fixed.bam -
   ```

   `FS="\t"` is the whole trick — a bare `awk` splits on spaces and shifts every field, which
   is presumably why v2 concluded none existed. It does **not** generalise: tab and LF names
   fail SAM round-trip with `Parse error at line 3`, because SAM text cannot represent either.
   The recipe therefore goes in the docs (step 6) with that caveat, and the error text points
   at it rather than dumping a wrapping pipeline into stderr.

   Escape the name for display. Use `escape_debug` on the lossy string, which covers `\t`,
   `\n`, `\r` **and** `\f` for free — v2 listed only `\t`/`\n`, and a CR-only name would
   otherwise print as an unremarkable string that overwrites itself in a terminal. Use
   `String::from_utf8_lossy`, not `?` on UTF-8 — an invalid name should not turn into a
   confusing encoding error, and the check runs before the existing UTF-8 conversion at `:897`.

   The callers already wrap with `"{}: BAM record {}"` (§2.4), so the user sees
   `Error: file.bam: BAM record 12345` above the message. On a large file that index is the
   practical handle; no extra work needed, but worth knowing before someone adds it.

2. **`src/bam.rs` unit tests** — QNAME: one per rejected byte (space, tab, `\n`, `\r`, `\f`),
   plus accepted cases (plain ASCII, non-ASCII UTF-8, `@`, and **`\v`** as the §2.7 pin, so a
   later hand-rolled predicate that "fixes" VT goes red). Tag values: `\n`/`\r`/`\f`/`\t`
   rejected, space accepted. Assert the error text contains the escaped name, so the escaping
   is pinned rather than assumed.

   **These stay unit tests.** Do not promote the accepted classes to integration cases: `\v`,
   `@` and non-ASCII all fail at the noodles writer on uBAM output (§2.5), so such a test would
   go red and the natural-looking fix — widening the reject set — would silently undo a
   maintainer decision.

3. **Fixtures — five, and a committed generator.** Three are samtools-built; two need a
   byte-level writer. Note the samtools-built *inputs* in `test_files/` today are three
   (`ubam_test.bam`, `ubam_test_with_tags.bam`, `ubam_paired_test.bam`); the two
   `ubam_out_*_REFERENCE.bam` files are Trim Galore outputs, not hand-built.

   ```bash
   printf '@HD\tVN:1.6\n' > ws.sam
   printf 'name with space\t4\t*\t0\t0\t*\t*\t0\t0\tACGT…\tIIII…\tCB:Z:AAACCC\n' >> ws.sam
   samtools view -b -o test_files/ubam_ws_qname.bam ws.sam
   ```

   | Fixture | Shape | Pins |
   |---|---|---|
   | `ubam_ws_qname.bam` | single record, space in QNAME | the sanity-check entry point (`main.rs:243`) |
   | `ubam_ws_qname_late.bam` | 3 records, **clean first**, space at record 2 | the trimming loop — the only arm that reaches a writer (§2.4) |
   | `ubam_ws_qname_paired.bam` | interleaved, **clean first pair**, whitespace on a later pair, flags `77`/`141` | the `:261` de-interleaver |
   | `ubam_lf_qname.bam` | byte-level, LF in QNAME, 3 records | §2.3's load-bearing fact, end-to-end |
   | `ubam_lf_tagvalue.bam` | byte-level, LF inside a `CB:Z` value | §2.8 |

   **v2's "a newline fixture is not constructible" was false, and both reviewers proved it by
   building one.** SAM *text* cannot express an LF in QNAME, but BAM's read-name field is
   NUL-terminated bytes, so a ~25-line stdlib Python generator writes the record directly and
   wraps it in a BGZF block. samtools 1.21 reads the result back at exit 0 with the raw LF
   intact. Commit the generator (`test_files/mk_bam_fixture.py`) and record it in
   `test_files/README.md` beside the existing samtools recipes — a fixture reproducible only
   from a script is fine here, that is what the README already documents.

   Two traps, both paid for already: **BGZF `BSIZE` is block length − 1**, and an off-by-one
   yields htslib's `CRC32 checksum mismatch`; and a *genuinely aligned* fixture record needs an
   `@SQ` line, because htslib silently forces `FUNMAP` when `RNAME` is `*`. A working generator
   and both LF fixtures are preserved at
   `plans/08112026_ubam-qname-whitespace/revA_artifacts/` — start from `mkbam.py` rather than
   rediscovering the record layout.

   CI needs no samtools and no Python: all five fixtures are committed.

4. **Integration tests** — in `tests/integration_ubam.rs`. **v2's design was wrong and would
   have passed for the wrong reason:** every case used a single-record fixture, so the refusal
   fired at `sanity_check_any` before any writer existed, V6's "nothing written" held
   vacuously, and V7 never reached the de-interleaver at all. Name the entry point in each test
   name so this cannot silently regress:

   - **Sanity-check arm.** `ubam_ws_qname.bam`, `--output-format ubam --preserve-tags CB` →
     non-zero exit, message names the whitespace, **nothing written**. This is a real path, just
     not the only one.
   - **Trimming-loop arm, uBAM out.** `ubam_ws_qname_late.bam` with `-a AGATCGGAAGAGC
     --no_poly_g` → non-zero exit **and a partial `*_trimmed.bam` present**. Assert the residue
     that actually occurs (§8 A7), not its absence.
   - **Trimming-loop arm, FASTQ out.** Same fixture and flags → refused (§2.3's widened scope,
     asserted deliberately so the decision is visible in a test rather than implied). Note BAM
     input yields plain `.fq`, **not** `.fq.gz` — the gzip rule mirrors the input extension —
     so do not glob `*.fq.gz`.
   - **`:380` arm.** Same fixture, `--cores 2`, FASTQ out → refused. Nothing else covers the
     threaded single-stream reader.
   - **`:261` arm.** `ubam_ws_qname_paired.bam`, `--paired -a … -a2 … --no_poly_g` → refused,
     with the de-interleaver's one-line wrapper shape (`{path}: BAM record {n}: {msg}`,
     `bam.rs:268`) rather than the sanity check's `Caused by:` shape. The producer runs on a
     background thread, so assert on process exit and stderr.
   - **LF end-to-end.** `ubam_lf_qname.bam`, FASTQ out → refused after the fix. The pre-fix
     corruption (13 lines for 3 records) is recorded in §9 V1b rather than committed as a test
     that must then be deleted.
   - **`Z`-tag-value arm.** `ubam_lf_tagvalue.bam` with `--preserve-tags CB` → refused (§2.8).
   - **Acceptance twin.** A clean uBAM still round-trips with its tag preserved, so a check that
     rejected everything could not pass.
   - **The unreachability pin for §2.6:** a FASTQ input with a genuine description still
     emits the notice, and its wording ("FASTQ header text") is now always accurate.

5. **`--clump_only` guard.** Not a formality: on the space fixture,
   `--clump_only --output-format ubam` drops the tag **silently at exit 0** today, while
   `--clump_only` advertises byte-identical preservation — the same contract `cli.rs:853-855`
   cites when refusing `--rename`. Add the case; it is the strongest single argument in §2.2.

6. **Docs.** One sentence in the uBAM **input** section — that is
   `docs/src/content/docs/quickstart.md:55` ("## Unaligned BAM (uBAM)"); `guide/outputs.md:54`
   is output-only — noting that read names containing whitespace are refused **on every output
   format**, with the reason. The output-format independence is the part a user will be
   surprised by, so it belongs in the sentence rather than being left implicit. Add the samtools
   remediation recipe from step 1 with its "space-only" caveat. No `--help` change: this is not
   a new flag.

7. **`CHANGELOG.md`.** Every one of the last six commits on `dev` touched it, and #398's
   `### Unreleased` entry carries an explicit "refuses some invocations that previously ran"
   note — this change is exactly that shape. **One sentence**, in the file's existing format.
   v2 omitted this from both §5 and §7.

8. **File the two out-of-scope siblings** so they are not lost:
   - **Partial output on any per-record bail.** Decided: cured separately, covering all five
     triggers (aligned, reverse-complemented, secondary, supplementary, and now whitespace).
     Unlink-on-error, or write-to-temp-and-rename. Cite `main.rs:402` — #408 already declined a
     mid-stream refusal *for this exact reason* — and #409 as the precedent for treating
     partial output as a data-integrity bug.
   - **The opaque noodles write error.** `\v` / `@` / non-ASCII / >254-byte QNAMEs fail on uBAM
     output with `invalid input parameter`, no read name, and a partial file. Same defect shape
     as #415 with a different trigger.

---

## 6. Efficiency

One pass over the read name per record — `name.iter().any(u8::is_ascii_whitespace)`, a few
bytes, no allocation on the happy path. The tag scan adds one pass over the preserved `Z`
values, only when `--preserve-tags` is given. The escaped-name formatting happens only on the
error path. Negligible against the existing per-record decode work.

The path form, not a closure, or clippy's `redundant_closure` fails the `-D warnings` build.

Worth recording because it explains §2.4 and §9: on a default run the QNAME predicate executes
**2–4 times per record** — sanity check (record 1), adapter detection (≤ 1 M), poly-G scan
(≤ 1 M), then the loop. Still free, and it is precisely why the refusal usually surfaces before
any writer opens, which is what makes the untested regime (past 1 M records, or under
`-a --no_poly_g`) the dangerous one.

---

## 7. Integration

- **Reads/writes:** `src/bam.rs` only, plus tests, five committed fixtures, a committed fixture
  generator, one docs sentence, and one `CHANGELOG.md` line.
- **Not touched:** `parse_name_and_data`, `append_to_id`, `fastq::read_id_prefix` (v1 called
  this `id_prefix` and cited `src/fastq.rs:184-199`, which is its doc comment — the function is
  at `:208`), and the FASTQ-input path. In particular the splitter keeps its current shape —
  the fix removes the input that defeats it rather than rewriting it.
- **Known residue, deliberately not cured here:** a refused run can leave a partial output
  file, and the partial uBAM is structurally valid and samtools-readable (§8 A7). Pre-existing
  for four other per-record bails; this change adds a fifth trigger. Tested and documented
  here, cured in the issue from §5 step 8.
- **Also uncured, deliberately:** a multi-input SE run whose *second* input holds the bad
  record leaves input 1's outputs complete and input 2's partial — `sanity_check_any` only
  inspects `cli.input[0]`. Same class as the above; one V row records it.
- **CI:** the `validation` job compares against Perl 0.6.11, which has no BAM support, so it
  is unaffected. `--fastqc` and reproducibility are untouched.
- **#422's tripwire:** unaffected — its arm cases all use clean fixtures, and a refused run
  produces no *new* candidate paths.
- **Unblocks:** #408's A5 becomes an enforced invariant. Worth a note on #408 when this
  lands, since its justification currently rests on the unenforced version.
- **Incidental benefit worth claiming:** an LF QNAME on the uBAM-output path currently fails
  with the opaque `invalid input parameter` plus a partial BAM. After this change it fails
  read-side with the name in the message.

---

## 8. Assumptions

| | Assumption | Status |
|---|---|---|
| A1 | `bam_record_to_fastq` is the only BAM→FASTQ conversion point — 3 callers, 4 entry points (§2.4) | **Verified** by both reviewers, independently, by grep over every `BamReader::open*` and `open_sync_reader` construction in `src/` |
| A2 | The QNAME reject set is `u8::is_ascii_whitespace`. `\v`, non-ASCII and `@` are accepted | **Decided, not derived.** v2 claimed the splitter union was *exactly* this set; both reviewers falsified that — `append_to_id`'s `trim_end` makes the true union Unicode whitespace (§2.7). VT acceptance is knowingly-scoped collateral |
| A3 | Rejecting is right even where nothing is currently lost (space/tab + FASTQ output, §2.3) | **Upheld at review**, Felix, 2026-08-11. Strengthened: the new-failure class is narrower than v2 thought, because `\v`/`@`/non-ASCII already fail at the writer |
| A3b | Header-in-QNAME uBAMs are a rare conversion mistake, so A3's over-rejection cost is small | **Judgement, now with one data point.** Both reviewers ran `samtools import` in every mode: QNAME is the first whitespace-delimited token only, so the canonical converter **cannot** produce this shape. Picard `FastqToSam` still unchecked |
| A4 | samtools can build the space fixtures and CI needs no samtools to consume them | **Verified**, samtools 1.21, by both reviewers |
| A5 | `noodles` rejects a whitespace QNAME on write, so we cannot produce one | **Verified but understated** — the real rule is printable-ASCII-minus-`@`, 1..=254 (§2.5) |
| A6 | After the fix the description notice is reachable only from FASTQ input, so its "FASTQ" wording is always accurate | **Verified** by both reviewers; structurally sound, not probabilistic. Pinned by a test |
| A7 | ~~A refused run writes nothing.~~ **FALSE.** A refused run leaves a partial output on every path when the offender is past the pre-scan window, and the partial uBAM is valid, BGZF-EOF-terminated and samtools-readable — nothing downstream flags it | **Falsified by both reviewers, independently, by running it.** Measured residues: 319–329 B BAM, 1-of-2-record FASTQ, 282–293 B paired BAM, and 0-byte files under `--cores 2`. Mechanism confirmed: `main.rs:2147-2148` opens the reader then `BamWriter::create`, which does `File::create` + `write_header` immediately (`bam.rs:565-570`) |
| A8 | Moving the check earlier would not fix A7 | **Verified** — a sanity-check-time peek sees only record 1 (`main.rs:38-46`). Only output cleanup fixes it, hence §5 step 8 |
| A9 | Space is legal in a `Z` tag value and must stay accepted, so §2.8 needs its own predicate | Fixed by the `Z` grammar `[ !-~]*` |

---

## 9. Validation

| | Verify | How | Expected |
|---|---|---|---|
| V1 | The defect exists before the fix | Run the space fixture on `dev` | exit 0, tag lost — observed, §2.1 |
| V1b | §2.3's load-bearing fact exists before the fix | `ubam_lf_qname.bam`, FASTQ out, on `dev` | exit 0, **13 lines for 3 records** — observed by both reviewers |
| V2 | Refused after, sanity-check arm | Integration test, single-record fixture, uBAM out | non-zero exit, message names the whitespace, nothing written |
| V3 | Refused after, trimming-loop arm | `ubam_ws_qname_late.bam`, `-a … --no_poly_g`, uBAM out | non-zero exit **and a partial `*_trimmed.bam` present** (A7) |
| V4 | A3's widened scope is real, not implied | Same fixture and flags, FASTQ out | refused. Do not glob `*.fq.gz` — BAM input yields plain `.fq` |
| V5 | The `:380` threaded reader is covered | Same fixture, `--cores 2` | refused |
| V6 | The `:261` de-interleaver is covered | `ubam_ws_qname_paired.bam`, `--paired -a … -a2 … --no_poly_g` | refused, with the `:268` wrapper shape — **not** the sanity check's |
| V7 | LF in QNAME, end-to-end | `ubam_lf_qname.bam`, FASTQ out | refused |
| V8 | §2.8's tag-value scan | `ubam_lf_tagvalue.bam`, `--preserve-tags CB` | refused |
| V9 | Refusal is not over-broad | Clean-fixture acceptance twin asserting the tag round-trips | passes |
| V10 | Every QNAME boundary byte | Unit tests: space, `\t`, `\n`, `\r`, `\f` | all rejected |
| V11 | Accepted classes stay accepted | **Unit tests only** (§5 step 2): `\v`, non-ASCII UTF-8, `@`, plain ASCII | accepted — `\v` is the §2.7 pin |
| V12 | The tag predicate is not the QNAME predicate | Unit tests: space in a `Z` value accepted; `\t`/`\n`/`\r`/`\f` rejected | as stated |
| V13 | `--clump_only` inherits the refusal | One case on the space fixture | refused — and it is lossless-by-contract, §5 step 5 |
| V14 | Multi-input partial-run shape is recorded | Two SE inputs, second one bad | input 1 complete, input 2 partial; documented in §7, not cured |
| V15 | §2.6's unreachability | FASTQ-input description case | notice still fires, wording accurate |
| V16 | No regression | `cargo test`, `cargo test --release`, `cargo nextest run`, clippy `-D warnings`, `cargo fmt --check` | 643 + new tests pass, zero warnings |
| V17 | Compile check after the signature change | `cargo test --no-run` — `cargo build --release` does not build test targets | clean |

**How V2 and V3 assert output residue.** v2's phrasing — "the output directory holds only the
input" — is not expressible here: `tests/integration_ubam.rs::fresh_tmpdir` hands each test an
empty directory and the fixtures live in `test_files/`. Assert on the absence (V2) or presence
(V3) of `*_trimmed.*` / `*_val*` instead.

---

## 10. Questions or ambiguities

Two were answered before drafting: hard error always, and all ASCII whitespace.

**Resolved at manual review (Felix, 2026-08-11):**

- **O1 — §2.3. Upheld: the check stays read-side on every output format.** The reason changed
  even though the answer did not. v1 defended it with three arguments; review retired all three
  (§2.3's table) and replaced them with one: narrowing by output format was not actually
  available, because symptom 2 — an LF in a QNAME desynchronising the entire output FASTQ —
  lands on the FASTQ path, so an output-format gate would close a lost-tag hole and leave a
  lost-file hole. The remaining live alternative was a class split, priced and declined in §2.3.
- **O4 — the `\v` contradiction, found at review.** v1's §6 predicate accepted VT while §3 and
  V4 demanded it be rejected; an implementer following the plan would have written a failing
  test. Resolved by **accepting** `\v` and reframing the reject set as splitter-protection
  (§2.7) — which dual review then corrected again: the reframing is a decision, not a
  derivation, because a fourth site strips VT.

**Resolved after dual plan review (Felix, 2026-08-11):**

- **O2 — A7. Settled by running it: A7 is false.** Both reviewers measured a partial output on
  every path. Decision: **test and document the residue here; cure it in a separate issue**
  covering all five per-record bails (§5 step 8). Rationale: the residue is pre-existing and
  orthogonal, and a write-to-temp-then-rename change touches the uBAM writer's lifecycle,
  which deserves its own review and its own golden-fixture re-verification (`@PG`,
  `assert_ubam_eq`). What this plan must not do is keep asserting the invariant backwards.
- **O5 — Reviewer A's C5, the `Z`-tag-value corruption. Folded in** (§2.8). Stopping at QNAME
  would ship a half-applied rule for a corruption the user cannot distinguish.
- **O6 — sanity-check-time refusal, considered and declined.** It is the only design under
  which A7 would be true as written, but a *complete* pre-flight means a full extra pass over
  the input, and a partial one (first N records) recreates exactly the false-confidence regime
  in §2.4. Recorded so the next reader does not have to rediscover why the check is per-record.

**Open, non-blocking:**

- **O3.** Should `--rename`'s help text stop naming the uBAM restriction as spec-derived
  once it is code-enforced? Cosmetic, and #408 owns that wording.

---

## 11. Self-Review

**Logic.** The fix removes the input that defeats the splitter rather than changing the
splitter — deliberate, because `parse_name_and_data`'s first-whitespace rule is correct for
FASTQ input, where a space genuinely introduces a description. Rewriting it to prefer the
tab would make the tag tail survive but leave a spec-invalid name flowing onward into
`read_id_prefix` and, for uBAM output, into a writer that noodles will reject anyway.

**Edge cases addressed.** Newline in QNAME (structurally corrupts the output FASTQ, not in the
issue, and now the load-bearing reason for §2.3); the same corruption via a `Z` tag value
(§2.8); leading/trailing whitespace (no silent trimming — that changes read identity); `\v`,
non-ASCII and `@` accepted, with their write-side fate now stated rather than called harmless;
a space inside a `Z` value, which stays legal; a name that is both invalid UTF-8 and
whitespace-bearing, which reports the whitespace because the new check precedes `:897`.

**Verified rather than assumed.** The mechanism, both output paths' current behaviour, the
notice text, samtools' willingness to build the fixture, and noodles' refusal to write one
were established by running code at `1fc00f6`. Dual review then established by running code
that A7 is false, that an LF fixture is constructible, that the samtools one-liner exists for
space-only names, that `trim_end` is a fourth boundary site, and that `\v`/`@`/non-ASCII
already fail at the writer.

**What review caught that drafting did not.**

*Manual review (v1 → v2).*
1. The predicate and the test matrix contradicted each other — `is_ascii_whitespace` excludes
   VT, while §3 and V4 demanded VT be rejected. Both halves were written from the phrase "all
   ASCII whitespace" without checking what the standard library means by it.
2. Two of §2.3's three justifications did not survive scrutiny, and the fact that does settle
   it — LF corrupting the FASTQ framing — was present in v1 but never connected to the scope
   argument. The conclusion held; the reasoning under it was the weaker two-thirds.
3. The error message carried the defect it was fixing: it asserted aux-tag loss unconditionally
   while the check fired on a path with no tag loss.
4. `fastq::id_prefix` does not exist; the function is `read_id_prefix` at `src/fastq.rs:208`.

*Dual plan review (v2 → v3). Both reviewers ran code; the pattern is unmistakable.*
5. **A7 was false** — and the plan's own V6 would have confirmed it as true, because every
   fixture was single-record and therefore intercepted at the sanity check. A plan that
   proposes a validation which cannot fail has not validated anything.
6. **§2.7's derivation was incomplete** — `append_to_id`'s `trim_end` is a fourth boundary
   site using Unicode whitespace, so the union is not `is_ascii_whitespace` and VT's acceptance
   is a judgement. This was the strongest-feeling argument in v2, and it was wrong on a line
   one above code I had already read.
7. **"Not constructible" was false** — both reviewers built an LF-in-QNAME BAM from struct
   bytes and demonstrated the plan's load-bearing fact end-to-end.
8. **"No clean samtools one-liner" was false** — v2 went further and instructed reviewers to
   demand one be shown running. Reviewer B showed it running.
9. **The `:261`/`:380` labels were swapped**, inherited from the session handoff and not
   caught by Reviewer A either, which sent v2's paired test at the wrong call site.
10. A whole class of missing coverage: the `Z`-tag-value corruption, `CHANGELOG.md`, the
    `--clump_only` lossless-contract violation, the multi-input partial shape, and a docs
    target that does not exist under the name v2 gave it.

**The pattern worth carrying forward:** every false claim in v1 and v2 came from a
*derivation* — reading code and reasoning to a conclusion. Every claim that survived came from
*executing* something. The reviewers' technique is the transferable part: to test whether a
*proposed* error behaves correctly, find an *existing* error on the identical path and run
that. Both used the aligned-BAM bail inside `bam_record_to_fastq` to settle A7 before a line
of the new check existed.

**Remaining risks.**
1. **A3b is the live risk, and it is a judgement** — that header-in-QNAME uBAMs are rare. Now
   supported by one converter's measured behaviour (`samtools import` cannot produce the shape);
   Picard `FastqToSam` unchecked. §2.3 names the fallback if it is wrong.
2. **The partial-output residue ships with this change.** Known, measured, tested, documented,
   and cured elsewhere by decision. The risk is that the follow-up issue does not get filed —
   §5 step 8 exists so it does.
3. VT's acceptance is knowingly-scoped collateral, not a proven-safe derivation (§2.7). The
   distinction is written down because V5 pins the behaviour, and an undocumented pin invites a
   future reader to "fix" it.
4. The error message must survive an invalid-UTF-8 name without turning into an encoding
   error — hence `from_utf8_lossy` rather than `?`.

---

## Implementation notes (2026-08-12)

Implemented on `dev` at `1fc00f6`. Order followed §5, with the pre-fix evidence captured
first: the release binary was built at the clean tree and V1/V1b run **before** any source
change, so the tests that now assert refusal are known to have failed beforehand.

**Pre-fix evidence, measured at `1fc00f6` (binary stamped `1fc00f6`):**

| | Run | Result |
|---|---|---|
| V1 | space QNAME, uBAM out, `--preserve-tags CB` | exit 0, QNAME `name`, `CB` gone |
| V1b | LF QNAME, FASTQ out | exit 0, **13 lines for 3 records** |
| V1c | LF inside a `CB:Z` value, FASTQ out | exit 0, **13 lines for 3 records** |

V1c is new. Reviewer A found the `Z`-value defect by reading `append_tag_type_and_value`;
nobody had measured its end-to-end effect. It is identical to V1b, which is the argument for
having folded §2.8 in rather than deferring it.

**What was built.** `src/bam.rs`: `reject_whitespace_in_qname` and
`reject_framing_in_tag_value`, called from `bam_record_to_fastq` (after the empty-name check,
before the UTF-8 conversion) and from `append_tag_type_and_value`'s `Value::String` arm. Six
unit tests, ten integration tests (nine in `tests/integration_ubam.rs`, one in
`tests/integration_clump_only_ubam.rs`), six committed fixtures, a committed fixture generator,
recipes in `test_files/README.md`, a docs subsection, and one `CHANGELOG.md` entry.

**Deviations from the plan, and why.**

1. **`reject_framing_in_tag_value` takes only the value, not `(tag, value)` as §4 specified.**
   The tag name is supplied by caller context — `.with_context(|| format!("aux tag '{}'",
   tag_name))` at the call site — so the message names the tag without the predicate having to,
   and the check sits where the value has already been converted to `&str`. Calling it from the
   loop instead, as §3 item 3 said, would have meant duplicating the `from_utf8` conversion and
   producing two error messages for one failure.
2. **The tag predicate is derived, not enumerated.** `b != b' ' && b.is_ascii_whitespace()`
   rather than `matches!(b, b'\t' | b'\n' | b'\r' | b'\x0C')` — it is the QNAME set minus the one
   byte `Z`'s `[ !-~]*` grammar makes legal, so the two predicates cannot drift apart. **This
   corrects §2.8's own table**, which justified rejecting `\f` as "breaks FASTQ line framing".
   It does not; it is rejected because it is a boundary byte with no legitimate place in a tag
   value.
3. **Six fixtures, not five.** The first draft of the space-in-`Z`-value acceptance case shelled
   out to samtools at test time with a skip when absent — a silent skip, and against the
   committed-fixture convention. Replaced by `test_files/ubam_tagvalue_space.bam`.
4. **The acceptance twin was already there.** §5 step 4 asked for "a clean uBAM still
   round-trips with its tag preserved"; `preserve_tags_roundtrip_matches_golden` is exactly
   that, against a committed golden. Not duplicated — the full-suite run covers it. Checked that
   the golden fixture's tags (`CB:Z:ATCGATCG-1`, `UB:Z:GCTAGCTA`) contain no whitespace, so the
   new tag predicate cannot disturb it.
5. **`\v` (VT) is untested at the integration level, deliberately** (§5 step 2). It is accepted
   read-side but fails at the noodles writer on uBAM output, so an integration case would be red
   for an unrelated reason.

**Iteration log.**

- **#1** — the first unit test used `unwrap_or_else` on a `Result`, whose closure fires on
  `Err`; it would have panicked on precisely the bytes it was asserting are rejected. Rewritten
  to `is_err()`. Caught by re-reading, not by the compiler.
- **#2** — the space-in-`Z`-value case was rewritten from a samtools-at-test-time skip to a
  committed fixture (deviation 3).
- **#3** — `rustfmt` reflowed the `--clump_only` `.args([…])` call to one argument per line;
  applied with `cargo fmt --all`.
- **#4** — clippy reported clean in 6 s with a single `Compiling` line, which proves nothing
  about test-target coverage. Planted a `clippy::len_zero` violation in
  `tests/integration_ubam.rs`; clippy failed **exit 101** naming `test "integration_ubam"` at
  `:517`. Probe removed and clippy re-verified clean. The earlier pass was real.

**Found while building the fixtures, and not in the plan: the two checks have different
entry-point exposure.** `sanity_check_any` opens the reader with `BamReader::open(path)`
(`main.rs:39`) — **no tag list**. So `reject_framing_in_tag_value` can never fire at the
sanity-check entry point, because no tags are requested there. A malformed `Z` value therefore
*always* reaches the trimming loop, which means it always leaves a partial output — where a
malformed QNAME in record 1 is caught before any writer exists. The §2.4 entry-point analysis
holds for QNAME but is one entry point short for tag values. This makes the §5 step 8
partial-output issue slightly more urgent than §7 implies: for the tag-value class there is no
benign regime at all.

**Still open, by decision, not oversight.** The two follow-up issues in §5 step 8 are not
filed — filing them is outward-facing and was left to Felix. The partial-output residue ships
as documented (§7), asserted by V3 so it cannot change unnoticed.

### Round 2 — dual code review + coverage audit (2026-08-12)

`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`, `COVERAGE.md`. Both reviewers worked in isolated
worktrees and executed rather than derived; the coverage audit re-ran every claimed gate.
Felix elected the full sweep, so everything below is applied.

**The one real code defect — §2.8 was half-applied.** Both reviewers independently found that
`Value::Character` (the `A` type) writes its raw byte into the same FASTQ id through a different
match arm, unguarded, and both reproduced **exit 0, 13 lines for 3 records** — the V1b signature
this change exists to close. Reviewer A additionally measured that on uBAM output *Trim Galore
itself emits a spec-invalid `A` value*, which **falsifies §2.5's claim** that we can never
produce such a BAM: true for read names (noodles validates them) and false for `A` tag values.

Fixed by hoisting the scan to cover every arm, present and future:
`append_tag_type_and_value` captures `start = out.len()`, and the final statement is
`reject_framing_in_tag_value(&out[start..])`. **Reviewer B proposed the same idea with a bug** —
its sketch captured the offset *before* `id.push('\t')`, so the scanned slice began with a tab
and the predicate would have rejected every tag. Reviewer A's placement, inside the emitter, is
the correct one and is what shipped. Verified load-bearing: removing the hoisted call turns
`lf_in_character_tag_value_refused_end_to_end` **and** the `Z` case red, and `src/bam.rs`
restored md5-identical afterwards.

**Reviewer B's sharpest catch, which A missed:** deleting
`.with_context(|| format!("aux tag '{}'", tag_name))` left **all 659 tests green**, because the
test asserted `stderr.contains("aux tag")` and the predicate's *own* message already opens with
"aux tag value cannot pass…". Deviation 1's whole justification rested on an unpinned line. Now
asserts `"aux tag 'CB'"`.

**One contradiction between the reviewers, adjudicated by running it.** B said the docs'
`FS="\t"` caveat named the wrong mechanism; A said it was correct under gawk and to leave it.
**B is right** — bare `awk` output is byte-identical to the input. `$1` stops at the first
space, so the substitution matches nothing and `$0` is never rebuilt; A's premise that gawk
always rebuilds after `gsub` does not hold when zero substitutions occur. The docs now describe
a silent no-op, which is the quieter and worse failure.

**Both reviewers corrected the flakiness mechanism I asserted.** V3's non-emptiness does not
come from the header write — both `bam::io::Writer` and `bgzf::Writer` buffer, so nothing is on
disk when `create()` returns. It comes from `noodles-bgzf::Writer`'s `Drop` running
`try_finish()`. Not flaky (25 and 105 runs), but the conclusion was right for the wrong reason.
V3 now asserts the **BGZF EOF marker**, which pins the property that actually makes the residue
dangerous: it is indistinguishable from a complete BAM.

**Everything else applied.** The docs sentence no longer claims FASTQ output loses nothing (false
for LF/CR/FF — the very fact §2.3 rests on); `test_files/README.md` drops a `had no targets in
header` warning samtools 1.21 never emits and replaces a probe that *could not fail* with a
record count that does; the QNAME docstring and the VT test comment present the set as chosen
rather than derived (§2.7); the tag-value message borrows the description notice's 60-char cap,
since `Z` values are unbounded; all seven samtools/byte-level fixtures were rebuilt with
`--no-PG` so the README recipes reproduce the committed bytes exactly and no scratchpad path
survives in a committed artefact; the CHANGELOG is one sentence and says "whitespace other than
a space"; and the space-bearing `Z` value is now asserted on the uBAM→uBAM path, where the tail
has to be re-split around the space.

**Two fixtures added, taking the count to eight:** `ubam_lf_atag.bam` (newline in an `A` value)
and `ubam_atag_ok.bam` (printable `A` values). The second exists because my first acceptance test
for the `A` arm asserted `!success` — it depended on the refusal, so mutation turned it red along
with the refusal tests. An acceptance twin must stay green under mutation; that one didn't, which
made it a refusal test wearing the wrong name.

**Coverage audit gap closed:** the §2.6 comment now sits above `emit_description_dropped_once`,
recording why "FASTQ" is unconditional. It was the only plan requirement that was neither met
nor documented as a deviation.

### Round 3 — the fixture generator moved to Rust (2026-08-12)

Felix asked why a Python script belonged in this repo. It didn't, and the reason it was there
was **a misreading that this plan introduced and that nobody caught in five reviews**.

§2.5 said "`noodles` **rejects** a whitespace QNAME on write, so we cannot produce one." That is
true of noodles' *record encoder* (`write_alignment_record`, which validates against
`[!-?A-~]{1,254}`) and false of `bgzf::Writer`, which is a plain `Write` sink with no
validation. Everyone — the plan, both plan reviewers, both code reviewers, the coverage
auditor — read "noodles refuses" as "Rust cannot produce this file", and Reviewer B stated the
conclusion explicitly: *"a stdlib generator is the only route."* It never was. `src/format.rs:526`
has been writing a BGZF file from Rust in a test since the uBAM work landed.

So the generator is now `examples/mk_bam_fixture.rs`: the records are packed by hand, exactly as
the Python did, and only the BGZF framing comes from noodles. What this buys:

- **One language.** The script would have been the first Python the *test corpus* depended on for
  reproducibility — `src/`, `tests/` and `test_files/` had none.
- **No `python3`, no `bgzip`** for regeneration; the Rust toolchain is already required.
- **Reproducibility is stronger, not weaker.** All four fixtures regenerate **byte-identically**
  from their documented `cargo run --example` command; the Python route depended on bgzip's
  compression settings matching.
- **`cargo clippy --all-targets` lints it.** Proven by planting a `len_zero` violation: exit 101
  naming `example "mk_bam_fixture"`.
- An `examples/` target is not a `[[bin]]`, so CLAUDE.md's single-binary invariant is untouched.

**Provenance independence is preserved.** The one argument for an external generator is that a
fixture built by the crate under test could mask a bug in that crate's reader. Two things answer
it: samtools reads all four back with the expected records (independent validation of the bytes),
and the packing here deliberately avoids noodles' record encoder, so it shares no code with the
path being tested. Verified before the swap that the Rust and Python outputs carry **identical
stored headers and identical records** — they differ only in deflate byte layout.

**Also removed:** a stray `test_files/ubam_lf_qname.bam.bai` that one of the code reviewers left
in the main tree at 12:17 while probing. My earlier "the reviewers modified nothing" was checked
with `git status` filtered to modified files, which cannot show a new untracked one.

---

## Revision History

- **v1** (2026-08-11) — first draft. Mechanism, both output paths, and the fixture
  constraints established by running code rather than reading it.
- **v2** (2026-08-11) — manual review. **O1 upheld, its reasoning replaced:** two of three
  justifications retired, and the LF-corrupts-FASTQ-framing fact (v1's symptom 3, verified at
  `src/fastq.rs:78-96`) promoted to carry the scope decision, which also removed the
  output-format gate from the option space. **`\v` now accepted** — new §2.7 derives the reject
  set from the three splitters rather than the SAM spec, resolving a contradiction between v1's
  §6 predicate and its own test matrix. Error message rewritten for the widened scope. A3b added
  as the plan's one unmeasured judgement. Symptom order in §2.2 changed to severity, with
  symptom 3's weakness stated. `id_prefix` → `read_id_prefix` throughout.
- **v3** (2026-08-11) — dual plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`), both
  reviewers running code in isolated worktrees. **A7 falsified:** a refused run leaves a valid,
  samtools-readable partial output; tested and documented here, cured in a separate issue by
  decision. **§9 redesigned** — v2's cases were all intercepted at `sanity_check_any`, so V6
  and V7 were false passes; fixtures now put the offender past the pre-scan window and each
  case names its entry point. **§2.7 corrected** — a fourth boundary site (`append_to_id`'s
  `trim_end`, Unicode whitespace) makes VT's acceptance a judgement rather than a derivation.
  **§2.8 added** — the `Z`-tag-value corruption, folded into scope. **§2.4 corrected** — the
  `:261`/`:380` labels were swapped, and the four entry points are now named. **§2.5 corrected**
  — noodles' real rule is printable-ASCII-minus-`@`, so the three "harmless" classes already
  fail on uBAM output. "Not constructible" and "no samtools one-liner" both retracted, each
  falsified by a reviewer building the thing. Added: `CHANGELOG.md`, the `--clump_only`
  lossless-contract case, the multi-input shape, `escape_debug`, and the real docs target.
