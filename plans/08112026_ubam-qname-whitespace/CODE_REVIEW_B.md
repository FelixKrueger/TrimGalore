# CODE REVIEW B — #415, whitespace in a uBAM QNAME

**Reviewer:** B (independent, fresh context, no shared state with Reviewer A)
**Target:** uncommitted working tree on `dev`, base `1fc00f6`
**Method:** isolated worktree at `$TMPDIR/codereview-B-415` (detached `1fc00f6` + `git apply` of the
working diff + the six untracked fixtures copied in). Every claim below was **executed**, not derived;
where a probe was confounded I say so and give the corrected measurement.

---

## Summary

The fix is well-placed, well-tested and correct as far as it goes. `bam_record_to_fastq` really is
the single chokepoint — I verified the three production call sites and the four entry points
independently, and then verified empirically that **every** mode inherits the refusal
(`--hardtrim5/3`, `--clock`, `--implicon`, `--clumpify --cores 2`,
`--clump_only --output-format ubam`, `--demux`, `--cores 2/4`, `--fastqc`, `--paired` interleaved,
multi-input SE). Mutation testing kills 5 of 6 mutants. Full suite in the worktree: **659 passed, 0
failed**; clippy `-D warnings` green **and proven able to fail** (planted `len_zero` → exit 101);
`cargo fmt --check` clean. The `mk_bam_fixture.py` recipe reproduces both LF fixtures **byte-identically**.
The documented samtools/awk remediation works verbatim.

One finding is substantive: **the rule is still half-applied.** The tag guard sits in
`append_tag_type_and_value`'s `Value::String` arm, but `Value::Character` (the `A` type) writes its
raw byte into the same FASTQ id with no guard — and an `A` value holding `\n` reproduces the plan's
own load-bearing corruption signature exactly (**exit 0, 13 lines for 3 records**). §2.8's argument
for folding in the `Z` case transfers to `A` word for word.

Four documentation/README statements are inaccurate, two of them in ways that invert the reasoning
they are meant to convey.

Nothing here blocks the *approach*. H1 should be fixed before merge; M1–M4 are cheap.

### Answers to the nine questions asked

| # | Verdict |
|---|---|
| 1 | Reject set is right; message accurate on all four paths; `escape_debug` sane; **no** lifetime hazard. One docstring inaccuracy → **M5** |
| 2 | Derivation correct for `Z`; rejecting `\f` is **justified**, not over-rejection. But applied to too narrow a scope → **H1** |
| 3 | 3 call sites / 4 entry points **confirmed independently**; no BAM input path bypasses the QNAME check. The **tag** check is bypassed by the `A` arm → **H1** |
| 4 | `with_context` breaks no existing message, test or doc. But it is itself **untested** → **M1** |
| 5 | **Not flaky** (25 + 40 + 40 runs). Your stated mechanism is wrong though — see H/M-free note in §5 below |
| 6 | 5 of 6 mutants killed; **M3 (the `with_context`) survives** the whole 659-test suite |
| 7 | Fixtures minimal and correct. LF pair reproduces **byte-identically**. The four samtools ones do **not** → **L3** |
| 8 | Recipe verified working; `FS` caveat mechanism wrong (**M4**); README warning claim false (**M3**); one docs sentence inverts §2.3 (**M2**) |
| 9 | **Nothing to flag.** New comments are 1–3 lines; the files' pre-existing norm already runs to 3/4/6/10/20-line blocks. No "previously"/"no longer"/"before this change" anywhere |

---

## Issues by area

### HIGH

#### H1 — Logic: the identical corruption survives through the `A` (Character) aux-tag arm

`append_tag_type_and_value` (`src/bam.rs:1085`) has three arms that emit *text* into the FASTQ id.
Only one is guarded:

```rust
Value::Character(c) => write!(out, "A:{}", *c as char).unwrap(),   // src/bam.rs:1088 — UNGUARDED
…
Value::String(bs) => {
    let s = std::str::from_utf8(bs).context("BAM string tag is not valid UTF-8")?;
    reject_framing_in_tag_value(s)?;                                // src/bam.rs:1098 — guarded
```

`Value::Character` writes the BAM's raw `A`-field byte verbatim. Measured with a byte-level BAM
(3 records, `XS:A:\n` on record 2), `--preserve-tags XS -a … --no_poly_g`:

| Run | Result |
|---|---|
| `XS:A:\n`, FASTQ out | **exit 0, 13 lines for 3 records** — byte-for-byte the plan's V1b / §2.3 load-bearing signature |
| `XS:A:\n`, uBAM out | **exit 0**, the raw LF round-trips into the output BAM's `A` value |
| `XS:A:\t`, uBAM out | exit 1, `malformed FastqRecord.id: "@readB\tXS:A:\t"` + a **partial 1-record BAM** — the opaque-message shape §2.8 calls a defect when it happens for `Z` |
| `XS:A:\t`, FASTQ out | exit 0, the tab forges a field in the id |

The plan's §2.8 rationale applies verbatim with `A` substituted for `Z`: *"stopping at QNAME would
ship a half-applied rule: `LF` in a read name refused, `LF` in a preserved tag value still silently
corrupting the file. A user cannot tell the two apart."* And the reachability argument is identical —
`A` tags are common in real BAMs (`XS:A:+`, `ts:A:+`), and a framing byte in one needs the same
byte-level construction the `Z` fixture needed. §2.8 explicitly refuses that as an excuse: *"that is
the same reasoning that failed for QNAME."*

`Value::Float` and the six integer arms are safe — I verified `f32`'s `Display` yields `NaN`/`inf`/
`-inf`, no whitespace. `Hex` and `Array` already bail. So `Character` is the **only** gap.

**Recommended change.** Don't add a second call inside the `Character` arm — that repeats the drift
risk deviation 2 was written to avoid. Instead move the scan to where §3 item 3 originally put it,
the loop in `bam_record_to_fastq`, so it covers whatever the emitter wrote regardless of arm:

```rust
// src/bam.rs:953-957, sketch
let tail_start = id.len();
id.push('\t');
id.push_str(tag_name);
id.push(':');
append_tag_type_and_value(&mut id, &value)
    .with_context(|| format!("aux tag '{}'", tag_name))?;
reject_framing_in_tag_value(&id[tail_start..])
    .with_context(|| format!("aux tag '{}'", tag_name))?;
```

That covers `A`, every future arm, and the `TAG:` prefix itself, at one pass over a few bytes. It also
keeps deviation 1's benefit (the tag name comes from caller context) and preserves the space
exemption, which `Z` needs.

Note that this trades away one thing: with the scan after the emitter, a `Z` value that is *both*
invalid UTF-8 and framing-bearing keeps reporting the UTF-8 error first (see **L2** — that is already
today's behaviour, so nothing regresses).

Add: a unit test on the `A` arm mirroring `read_side_emitter_rejects_newline_in_string_tag`, and an
`A`-type LF fixture. The `mk_bam_fixture.py` generator needs ~4 lines for an `a_tag` helper — I built
one to produce the measurements above, so the fixture cost is small.

**Not in scope but worth knowing:** `XS:A:<space>` also round-trips at exit 0 into a spec-invalid `A`
value (`A` is `[!-~]`, space excluded). I would **not** distinguish space for `A` — a framing-only
rule is sufficient and keeps one predicate. Flagged only so it is not mistaken for a regression
introduced by the fix above.

---

### MEDIUM

#### M1 — Errors/tests: the `with_context` tag name is untested; mutant M3 survives all 659 tests

Removing `.with_context(|| format!("aux tag '{}'", tag_name))` from `src/bam.rs:957` leaves the
**entire suite green** — 419 lib + 16 clump_only_ubam + 17 integration_ubam, 0 failures.

The reason: `lf_in_tag_value_refused_end_to_end` (`tests/integration_ubam.rs:445`) asserts
`stderr.contains("aux tag")`, and the *predicate's own* message already begins
`"aux tag value cannot pass through the FASTQ header"`. The assertion's comment says "message must
name the tag" — it does not check that. Deviation 1's whole justification ("the message names the tag
without the predicate having to") therefore rests on an unpinned line.

The real message chain, verified with the release binary:

```
Error: test_files/ubam_lf_tagvalue.bam: BAM record 2
Caused by:
    0: aux tag 'CB'
    1: aux tag value cannot pass through the FASTQ header: "AAA\nCCC". …
```

**Change:** `stderr.contains("aux tag 'CB'")` in place of `stderr.contains("aux tag")` at
`tests/integration_ubam.rs:445`. That converts M3 from a surviving mutant to a killed one.

#### M2 — Docs: the sentence at `quickstart.md:72` inverts the reason the scope was widened

> "A read name containing whitespace is refused, **on every output format** — including FASTQ output,
> where nothing is currently lost."

For LF/CR/FF that is false, and the falsity is exactly §2.3's load-bearing fact: on FASTQ output an LF
in a QNAME emits a five-line record and desynchronises the whole file (measured pre-fix by both plan
reviewers, and re-measured by me today via the unguarded `A` arm — 13 lines for 3 records, exit 0).
§2.3 rejected an output-format gate *because* of that, so a docs sentence claiming FASTQ output loses
nothing tells the reader the refusal there is pure over-rejection.

The test comment gets this right (`// The widened scope: FASTQ output loses no tag **text** today`) —
the docs generalised it. Suggested: "…including FASTQ output, where a space or tab loses nothing today
but a newline silently breaks the record structure."

#### M3 — Docs: `test_files/README.md` documents a samtools warning that samtools does not emit

> "Verify (both report `had no targets in header` — expected for unmapped-only BAMs with no `@SQ`,
> same as `ubam_test.bam`)"

Measured with stdout and stderr redirected to **separate files** (my first probe used
`2>&1 >/dev/null`, which zsh's MULTIOS silently turns into a tee — a false negative I corrected):

| Fixture | `samtools view` exit | stderr bytes |
|---|---|---|
| `ubam_lf_qname.bam` | 0 | **0** |
| `ubam_lf_tagvalue.bam` | 0 | **0** |
| `ubam_test.bam` | 0 | **0** |
| `ubam_ws_qname.bam` | 0 | **0** |

`flagstat`, `idxstats` and `stats` are silent too. samtools 1.21 emits nothing. A reader following the
README will hunt for a warning that never appears and may conclude the fixture is malformed.

Second inaccuracy in the same block: the inline comment `# shows readB \n EVIL`. The documented
command `samtools view … | sed -n 2p | cut -f1 | od -c` yields exactly `r e a d B \n` — `EVIL` sits on
the *next* text line, which `sed -n 2p` excludes. `sed -n '2,3p'` would show both.

#### M4 — Docs: the `FS="\t"` caveat names the wrong failure mode

> "`FS="\t"` is essential — a bare `awk` splits on spaces and shifts every field."

The conclusion is right, the mechanism is not. Measured: the bare-awk output is **byte-identical to
the input**.

```
orig vs bare-awk : IDENTICAL  -> silent no-op, not a field shift
orig vs FS-awk   : name with space -> name_with_space   (one field changed, nothing shifted)
```

Why: under the default `FS`, `$1` is `name` — it *stops* at the first space by definition — so
`gsub(/ /,"_",$1)` matches nothing, `$0` is never rebuilt, and the line is reprinted verbatim. Field
shifting is unreachable through this particular pipeline. The genuine hazard is quieter and worse to
leave undescribed: the user gets an unchanged file, `samtools view -b` accepts it at exit 0, and Trim
Galore still refuses it. Suggested: "…a bare `awk` splits on spaces, so `$1` stops at the first space,
the substitution silently does nothing, and the file comes out unchanged."

(The plan §5 step 1 carries the same claim, hedged with "presumably which is why v2 concluded none
existed" — the tell that it was reasoned rather than run. Worth noting given §11's own conclusion that
every false claim in v1/v2 came from a derivation.)

The recipe **itself is correct**, verified end-to-end: `name with space` → `name_with_space`,
`CB:Z:AAACCC` intact, and the fixed BAM then round-trips through the built binary at exit 0. The
non-generalisation caveat is also correct — the LF fixture through the same pipeline dies with
`[W::sam_read1_sam] Parse error at line 4`, pipeline exit 1. (The plan said line 3; the docs wisely
quote no line number.)

#### M5 — Structure: the `reject_whitespace_in_qname` docstring restates the derivation the plan falsified

```rust
/// The set is `u8::is_ascii_whitespace`: `parse_name_and_data` splits on space/tab,
/// `fastq::read_id_prefix` on ASCII whitespace, and a FASTQ id occupies one line.
```

§2.7 established there are **five** boundary sites, not three, and that
`FastqRecord::append_to_id`'s `trim_end` uses **Unicode** `White_Space` — so the union is *not*
`is_ascii_whitespace` and VT's acceptance is a knowingly-scoped judgement, not a derivation. The
docstring presents the three-site version as if it derived the set. `qname_non_boundary_bytes_accepted`
repeats it (`// None of the three splitters named on reject_whitespace_in_qname treat these as a
boundary`).

Plan §11 "Remaining risks 3" required precisely this distinction be written down: *"an undocumented
pin invites a future reader to 'fix' it."* The pin exists and works — mutant M6 (widen to include VT)
is killed by `qname_non_boundary_bytes_accepted` — but the recorded *reason* is the falsified
derivation, so a future reader who finds `trim_end` will conclude the code is wrong rather than
deliberate, which is the exact failure mode §2.7 predicted.

**Change:** one line on the test, not the docstring (the docstring is fine as a statement of the rule
if it drops the causal "The set is X **because**" framing). E.g. `// VT is accepted deliberately: not
in is_ascii_whitespace, and rejecting it would over-reject.` Keeping it to one line satisfies
`CLAUDE.md` and makes the pin self-explaining.

---

### LOW

**L1 — Efficiency/UX: the tag-value error message is unbounded.** A 6001-byte `Z` value carrying one
newline produced **6267 bytes** on stderr. `emit_description_dropped_once` truncates at 60 chars 30
lines above in the same file (`src/bam.rs:1071-1075`), so the precedent is right there. QNAME is
bounded (BAM's `l_read_name` is a `u8`, so ≤ 254 bytes → ≤ ~1 KB escaped; measured 757 B for a
250-byte name), but tag values are not. Suggest reusing the 60-char truncation for both.

**L2 — Logic: error precedence is inverted between the two predicates, undocumented.** For QNAME the
whitespace check deliberately precedes the UTF-8 conversion — verified: `read\xffbad one` reports the
whitespace, and `read\xffbad` (no whitespace) still reports the pre-existing
`BAM record name is not valid UTF-8` unchanged. Both match §3's edge-case table. For tag values the
order is reversed, because the check sits after `from_utf8`: `CB:Z:AAA\xff\tCCC` reports
`BAM string tag is not valid UTF-8`, not the framing byte. Defensible, but the plan documents one
ordering and ships the opposite for the sibling predicate. One line in §3's table, or in the tag
docstring.

**L3 — Structure: four committed fixtures embed a session-specific absolute path, and the README
recipe does not reproduce them.** All four samtools-built fixtures carry
`CL:samtools view -b -o … /private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/d03b1224-…/scratchpad/fx/{ws,late,pair,spacetag}.sam`
in `@PG`. The README documents a *stdin* recipe, which yields 214 B against the committed 329 B — the
records are identical, only that header line differs. `samtools view --no-PG -b` is byte-identical
across three runs (162 B) and would make all six fixtures reproducible from the README, matching the LF
pair. No test reads these headers, so this is hygiene, not correctness — but the project's
reproducibility discipline is otherwise strict, and a 115-byte scratchpad path is odd in a committed
artefact.

**L4 — Docs: `test_files/README.md` intro overgeneralises.** "The offending record sits **second** in
the multi-record fixtures" — in `ubam_ws_qname_paired.bam` it is record **3** (the second *pair*),
which is what `ws_qname_refused_via_interleaved_deinterleaver` asserts. The table below is accurate;
the intro is not. "…sits after the first record/pair" fixes it.

**L5 — CHANGELOG: two sentences, and one incomplete enumeration.** `CLAUDE.md` says one sentence, no
exceptions; the entry is a bold headline plus one long compound sentence. Reported factually rather
than as a demand: it is **by far the shortest entry in `### Unreleased`** — #398 runs to four
paragraphs and #408 to one 20-line block — so tightening #415 alone changes little for a reader. Also
"a tab or newline inside a preserved `Z` tag value" omits `\r` and `\f`, which the predicate also
rejects; "whitespace other than a space" would be exact and no longer.

**L6 — Coverage: the space-in-`Z` acceptance is only tested on the FASTQ-output path.**
`clean_tag_value_with_space_still_accepted` reads `*_trimmed.fq`. I verified the more interesting
direction by hand — `--output-format ubam --preserve-tags CB` on `ubam_tagvalue_space.bam` exits 0
with `CB:Z:has a space` intact in the output BAM. Worth one extra assertion, since uBAM→uBAM is where
`parse_name_and_data` has to re-split the tail around the space.

**L7 — Optional: the A7 residue assertion pins less than it could.**
`ws_qname_refused_in_trimming_loop_ubam_out_leaves_partial` asserts `exists() && len() > 0`. Because
the Drop chain finalises unconditionally (see §5), that assertion would also pass if **zero** records
had leaked. The dangerous property A7 documents is that the residue is *indistinguishable from a
complete BAM*; asserting the record count (1) or that the BGZF EOF marker is present would pin that.
Not required — the test does catch the mutation — but it would make V3 assert the thing that makes the
residue a data-integrity problem.

---

## §5 — is the partial-output assertion flaky? No. But the mechanism is not the one you named.

You wrote: *"I believe the header write in `BamWriter::create` guarantees non-empty regardless of
whether record 1 was flushed."* **The guarantee holds; that reasoning does not.** `BamWriter::create`
writes the header into `bam::io::Writer<bgzf::Writer<BufWriter<File>>>` (`src/bam.rs:567-570`). Both
layers buffer, and one BAM header is far under either threshold, so **nothing is on disk when
`create()` returns.** Non-emptiness comes from the Drop chain: `bgzf::Writer`'s drop runs
`try_finish`. I decoded the residue to confirm:

```
total bytes: 407
ends with BGZF EOF marker: True
  block 0: offset 0,   BSIZE 379
  block 1: offset 379, BSIZE 28     <- the 28-byte BGZF EOF terminator
blocks: 2 ; samtools view -c => 1 record ; @PG chain intact (samtools -> trim_galore)
```

So the residue is a **fully finalised, valid, samtools-readable BAM** — which is A7's point, stated
more sharply than A7 states it, and the reason L7 above is worth considering.

Stability, measured:

| Scenario | Runs | Result |
|---|---|---|
| V3, `--output-format ubam`, partial present + non-empty | 25 | 25/25 present, 407–408 B (the 1-byte spread is only the output-dir path length entering `@PG CL:`) |
| `--cores 2` threaded reader, refused with the message | 40 | 40/40 |
| `--paired` de-interleaver, refused + `BAM record 3` | 40 | 40/40 |

No flakiness found. The two background-thread cases were the ones I expected to be racy and they are
not — the error is sent on the channel before the producer returns, and the consumer always drains it.

---

## §6 — mutation results

Baseline in the worktree: 419 lib / 16 clump_only_ubam / 17 integration_ubam, all green.
Run with `--no-fail-fast` (without it cargo stops at the first failing target and hides the rest).

| Mutant | Killed by |
|---|---|
| **M1** — delete `reject_whitespace_in_qname(name)?` (`:933`) | **7 tests.** `ws_qname_refused_at_sanity_check_entry_point`, `…_in_trimming_loop_ubam_out_leaves_partial`, `…_in_trimming_loop_fastq_out`, `…_via_threaded_reader`, `…_via_interleaved_deinterleaver`, `lf_in_qname_refused_end_to_end`, `clump_only_ubam_refuses_whitespace_qname`. Lib tests unaffected (by design — they test the predicate, not the wiring) |
| **M2** — delete `reject_framing_in_tag_value(s)?` (`:1098`) | `bam::tests::read_side_emitter_rejects_newline_in_string_tag`, `lf_in_tag_value_refused_end_to_end` |
| **M3** — delete `.with_context(\|\| format!("aux tag '{}'", tag_name))` (`:957`) | **NOTHING — survives all 659.** See **M1** above |
| **M4** — widen the tag predicate to reject space | `bam::tests::tag_value_space_accepted`, `clean_tag_value_with_space_still_accepted` |
| **M5** — narrow the QNAME predicate to space only | `bam::tests::qname_boundary_bytes_rejected`, `bam::tests::qname_error_escapes_invisible_bytes`, `lf_in_qname_refused_end_to_end` |
| **M6** — widen the QNAME predicate to include VT | `bam::tests::qname_non_boundary_bytes_accepted` — the §2.7 pin fires as designed |

5 of 6 killed. Every new test earns its place; the only unpinned line is the one M3 removes.

---

## §3 — call-site coverage, verified independently

**Static.** `bam_record_to_fastq` has exactly three production callers — `src/bam.rs:261`
(de-interleaver, wrapper `"{path}: BAM record {n}: {msg}"` at `:268`), `:380` (threaded single-stream,
wrapper at `:394`), `:451` (sync/Direct, wrapper at `:452`). Confirmed by grepping every
`read_record`, `BamReader::`, `open_sync_reader` and `open_threaded_reader` in `src/`, and by
confirming the only other noodles users are `src/format.rs` (which constructs `BamReader`, `:403`/
`:420`) and `src/clump_only.rs` (whose direct noodles use is inside a `#[cfg(test)]` block at
`:1531-1590`). `parse_name_and_data` has one caller, `BamWriter::write_record` (`:589`).
`append_tag_type_and_value` has one production caller (`:956`). Four entry points confirmed:
`main.rs:39` (sanity peek), `adapter.rs:115`, `adapter.rs:268`, then the loop.

**Dynamic.** All refuse, all with the whitespace message (`whitespace_msg=1` = the message appeared):

```
--hardtrim5 5                    exit=1  msg=1        --cores 4                       exit=1  msg=1
--hardtrim3 5                    exit=1  msg=1        --fastqc                        exit=1  msg=1
--clock (two BAMs)               exit=1  msg=1        --demux (SE BAM)                exit=1  msg=1
--implicon (two BAMs)            exit=1  msg=1        --paired (interleaved)          exit=1  msg=1
--clumpify --cores 2             exit=1  msg=1        --clump_only --output-fmt ubam  exit=1  msg=1
```

(My first pass showed `msg=0` for four of these; that was my flag combinations being rejected at the
CLI layer first — `--clock` needs an even file count, `--clumpify` needs `--cores ≥ 2`, `--clump_only`
+ BAM needs `--output-format ubam`. Corrected above.)

**V14 (multi-input, bad second input), the shape §7 documents but does not test:** confirmed real —
`ubam_test_trimmed.fq` 1479 B (complete), `ubam_ws_qname_late_trimmed.fq` 91 B (partial), exit 1.

**Sanity-check arm residue:** the output directory is genuinely **empty** after the refusal — not even
a trimming report. The test's `_trimmed`/`_val` glob is conservative and correct.

---

## §1 — predicate and message details

- **Reject set.** Verified with `rustc`: space/tab/LF/FF/CR → `is_ascii_whitespace() == true`;
  **VT (0x0B) → false**. Matches §2.7's table. Path form (`any(u8::is_ascii_whitespace)`) avoids
  `clippy::redundant_closure`; clippy green.
- **No false positives on non-ASCII.** Both predicates scan bytes, and UTF-8 continuation bytes are
  all ≥ 0x80, so a multi-byte character can never be misread as whitespace. `qname_non_boundary_bytes_accepted`
  covers this with `read\u{e9}one`.
- **`escape_debug` on invalid UTF-8.** `String::from_utf8_lossy` maps bad bytes to U+FFFD, which is a
  *printable* char, so `escape_debug` leaves it as a literal `�` rather than escaping it. Measured:
  `read\xFF\xFEbad\tone` prints as `read��bad\tone` — the whitespace is escaped, which is the part that
  matters. It also escapes `"` as `\"`, so an embedded quote cannot be confused with the wrapping
  quotes in the message. Sensible.
- **No temporary-lifetime hazard.** `shown` is a named binding, and I confirmed the inline form
  (`String::from_utf8_lossy(name).escape_debug()` directly in the format args) compiles and behaves
  identically. The named binding is the clearer of the two.
- **Message accuracy across paths.** The aux-tag clause is qualified twice — by `--preserve-tags` and
  by "on uBAM output" — so it reads as conditional context on the paths where no tag is lost, and the
  "a newline would break the FASTQ record structure" clause carries the FASTQ-output case. Accurate on
  all four. It cites the SAM spec as context while the predicate is deliberately *not* the spec's set
  (`@`, VT and non-ASCII are spec-forbidden but accepted) — consistent with §2.7's framing, and the
  message never claims the predicate *is* the spec.

## §2 — is rejecting `\f` in a `Z` value over-rejection? No.

Two independent reasons, both checkable: `\f` is in `split_ascii_whitespace`'s set, so
`fastq::read_id_prefix` would truncate the paired sync key at it; and `\f` (0x0C) is outside `Z`'s
`[ !-~]` grammar, so a conforming `Z` value cannot contain one — rejecting it cannot cost a legal
input. Deviation 2's correction of §2.8's own table ("it is rejected because it is a boundary byte
with no legitimate place in a tag value", not because it breaks line framing) is right, and deriving
the predicate as *QNAME-set-minus-space* rather than enumerating four bytes is the better choice —
mutants M4 and M5 confirm the two predicates cannot drift without a test going red.

## §4 — does `with_context` disturb any existing message?

No. The pre-existing Hex and Array bails are asserted only in `bam::tests` (`src/bam.rs:1438`,
`:1521`), which call `append_tag_type_and_value` **directly** and therefore never see the call-site
context. Grep found no test, fixture or doc asserting those strings through the pipeline. Full suite
green. The user-visible effect is that a Hex/Array refusal now gains an `aux tag 'XX'` layer naming
which tag caused it — an improvement, and untested in the same way M1 describes.

## §7 — fixtures

Minimal, correct, deterministic. Records verified individually with `samtools view | cat -A`:
offender placement is right in all six, and `ubam_ws_qname_paired.bam` uses flags 77/141 with the
offender as record 3, matching the de-interleaver test's assertion.

**`mk_bam_fixture.py` reproduces both LF fixtures byte-identically** with default bgzip
(htslib 1.21) — md5 match on both the compressed file and the uncompressed payload:

```
lf_qname     committed d60f1170…  regenerated d60f1170…  IDENTICAL
lf_tagvalue  committed 155a1f27…  regenerated 155a1f27…  IDENTICAL
```

So the documented recipe is exact, not merely equivalent. The generator itself is clean — the two
traps §5 step 3 flags (BGZF `BSIZE` = length − 1, and `@SQ` for genuinely aligned records) are both
avoided, and the docstring says why samtools cannot build these. The four samtools fixtures reproduce
their *records* exactly but not their bytes (**L3**).

## §9 — comment style

Nothing to flag. Every comment the diff adds is 1–3 lines plus one 8-line section banner, and the
pre-existing style in the same two files already contains 3-, 4-, 6-, 10-, 12-, 14-, 19- and 20-line
blocks — so the new comments are *shorter* than the files' norm. None reference the state being fixed:
no "previously", "no longer", "before this change", "rather than". Two phrases explain why an
assertion has a surprising shape rather than what the next line does (`// Asserted, not desired — the
cure spans all per-record bails and is tracked separately.`, `// The widened scope: FASTQ output loses
no tag text today, and is refused anyway.`); both are in tests, where that is the legitimate use of a
comment, and the second is more precise than the docs sentence it parallels (see **M2**).

---

## Recommendations, in order

| Priority | Action | Where |
|---|---|---|
| **High** | Guard the `A` (Character) arm — move the scan into the loop over the emitted tail so every arm is covered. Add a unit test + an `A`-type LF fixture | `src/bam.rs:953-957`, `:1088`; `test_files/mk_bam_fixture.py` |
| **Medium** | `stderr.contains("aux tag 'CB'")` so deviation 1's `with_context` is pinned | `tests/integration_ubam.rs:445` |
| **Medium** | Fix "where nothing is currently lost" — a newline *does* lose the whole FASTQ | `docs/…/quickstart.md:72` |
| **Medium** | Drop the `had no targets in header` claim; fix `# shows readB \n EVIL` (needs `sed -n '2,3p'`) | `test_files/README.md` |
| **Medium** | Restate the `FS` caveat as a silent no-op, not a field shift | `docs/…/quickstart.md:82` |
| **Medium** | One line recording that VT acceptance is a judgement, not a derivation from three sites | `src/bam.rs:1467-1468` |
| Low | Truncate the shown tag value (reuse the 60-char precedent at `:1071`) | `src/bam.rs:889-899` |
| Low | Note the inverted UTF-8-vs-framing precedence for tag values | plan §3 table or the tag docstring |
| Low | Rebuild the four samtools fixtures with `--no-PG` so the README recipe reproduces them byte-for-byte | `test_files/` |
| Low | "sits after the first record/pair", not "second" | `test_files/README.md` |
| Low | Trim the CHANGELOG entry to one sentence; say "whitespace other than a space" | `CHANGELOG.md` |
| Low | Assert the space-bearing `Z` value on the uBAM→uBAM path too | `tests/integration_ubam.rs` |

## Verification performed

| Check | Result |
|---|---|
| `cargo test --no-fail-fast` (full, worktree) | **659 passed, 0 failed** — independently confirms the brief |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0 — **and proven able to fail**: planted `clippy::len_zero` in the new test → exit 101 naming `tests/integration_ubam.rs:297` |
| `cargo fmt --all -- --check` | clean |
| Mutation testing, 6 mutants | 5 killed, 1 survived (**M1**) |
| Flakiness, 105 binary runs across 3 scenarios | 0 failures |
| Fixture byte-reproducibility | LF pair identical; samtools four differ in `@PG` only |
| Docs samtools/awk recipe, run verbatim | works; fixed BAM round-trips at exit 0 |
| Every specialty/parallel/paired/demux path | all inherit the refusal |
| BGZF-level decode of the partial-output residue | valid 2-block BAM, EOF marker present, 1 record |

Worktree: `$TMPDIR/codereview-B-415` (source files restored to the reviewed state after mutation).
**No file under `/Users/fkrueger/Github/TrimGalore` was modified except this report.**
