# PLAN REVIEW A — #415 uBAM QNAME whitespace (plan v2)

**Reviewer:** A (independent; no shared state with Reviewer B)
**Reviewed:** `plans/08112026_ubam-qname-whitespace/PLAN.md` v2
**Tree:** `dev` @ `1fc00f6`, read from the main tree; all code execution in a detached
worktree at `1fc00f6` (`git worktree add --detach`). Main tree untouched.
**Method:** every claim below marked **[ran]** was established by executing code, not by
reading it. Binary built once in the worktree; fixtures built with samtools 1.21 and with a
hand-written byte-level BAM writer. Predicate facts checked with `rustc` directly.

**Verdict:** the plan's diagnosis is correct and its central decision is sound. Two of its
stated facts are false and one is materially incomplete — all three are things v2 explicitly
asked a reviewer to test rather than assume, and all three change what gets implemented:

- **A7 is false** (C1). A refused run *does* leave a partial output — I have the bytes. V6 as
  written will pass for the wrong reason and pin a false invariant.
- **§5 step 3's "not constructible" is false** (C4). I built an LF-in-QNAME BAM and used it to
  demonstrate §2.3's load-bearing fact end-to-end. The stated coverage gap can be closed.
- **§2.7's splitter union is incomplete** (C2). There is a fourth boundary-consuming site, and
  it strips `\v`.

---

## 0. The four load-bearing facts, checked

| # | Plan's claim | Result |
|---|---|---|
| 1 | `noodles` refuses to write a whitespace QNAME (`InvalidInput`) | **True but understated** → C3. The real predicate is far stricter than "whitespace" |
| 2 | The two output paths differ; uBAM loses the tag, FASTQ keeps the text | **True, reproduced byte-for-byte** including the notice string. §2.1's table is accurate |
| 3 | `find(['\t', ' '])` always sees the tab first, so §2.6's unreachability holds | **True.** Verified by reading; the argument is sound. Detail below |
| 4 | A7/O2 genuinely unverified — settle by running | **Settled by running: A7 is FALSE** → C1 |

### Fact 2 — reproduced exactly [ran]

Single-record samtools fixture, `QNAME = "name with space"`, one `CB:Z:AAACCC`:

```
--output-format ubam --preserve-tags CB  → exit 0; written QNAME "name"; CB absent
   stderr: NOTE: BAM read names cannot contain whitespace, so FASTQ header text after
           the first space is not carried into uBAM output. First dropped:
           "with space<TAB>CB:Z:AAACCC"
--preserve-tags CB (FASTQ out)           → exit 0; header bytes:
   @ n a m e SP w i t h SP s p a c e TAB C B : Z : A A A C C C \n
```

§2.1's table and the §2.6 notice quotation are both correct as written. One incidental
observation the plan does not mention: BAM input produces **plain `.fq`, not `.fq.gz`** (the
input-mirrors-output gzip rule sees no `.gz` extension on `.bam`). Irrelevant to this fix, but
an integration test that globs `*_trimmed.fq.gz` on a BAM input will find nothing.

### Fact 3 — §2.6's unreachability argument holds

`str::find` with a `char` array returns the earliest match of *either* needle, so the argument
is structurally sound, not merely probable. Post-fix the id is `@{QNAME}` + `("\t" + tag)*` with
a whitespace-free QNAME, so the first `\t` (start of the first tag) necessarily precedes any
space inside a tag value. `--rename` does not disturb this: `append_to_id` (`fastq.rs:176-177`)
splices the suffix *before* the first tab. §2.6's decision to delete the possibility rather than
reword the message is well-founded. **No finding here** — the plan is right and its reasoning is
right.

---

## 1. Logic review — findings

### C1 (Critical) — A7 is false; V6 as specified is a guaranteed false pass

**Settled by running it, per the brief.** Three-record BAM, offender at record 2 (flag 20 =
unmapped + reverse, which bails inside `bam_record_to_fastq` at `bam.rs:880-882` — the same
function, same call site, same propagation path the new check will use):

| Invocation | Where it was caught (stderr) | Left on disk |
|---|---|---|
| default flags | `Auto-detecting adapter type...` | **nothing** |
| `-a AGATCGGAAGAGC` | `Scanning for poly-G content...` | **nothing** |
| `-a … --no_poly_g` | `Trimming: mid.bam` / `Output: …` | **`mid_trimmed.bam`, 329 bytes, 1 valid record (`r1`)** |
| `-a … --no_poly_g` (FASTQ out) | same | **`mid_trimmed.fq`, 88 bytes, 1 record** |
| `-a … --no_poly_g --paired` (interleaved) | producer thread, record 3 | **`il_val.bam`, 282 bytes** |

Exit was 1 in every case, and the partial BAM is *readable and structurally valid* — `samtools
view -c` reports 1. A pipeline that checks `$?` is safe; one that globs outputs, or that runs
under `set +e`, silently inherits truncated data.

**Mechanism.** `run_ubam_output_single` opens the reader at `main.rs:2147` and then calls
`BamWriter::create` at `:2148`, which does `File::create` + `write_header` immediately
(`bam.rs:565-570`). The writer therefore exists before the first record is read. What masked
this in my first two runs is that **both pre-scans read the input first** — adapter
auto-detection (skipped by `-a`/a preset) and the poly-G scan (skipped by `--poly_g`/
`--no_poly_g`) — and each reads at most `MAX_SCAN_READS = 1_000_000` (`adapter.rs:84`). So the
partial-output regime is entered whenever the offending record is past read 1 000 000, *or*
whenever both pre-scans are skipped. Neither is exotic.

**Why V6 cannot catch this.** `sanity_check_any` (`main.rs:38-47`) peek-reads record 1 through
`BamReader::open` → `next_record` → `bam_record_to_fastq` before anything else happens. With the
plan's **single-record** fixtures the new check therefore fires at sanity-check time on every
path, before any pre-scan and before any writer — so "the output directory holds only the input"
will be true, V6 goes green, and the plan records `A7: verified`. The invariant it pins is false.
Confirmed [ran]: a first-record offender leaves the output directory completely empty.

**In-repo evidence the plan did not cite.** `main.rs:402`, from #408: *"The whole run is refused
rather than decided per record, because a mid-stream refusal would leave a partial BAM behind."*
#408 hit this exact question three PRs ago and answered it by refusing up front. This plan
proposes precisely the per-record mid-stream refusal that #408 declined to permit, and inherits
the consequence #408 named.

**What to change.** Pick one, but the plan must stop asserting "nothing is written":

1. **Cheapest and honest:** restate A7 and V6 as "non-zero exit; a partial output file may
   remain when the offender is past the pre-scan window", add the multi-record fixture with the
   offender at record ≥ 2 run under `-a --no_poly_g`, and assert the *real* behaviour. Document
   the residue in §7 as a known, pre-existing property of every per-record read-side bail.
2. **Matches house style (#408):** hoist the QNAME scan to sanity-check time — but a *complete*
   pre-flight means a full pass over the input, which is a real cost the plan has not priced.
   A partial pre-flight (first N records) reintroduces exactly the false-confidence regime above.
3. **Smallest genuine fix for the residue:** unlink the output on error, or write to a temp path
   and rename on success. Broader than this issue, and it would also fix the four other
   per-record bails (aligned, reverse, secondary, supplementary) — arguably belongs in its own
   issue, in which case say so in §7 rather than claiming A7.

Whichever is chosen, note that the partial-output residue **already exists today** for the four
pre-existing per-record bails; this plan does not create it, it adds a fifth trigger. That is a
fair reason to scope the cure out — but not a reason to state the invariant backwards.

### C2 (Important) — §2.7's splitter union misses a fourth site, and that site strips `\v`

The brief asks whether the derivation is complete. It is not. There is a fourth place a read
name is consumed at a whitespace boundary:

| Missing site | Boundary set |
|---|---|
| `FastqRecord::append_to_id`, `src/fastq.rs:175` — `let id = self.id.trim_end();` | **Unicode** whitespace (`char::is_whitespace`) |

`str::trim_end` uses `char::is_whitespace`, i.e. Unicode `White_Space`, which is a **strict
superset** of `u8::is_ascii_whitespace` — it includes VT. Verified with `rustc` [ran]:

```
   VT 0x0B  u8::is_ascii_whitespace=false   char::is_whitespace=true
trim_end("a\u{b}")    = "a"
trim_end("a\u{a0}")   = "a"      (NBSP)
trim_end("a\u{2028}") = "a"      (LINE SEPARATOR)
```

Demonstrated end-to-end [ran]: hand-built fixture with `QNAME = "readtrail\v"`, run with
`--rename --clip_R1 3`, FASTQ output → **exit 0, output header `@readtrail:clip5:ACG`**. The
trailing VT is silently gone. That is precisely the harm the plan's own edge-case row invokes to
reject trailing space: *"no trimming; silently repairing a name changes read identity."*

Reachability from BAM input is real: `main.rs:383-393` refuses `--rename` + uBAM output **only
when an input is FASTQ**, so `--rename` + BAM input is admitted on both output formats.

So two sentences in the plan are false as written:
- §2.7: "their union is exactly `u8::is_ascii_whitespace`" — the union is Unicode whitespace.
- Edge-case table: "`\v` … Not in `u8::is_ascii_whitespace`, and none of the three splitters
  treat it as a boundary" — there is a fourth, and it does.

**This does not have to overturn the settled decision.** The harm needs `--rename` **and** a
clip flag **and** a *trailing* (not interior) VT — vanishingly narrow, and `\v` in a QNAME is
already pathological. Accepting VT remains defensible. What has to change is the *status* of that
acceptance: after adding the fourth row, VT's exclusion is a **judgement that a narrow residual
risk is tolerable**, not a derivation. Right now the plan presents it as derived and
drift-proof — and V5's `\v` test is explicitly framed as an anti-drift pin, so a future reader
will trust a guarantee the code does not provide. Add the row, say "VT survives every splitter
except `append_to_id`'s `trim_end`, which only strips it under `--rename` + a clip flag;
accepted", and the reasoning is sound again.

(The same row also means a trailing NBSP or U+2028 in a non-ASCII QNAME is silently trimmed
under `--rename`. Same class, same verdict — worth a clause, not a decision.)

### C3 (Important) — the write-side rule is much stricter than §2.5 says, and all three
"accepted" classes already fail hard on uBAM output

§2.5 says noodles rejects "a whitespace QNAME". The actual predicate
(`noodles-bam-0.89.0/src/record/codec/encoder/name.rs:63-67`) is:

```rust
fn is_valid(buf: &[u8]) -> bool {
    (MIN_LENGTH..=MAX_LENGTH).contains(&buf.len())   // 1..=254
        && buf != MISSING                            // != b"*"
        && buf.iter().all(|&b| b.is_ascii_graphic() && b != b'@')
}
```

`is_ascii_graphic()` is `0x21..=0x7E`. So the write side rejects **everything except printable
ASCII minus `@`** — which is all whitespace, *and* VT, *and* every non-ASCII byte, *and* `@`,
*and* names over 254 bytes. Measured on the uBAM-output path with `-a … --no_poly_g
--preserve-tags CB` [ran]:

| QNAME | FASTQ output | uBAM output |
|---|---|---|
| `read\vone` | exit 0, works | **exit 1, `invalid input parameter`**, header-only BAM left behind |
| `read@one` | exit 0, works | **exit 1**, same |
| `readéone` (valid UTF-8) | exit 0, works | **exit 1**, same |

Consequences for the plan:

1. The edge-case table's justification for all three — "**accepted** — spec-invalid but harmless
   here, and rejecting it would break files that work today" — is **false on the uBAM-output
   path**. Those files do not work there today; they abort with `invalid input parameter`, an
   error that never mentions the read name, after leaving a header-only BAM. The three rows
   should say "accepted by this check; uBAM output still fails at the writer".
2. **V5 must not become an integration test.** As unit tests over the predicate it is correct. If
   an implementer reads "accepted classes stay accepted" and writes a `\v`/`@`/non-ASCII
   uBAM→uBAM round-trip case, it will fail, and the natural-looking fix is to widen the reject
   set — undoing a maintainer decision by accident. Say "unit tests only" explicitly.
3. §2.5's conclusion ("the defect is strictly read-side") still holds, and the sentence "neither
   fixture can be built with our own code" is right — just for a broader reason than whitespace.
4. Adjacent cheap win, if wanted: the same call site could turn `invalid input parameter` into a
   real message for these classes. Out of scope is a fine answer; claiming they are harmless is
   not.

### C4 (Important) — an LF-in-QNAME fixture *is* constructible; the stated coverage gap can be closed

§5 step 3: *"A newline-in-QNAME fixture is not constructible, and that is a stated gap."* What is
not constructible is a newline in QNAME **via SAM text**, which is what the plan tested. The BAM
container has no such limitation, and a byte-level writer is short. I built one [ran] — ~20 lines
of Python emitting the uncompressed BAM stream, piped through `bgzip`:

```
python3 mkbam.py 'read\none' | bgzip > lf_qname.bam     # 121 bytes
samtools view lf_qname.bam | od -c
  → r e a d \n o n e \t 4 \t * \t 0 \t 0 …               # LF intact in QNAME
```

Record layout used: `block_size i32`, `refID -1`, `pos -1`, `l_read_name`, `mapq 0`, `bin 4680`,
`n_cigar_op 0`, `flag 4`, `l_seq`, `next_refID -1`, `next_pos -1`, `tlen 0`, NUL-terminated name,
4-bit-packed seq, raw-phred qual, aux. Header `BAM\1` + `l_text` + `@HD VN:1.6\n` + `n_ref 0`.
samtools reads it without complaint (`quickcheck` only notes the absent `@SQ`, expected for an
unmapped-only BAM).

**With that fixture, §2.3's load-bearing fact is now demonstrated rather than argued** [ran].
Three records — `clean_one`, `bad\nread`, `clean_two`, each with `CB:Z:AAACCC` — run with
`--preserve-tags CB`, FASTQ output, at `1fc00f6`:

```
exit=0
wc -l < lf3_trimmed.fq  →  13        # three records must be 12 lines
@ b a d \n r e a d \t C B : Z : A A A C C C \n   ← one record occupying five lines
```

Exit 0, silent, and every record after it is desynchronised for any 4-line-block reader. This is
the single fact the entire scope decision (§2.3, O1) rests on, and v2 supports it only by reading
`write_to` plus a proposed unit test — the plan says so honestly, which is why it is worth
closing.

**Recommendation:** commit a third fixture (`test_files/ubam_lf_qname.bam`, and ideally the
3-record `ubam_lf_qname_multi.bam` above, which doubles as C1's mid-stream fixture), record the
byte-level recipe in `test_files/README.md` next to the samtools recipes, and add a V row
asserting the LF case is refused end-to-end on **FASTQ output**. That converts the plan's most
important claim from inference to regression coverage, and it removes the one gap §5 step 3 asks
the reader to accept. Reframe the paragraph: SAM text is line-delimited, BAM is not.

### C5 (Important) — no validation on `Z` tag values, which land in the same id

`append_tag_type_and_value` (`bam.rs:1058-1061`) writes `Z:{s}` verbatim into the id with no
whitespace check. So under `--preserve-tags`:

- an **LF** inside a preserved `Z` value corrupts FASTQ framing exactly as C4's QNAME case does —
  same function, same id, same `write_to`;
- a **TAB** inside one splits the tag tail into a bogus field, which on uBAM output surfaces as
  `tag name must be 2 chars` (`bam.rs:748`) rather than anything about the input.

Both are spec-invalid (`Z` is `[ !-~]*`) and neither is reachable from SAM text — but that is the
same reasoning that failed for QNAME in C4, and the byte-level writer reaches both. No proposed V
row covers either, and the plan's §2.3 argument ("an LF desynchronises the whole output file")
applies verbatim.

Not a reason to widen the fix — a space in a `Z` value is legitimate and common, so the check
cannot simply run over the assembled id without care. But the plan should either extend the scan
to the appended tag text (still one pass, O(name + tags), no allocation) or name this in §7 as a
deliberate out-of-scope sibling with the reason. Silence is what turns it into the next issue.

### Minor logic notes

- **Check placement verified.** The empty-name check is at `bam.rs:894-896`, `from_utf8` at
  `:897`, `format!("@{}")` at `:898`. Inserting after `:896` is correct, and the edge-case row
  "whitespace + invalid UTF-8 reports the whitespace" follows. `rec.name()` yields `&BStr`, which
  derefs to `[u8]`, so `reject_whitespace_in_qname(name)` coerces without ceremony.
- **Caller coverage verified.** Three `bam_record_to_fastq` callers (`:261`, `:380`, `:451`);
  every BAM reader construction in `src/` goes through `BamReader::open*`
  (`format.rs:403`/`:420`, `main.rs:39`/`:1848`/`:1871`/`:2439`, `clump_only.rs:923`). Specialty
  modes use `format::open_sync_reader` (`specialty.rs:53, 96, 151, 203, 268-269, 381-382`);
  `--clump_only` uses it at `:766`/`:968` plus the interleaved reader at `:923`. The
  `FastqReader::open` calls at `clump_only.rs:290/430/1132` and `specialty.rs:551` are
  FASTQ-only/test paths. **§3 item 4 and A1 are correct.**
- **Add `sanity_check_any` to §2.4.** `main.rs:38-47` peek-reads record 1 through the `:451`
  reader, so it is the earliest place the new check can fire. It is not a fourth caller of
  `bam_record_to_fastq`, but it *is* a fourth entry point, and it is exactly what separates C1's
  two regimes. One sentence in §2.4 saves the next reader the confusion I had.
- **Producer-thread propagation is fine** [ran]. Interleaved fixture with a bad record at index
  3: exit 1 within a second, `Error: il.bam: BAM record 3: …`, no hang, no deadlock on the
  bounded channel. §5 step 4's note about asserting on process exit/stderr is right.
- **Error context comes free.** All three callers wrap with `"{}: BAM record {}"`, so the user
  sees `Error: file.bam: BAM record 12345` / `Caused by: BAM read name contains whitespace: …`
  [ran]. The record index is the practical remediation handle on a large file and is worth
  mentioning in §5 step 1 — it partly answers the "no remediation command" gap.

---

## 2. The §5 step 1 error message, clause by clause

The brief asks whether every clause is true on every path the check fires on. Close, but two
clauses are still over-scoped:

| Clause | True everywhere? |
|---|---|
| "BAM read name contains whitespace: `…`" | Yes |
| "A FASTQ header ends the read name at the first whitespace" | True of *our reader*; a user running FASTQ output can look at the file and see the whole name present (§2.1 row 2), so it reads as contradicted by their own output. Suggest "Trim Galore reads a FASTQ header as ending at the first whitespace" |
| "aux tags following it are lost on uBAM output" | Only with `--preserve-tags`. Without it there are no tags, so the clause describes a loss that is not occurring. Suggest "any `--preserve-tags` aux tags following it are lost on uBAM output" |
| "a newline would break the FASTQ record structure" | Yes, correctly hedged — and now demonstrable (C4) |
| "The SAM specification forbids whitespace in QNAME" | Yes |
| "Rename the reads before trimming." | Yes; the "no remediation command" decision is right, and I could not find a clean samtools one-liner either |

The v1 → v2 rewrite fixed the real defect (the unconditional tag-loss assertion). What remains is
two scope qualifiers, both one word.

---

## 3. Assumptions

| | Verdict |
|---|---|
| A1 | **Verified** [ran/grep]. Add `sanity_check_any` as the earliest entry point |
| A2 | **Falsified as a derivation** (C2). The union is Unicode whitespace; `append_to_id`'s `trim_end` is the missing site. Accepting VT stays defensible as a judgement |
| A3 | Not relitigated. C4 strengthens it: the fact it rests on is now demonstrable end-to-end |
| A3b | **Supported by the evidence I could find** — see below |
| A4 | **Verified** [ran]. samtools 1.21 builds and reads back the space-in-QNAME fixture |
| A5 | **Verified but understated** (C3) |
| A6 | **Verified** by reading `find`'s semantics; the argument is structurally sound, not probabilistic |
| A7 | **FALSE** (C1) |

**A3b, with evidence.** The brief asks whether anything can be found either way on
header-in-QNAME uBAMs being rare. What I could test [ran]: `samtools import` on
`@A00690:224:H7GC5DRX2:1:2101:1081:1016 1:N:0:CGATCG` yields `QNAME =
A00690:224:H7GC5DRX2:1:2101:1081:1016` — the description is dropped, not folded in. Same with
`-T '*'`. `samtools fastq` round-trips it back identically. So **the dominant converter cannot
produce the shape A3b calls rare**, which supports the judgement. I could not test Picard
`FastqToSam` or bedtools locally, so that is one converter checked, not the field. A3b remains a
judgement — but it is now a judgement with one supporting data point rather than none, and that
is worth a line in the assumption row.

---

## 4. Efficiency

Nothing to fault. `name.iter().any(u8::is_ascii_whitespace)` is one pass over ≤ 254 bytes with no
allocation on the happy path, against a per-record BAM decode that already touches the sequence
and quality arrays. Immeasurable.

Two small confirmations: the path form does type-check (`u8::is_ascii_whitespace` is
`fn(&u8) -> bool`, which is exactly what `Iterator::any` wants for `Iterator<Item = &u8>`), and
§6's note about `clippy::redundant_closure` under `-D warnings` is the right reason to prefer it.

If C5 is taken, scanning the assembled id after tag appending is still one pass over
name + tags — but it cannot reuse the same predicate, since a space in a `Z` value must stay
legal. Two predicates, or a scan restricted to the reject set minus space, would be needed. That
cost is a legitimate reason to defer C5; it should be stated rather than left unexamined.

---

## 5. Alternatives

- **The class split (§2.3's named fallback)** stays the right fallback and is priced correctly.
  C1 does not change that; C4 strengthens the case against it, since the LF harm is now
  demonstrated rather than argued.
- **Sanity-check-time refusal** (#408's house style) is the only option that makes A7 true as
  written, and it is the one the plan does not consider. It costs a full extra pass over the
  input — real, and not obviously worth it. Worth one line in §10 as considered-and-declined, so
  the next reader does not have to rediscover why the check is per-record.
- **Unlink-on-error / write-to-temp-then-rename** is the smallest change that removes the residue
  in C1 and would also cover the four existing per-record bails. Probably its own issue; the plan
  should say so rather than assert A7.
- **Rewriting the splitter to prefer the tab** — §11 already rejects this for the right reason
  (the FASTQ-input semantics are correct and would break). Agreed, no note needed.

---

## 6. Validation sufficiency

The 10+ V rows cover the predicate thoroughly. The gaps are all in the *call-site* dimension:

| Gap | Consequence |
|---|---|
| **Every fixture is single-record**, so `sanity_check_any` catches record 1 on every path | V6 passes for the wrong reason and pins a false invariant (C1). No row exercises the trim-loop regime at all |
| **No row varies the pre-scan flags** | The three regimes in C1's table are behaviourally different and only one is tested. `-a … --no_poly_g` is the arm that reaches the writer |
| **No end-to-end LF case** | The fact carrying §2.3 is tested only in a unit test — avoidable, per C4 |
| **No `Z`-tag-value case** | C5: the identical framing corruption via `--preserve-tags`, uncovered |
| **V5's "accepted" is ambiguous** | Read as an integration assertion it is false on uBAM output (C3) |
| **No long-QNAME case** | noodles rejects > 254 bytes on write with the same opaque `invalid input parameter`. Pre-existing, adjacent, not this plan's job — but it is the same message and the same call site |

The highest-risk failure mode not covered by any V row: **a whitespace QNAME beyond the pre-scan
window (or under `-a --no_poly_g`) on uBAM output, which fails with a partial BAM on disk.** Every
proposed row would go green with that behaviour in place.

---

## 7. Action items

### Critical

1. **Correct A7 and V6.** A refused run *can* leave a partial output — 329-byte BAM with one
   valid record, 88-byte FASTQ, 282-byte paired BAM, all measured. Add a multi-record fixture
   with the offender at record ≥ 2, run it under `-a … --no_poly_g`, and assert the real
   behaviour. Then pick one of C1's three resolutions and record it. Cite `main.rs:402` — #408
   already answered this question for its own case.
2. **Add the fourth splitter row to §2.7** (`fastq.rs:175`, `char::is_whitespace`) and demote
   VT's acceptance from derived to judged, keeping the decision. Also fix the two sentences that
   assert the union is exactly `u8::is_ascii_whitespace`.

### Important

3. **Rewrite the three "accepted" edge-case rows** (`\v`, non-ASCII, `@`): accepted by this
   check, still fatal at the noodles writer on uBAM output, with `invalid input parameter` and a
   header-only BAM. Quote noodles' real predicate in §2.5. State that V5 is unit-test-only.
4. **Replace §5 step 3's "not constructible" with the byte-level recipe**, commit
   `test_files/ubam_lf_qname*.bam`, and add an end-to-end V row for the LF case on FASTQ output —
   §2.3's load-bearing fact deserves regression coverage, and it is now cheap.
5. **Resolve C5**: either extend the check to the appended tag text (noting a space in a `Z`
   value must remain legal, so it is not the same predicate) or name it in §7 as a deliberate
   out-of-scope sibling with the reason.
6. **Two words in the error message**: "any `--preserve-tags` aux tags", and attribute the
   first-whitespace rule to Trim Galore's reader rather than to FASTQ headers in general.

### Optional

7. Name `sanity_check_any` (`main.rs:38-47`) in §2.4 as the earliest entry point.
8. Mention in §5 step 1 that the callers supply `file: BAM record N`, which is the practical
   remediation handle in place of the absent samtools one-liner.
9. Add A3b's one supporting data point: `samtools import` splits at the first whitespace and
   cannot produce a header-in-QNAME uBAM (verified; Picard/bedtools unchecked).
10. Record sanity-check-time refusal in §10 as considered-and-declined, with the full-extra-pass
    cost as the reason.
11. Note that BAM input yields plain `.fq`, not `.fq.gz`, so FASTQ-output integration assertions
    must not glob `*.fq.gz`.

---

## Appendix — reproduction

Worktree: `git worktree add --detach "$TMPDIR/review-A-wt-415" 1fc00f6`; `cargo build --release`.
Scratch (fixtures, stderr captures, `mkbam.py`, `mkbam3.py`, the `rustc` predicate probe) under
`<worktree>/revA_scratch/`. samtools 1.21. Main tree was never built in and never modified; the
worktree is removable with `git worktree remove --force`.
