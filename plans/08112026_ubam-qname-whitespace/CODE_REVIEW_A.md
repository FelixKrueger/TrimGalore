# CODE REVIEW A — #415, whitespace in a uBAM QNAME

**Reviewer:** A (independent, fresh context, no shared state with Reviewer B)
**Target:** uncommitted working changes on `dev`, base `1fc00f6`
**Method:** isolated worktree at `$TMPDIR/codereview-A-415` (base commit + `git apply` of the
working diff + the six untracked fixtures copied in). No file in the main tree was modified.
**Date:** 2026-08-12

---

## Summary

The fix is well-targeted, the QNAME predicate is placed at the one chokepoint that every BAM
input path funnels through, and the test suite is genuinely load-bearing — I mutated the
implementation five ways and every mutation went red, including both over-rejection directions.
Verification claims in the hand-off all reproduce independently.

**One High finding blocks a clean sign-off:** the framing check guards `Z` tag values but not
`A` (`Value::Character`) tag values, which `append_tag_type_and_value` writes unvalidated one
match arm above the new call. I built a fixture and reproduced, with the patched binary, the
*identical* corruption signature the patch exists to close — **13 output lines for 3 input
records at exit 0** — and, on the uBAM-output path, Trim Galore emitting a BAM that carries a
spec-invalid `A` tag value at exit 0. The plan's own argument for folding §2.8 in ("stopping at
QNAME would ship a half-applied rule for a corruption the user cannot distinguish") applies
verbatim one arm over.

Three Medium findings are documentation accuracy: one docs sentence contradicts §2.3's
load-bearing fact, and two `test_files/README.md` claims are false as written (the samtools
recipe does not reproduce the committed fixtures; the "verify" command cannot fail).

Nothing in the Rust implementation is wrong on the paths it covers. No temporary-lifetime
hazard, no bypassed call site for QNAME, no message regression from the `with_context`
addition, and the partial-output assertion is **not** flaky.

---

## Verification I performed independently

| Check | Result |
|---|---|
| `cargo test --release`, full log captured (no `tail`) | **659 passed, 0 failed** across 15 targets — matches the hand-off exactly |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0, **0 warnings** |
| Clippy really lints test targets | Planted a `clippy::len_zero` probe in `tests/integration_ubam.rs` → **exit 101**, `tests/integration_ubam.rs:471`. Probe removed, md5 restored. The clean pass is real |
| `cargo fmt --all -- --check` | exit 0 |
| `mk_bam_fixture.py \| bgzip` reproduces both LF fixtures | **byte-identical**, md5 `d60f1170…` and `155a1f27…` |
| Docs `samtools`/`awk` remediation recipe | **works as documented** — `name with space` → `name_with_space`, `CB:Z:AAACCC` intact, and the fixed BAM round-trips through Trim Galore with the tag in the header |
| `FS="\t"` caveat | Correct in substance (see L6 for a nuance) |
| Invalid-UTF-8 **and** whitespace QNAME | Whitespace error wins, as §3's table specifies; `from_utf8_lossy` renders the bad byte as U+FFFD, no panic, no encoding error |
| Partial-BAM residue, 25 repeats | Always present, always non-empty (403–405 B). **Not flaky** — see Q5 |
| A1 (single conversion point) | Verified — 3 production callers, and the `--passthrough` carrier uses `FastqReader` directly so it never becomes a BAM record |

### Mutation testing — every mutation is caught

| # | Mutation | Tests that went red |
|---|---|---|
| M1 | delete `reject_whitespace_in_qname(name)?;` | **7**: `ws_qname_refused_at_sanity_check_entry_point`, `…_in_trimming_loop_ubam_out_leaves_partial`, `…_in_trimming_loop_fastq_out`, `…_via_threaded_reader`, `…_via_interleaved_deinterleaver`, `lf_in_qname_refused_end_to_end`, `clump_only_ubam_refuses_whitespace_qname` |
| M2 | delete `reject_framing_in_tag_value(s)?;` | **2**: `bam::tests::read_side_emitter_rejects_newline_in_string_tag`, `lf_in_tag_value_refused_end_to_end` |
| M3 | widen tag predicate to the QNAME set (reject space) | **2**: `bam::tests::tag_value_space_accepted`, `clean_tag_value_with_space_still_accepted` |
| M4 | "fix" VT — also reject `0x0B` | **1**: `bam::tests::qname_non_boundary_bytes_accepted` — the §2.7 pin works exactly as the plan claimed |
| M5 | narrow QNAME predicate to LF only | **8**: the M1 seven minus the LF case, plus `qname_boundary_bytes_rejected` and `qname_error_escapes_invisible_bytes` |

Note the division of labour is correct and deliberate: under M1 the lib unit tests stay green
(419 passed) because they exercise the predicate directly, and the integration tests catch the
missing wiring. `read_side_emitter_rejects_newline_in_string_tag` is the one unit test that
pins a *call site* rather than a predicate, which is why M2 is detectable at the lib level.

---

## Issues

### Logic

#### **HIGH — H1. `A`-type aux tag values bypass the framing check; the same corruption ships**

`src/bam.rs:1088`:

```rust
Value::Character(c) => write!(out, "A:{}", *c as char).unwrap(),
```

The new check sits in the `Value::String` arm (`:1098`). The `Value::Character` arm — the other
textual arm of the same `match`, in the same function, writing into the same FASTQ id consumed
by the same `FastqRecord::write_to` — writes its byte verbatim with no validation. `A` is
`[!-~]` per the SAM spec, so an `A` value holding `\n` or `\t` is exactly as malformed as the
`Z` case §2.8 rejects, and exactly as indistinguishable to the user.

Reachability is **identical** to the `Z` case the plan accepted as in scope: a byte-level BAM
plus `--preserve-tags`. I wrote the generator in five minutes by changing three lines of
`test_files/mk_bam_fixture.py` (`tag + b"A" + bytes([0x0A])` instead of a NUL-terminated `Z`
value). samtools 1.21 reads the result at exit 0.

Measured with the **patched** release binary (`1fc00f6 — macos/aarch64 — built
2026-08-12T10:20:31Z`), 3-record fixture, record 2 carrying `XA:A:<LF>`:

| Run | Result |
|---|---|
| `--preserve-tags XA` (FASTQ out) | **exit 0**, output id `@readB\tXA:A:` + a stray blank line → **13 lines for 3 records**. Byte-for-byte the V1b signature |
| `--preserve-tags XA --output-format ubam` | **exit 0**, and the emitted BAM's record 2 aux field is `XA` `A` `0x0A` — a spec-invalid tag value **that Trim Galore produced**. `samtools view` on our own output yields 4 lines for 3 records |
| `--preserve-tags XA` with a `<TAB>` value instead | **exit 0**, id `@readB\tXA:A:\t` — a forged empty tag field |

Control, same binary: the `Z` fixture (`CB:Z:AAA\nCCC`) is correctly **refused**, exit 1,
`aux tag 'CB'`. So this is the patch's own gap, not an unpatched build.

The uBAM-output row is worth separating out, because it weakens a claim §2.5 leans on: "Trim
Galore can never *produce* such a BAM, so the defect is strictly read-side." That is true of
read *names* (noodles validates those) and false of `A` tag values (noodles does not).

**Recommended fix** — hoist the check so it covers every arm, present and future, rather than
adding a second inline call that can drift:

```rust
fn append_tag_type_and_value(out: &mut String, value: &Value<'_>) -> Result<()> {
    use std::fmt::Write;
    let start = out.len();
    match value { /* unchanged, minus the inline reject_framing_in_tag_value call */ }
    reject_framing_in_tag_value(&out[start..])
}
```

This is the "put the invariant on the thing it constrains" shape — the constrained thing is the
appended `TYPE:VALUE` text, not one type code. Three things make it cheap:

- **All four existing tests stay green unmodified.** The two unit tests calling the predicate
  directly pass bare values; the two that assert `AAA\nCCC` appears escaped in the message still
  match, because `Z:AAA\nCCC` contains `AAA\nCCC` as a substring.
- The numeric arms get scanned too — a handful of bytes, no allocation, and they can never fail.
- The docstring then reads "the appended `TYPE:VALUE` text" rather than "a `Z` tag value".

Add one unit test on the `Character` arm (mirroring
`read_side_emitter_rejects_newline_in_string_tag`) and, if an integration case is wanted, one
more `mk_bam_fixture.py` case — but the unit test is what M2-style mutation testing needs, and
the plan's §5 step 2 reasoning for keeping things at the unit level applies here too.

FASTQ-input → uBAM-output cannot reach this: `parse_tag_value(b'A', raw)` accepts any single
byte, but a tab splits the id and a newline ends the line, so neither can arrive from FASTQ.
Read-side only, same as the `Z` case.

---

#### Verified sound — no issue

- **Q3, call-site placement.** `bam_record_to_fastq` has exactly the three production callers
  the plan names — `src/bam.rs:261` (paired interleaved de-interleaver), `:380`
  (`open_threaded_inner`), `:451` (`Direct`) — and the QNAME check is *inside* the function, so
  all three inherit it unconditionally. Every reader construction in `src/` reaches BAM records
  only via `BamReader` (`format.rs:403`/`:420` factories, used by `adapter.rs:115`/`:268`,
  `specialty.rs` ×6, `clump_only.rs:766`/`:923`/`:968-969`, `main.rs` ×7). **No QNAME bypass.**
  The `--passthrough` carrier at `main.rs:1554`/`:1593` uses `FastqReader::open{,_threaded}`
  directly, bypassing format dispatch — so a BAM carrier would be misread as FASTQ, but it also
  means the carrier never produces a BAM record. That closes the last hole in A1 rather than
  opening one.
- **Q1, the QNAME predicate.** `name.iter().any(u8::is_ascii_whitespace)` type-checks against
  `FnMut(&u8) -> bool` and is the path form clippy wants. The set is space/`\t`/`\n`/`\f`/`\r`,
  VT excluded — the documented decision. `String::from_utf8_lossy(name)` is bound to a local
  before `.escape_debug()`, so there is **no temporary-lifetime hazard**; writing it inline
  inside the `bail!` argument list is what would have been fragile. Verified end-to-end on an
  invalid-UTF-8 name (see the table above).
- **Q2, the tag predicate.** `b != b' ' && b.is_ascii_whitespace()` is right. `Z` is `[ !-~]*`,
  so space is the one byte the QNAME set over-rejects for tag values, and deriving the predicate
  by subtraction rather than enumeration is the better call — it cannot drift from the QNAME set.
  Rejecting `\f` is **not** over-rejection: `fastq::read_id_prefix`'s `split_ascii_whitespace`
  treats it as a boundary, so an `\f` in a tag value would truncate the paired sync key.
  Deviation 2 is right to correct §2.8's table on *why* `\f` is rejected.
- **Q4, the `with_context` addition.** No message regression. The Hex (`:1108`), Array (`:1117`)
  and UTF-8 (`:1097`) bails are asserted only by `src/bam.rs:1438` and `:1521`, which call
  `append_tag_type_and_value` **directly** and therefore never see the caller's context. `grep`
  over `tests/` and `docs/` finds no reference to any of those strings. The added context is a
  strict improvement — it names the tag, which the predicate deliberately does not.
- **Q5, is the partial-output assertion flaky? No.** `BamWriter::create` (`bam.rs:565-570`)
  does `File::create` then `write_header` into `bam::io::Writer`, which
  `noodles-bam-0.73.0/src/io/writer.rs:143` wraps in `bgzf::Writer::new`. `noodles-bgzf 0.35.0`
  (`src/writer.rs:235-239`) implements `Drop` as `try_finish()`, which deflates the pending
  block and emits the BGZF EOF marker; the wrapped `BufWriter<File>` then flushes on its own
  drop. So the file is deterministically non-empty even when record 1 was never logically
  flushed. 25 consecutive runs: present every time, 403–405 B. The 1–2 byte spread is not
  non-determinism — it is the `-o` path length inside the `@PG CL:` line (`part_1` vs
  `part_10`), which is fixed for a given invocation. The plan's measured "0-byte files under
  `--cores 2`" is a different arm (the parallel FASTQ writer) and is not what this test asserts.
  Residual risk is only L5 below.
- **Q6, can the new tests fail?** Yes — all five mutations caught, table above.

---

### Documentation

#### **MEDIUM — M1. The docs sentence contradicts §2.3's load-bearing fact**

`docs/src/content/docs/quickstart.md:72`:

> A read name containing whitespace is refused, **on every output format** — including FASTQ
> output, **where nothing is currently lost**.

False for `\n`/`\r`/`\f`. An LF in a read name desynchronises the entire output FASTQ (13 lines
for 3 records, exit 0) — that is precisely the fact §2.3 promotes to carry the whole scope
decision, and the reason an output-format gate was ruled out rather than merely declined. The
sentence understates the risk in the one case where the risk is total, and it is the sentence a
user reads when their working pipeline starts failing.

Suggested rewrite of that clause: "…including FASTQ output, where a space or tab loses nothing
today but a newline silently desynchronises the whole output file."

The CHANGELOG's phrasing ("also refuses FASTQ-output runs that previously succeeded") is
neutral and fine.

#### **MEDIUM — M2. The README samtools recipe does not reproduce the committed fixtures**

`test_files/README.md:73-79` presents a pipe:

```sh
{ printf '@HD\tVN:1.6\n'; printf 'name with space\t…\n'; } | samtools view -b -o test_files/ubam_ws_qname.bam -
```

Running it verbatim gives md5 `bc7b260d…`, 307 B. The committed fixture is md5 `91f0982a…`,
329 B. The difference is the `@PG CL:` line: the committed fixtures were built from a temp
`.sam` **file**, not a pipe, so their header records
`CL:samtools view -b -o test_files/ubam_ws_qname.bam …/scratchpad/fx/ws.sam` while the recipe
produces `… -`. The two byte-level fixtures, by contrast, reproduce **byte-identically** — the
inconsistency is confined to the four samtools-built ones.

Either regenerate the four through the documented pipe (so recipe and artifact agree, and L1
disappears with it), or say in the README that the recipe reproduces the *records* but not the
`@PG` line.

#### **MEDIUM — M3. The README "verify" block is wrong twice, and the command cannot fail**

`test_files/README.md:96-101`:

> Verify (both report `had no targets in header` — expected for unmapped-only BAMs with no
> `@SQ`, same as `ubam_test.bam`):
> ```sh
> samtools view test_files/ubam_lf_qname.bam | sed -n 2p | cut -f1 | od -c  # shows readB \n EVIL
> ```

1. **No such message is emitted.** With samtools 1.21, stderr is *empty* for `view`, `flagstat`
   and `stats` on all of `ubam_lf_qname.bam`, `ubam_lf_tagvalue.bam` and `ubam_test.bam`. A
   reader who runs this and sees nothing will assume the fixture is broken.
2. **The command is non-probative.** `sed -n 2p` prints one line; that line is `readB`, and the
   `\n` in the `od -c` output is *sed's own line terminator*, not the byte inside the QNAME.
   Actual output is `r e a d B \n` and stops — `EVIL` is on line 3. A BAM whose record-2 QNAME
   were plain `readB` would produce byte-identical output, so the check cannot fail. That is the
   same defect the plan calls out in v2's own validation design ("a plan that proposes a
   validation which cannot fail has not validated anything").

A command that does fail on the wrong input:

```sh
samtools view test_files/ubam_lf_qname.bam | wc -l   # 4 lines for 3 records
samtools view test_files/ubam_lf_qname.bam | sed -n '2,3p' | od -c   # readB \n E V I L \t 4 …
```

---

### Structure and style

#### **LOW — L1. Committed fixtures embed an absolute session path**

All four samtools-built fixtures carry, in their `@PG CL:` header,
`/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/d03b1224-…/scratchpad/fx/ws.sam` —
a session-scoped temp path with a UUID. Harmless to the tests, but it is committed noise a
future maintainer will puzzle over, and `ubam_tagvalue_space.bam` is the one fixture whose run
*succeeds*, so under `--output-format ubam` that path would be preserved into the output header
by the `@PG`-chain rule. Regenerating through the README's pipe (M2) fixes both at once.

#### **LOW — L2. "the three splitters" under-counts, in the exact place §2.7 warns about**

`src/bam.rs:867-868` (docstring) and `tests`' `qname_non_boundary_bytes_accepted` comment both
say the set follows from three boundary sites. §2.7 exists *because* that claim was falsified in
plan v2 — there are five, and `append_to_id`'s `trim_end` uses Unicode `White_Space`, a strict
superset. A reader who checks the docstring's derivation finds the fourth site and concludes the
comment is wrong rather than that VT-acceptance was a decision. Suggest presenting the set as
chosen: e.g. "The set is `u8::is_ascii_whitespace` — the boundary bytes of `parse_name_and_data`
and `fastq::read_id_prefix`, plus the FASTQ id's one-line framing. VT is deliberately outside
it." One line longer, and it stops inviting the correction.

#### **LOW — L3. Four three-line comments against CLAUDE.md's two-line maximum**

`ws_qname_refused_in_trimming_loop_ubam_out_leaves_partial`,
`clean_tag_value_with_space_still_accepted`, `description_notice_wording_is_fastq_input_only`,
and the `clump_only_ubam_refuses_whitespace_qname` doc comment are each three lines; the `#415`
banner in `tests/integration_ubam.rs` is seven. All state current facts — no "previously" / "no
longer" / "before this change" anywhere in `src/` or `tests/` (grepped), so the *historical
narrative* rule is clean. This is length only, and the banner's content (the three pre-scan
entry points) is genuinely non-obvious, so I would trim rather than delete. Lowest-priority item
in this review.

#### **LOW — L4. CHANGELOG entry is two sentences, the second carrying rationale**

`CHANGELOG.md:8-12`. CLAUDE.md says one sentence, no rationale, point at the docs. The second
sentence carries the read-side rationale plus the `Z`-tag change plus the docs pointer. Worth
stating plainly that it is **already by far the shortest entry in the section** — #398 above it
runs 26 lines across four paragraphs — so this is a nudge toward the standing rule, not a
regression against local practice. A one-sentence form: "**A BAM read name containing whitespace,
or a tab or newline inside a preserved `Z` tag value, is now refused on every output format**
([#415](…)) — see the uBAM section of the docs for the samtools recipe."

#### **LOW — L5. The partial-output assertion is coupled to noodles' `Drop` behaviour**

Not flaky today (Q5), but `assert!(… .len() > 0)` holds only because
`noodles-bgzf::Writer::drop` calls `try_finish()`. A future noodles bump that removed the
Drop-flush would turn this red with a message ("a mid-stream refusal leaves a partial uBAM
behind; assert the real behaviour") that misdirects to #415. And when the partial-output
follow-up lands, this assertion must be inverted. The test's own comment already says the second
part, which is the right call — consider naming the mechanism in one clause so the first is
diagnosable too.

#### **LOW — L6. The `FS="\t"` caveat is right for gawk, unreproducible on macOS awk**

`quickstart.md:82` — "a bare `awk` splits on spaces and shifts every field." Correct under gawk,
where `gsub` on a field always rebuilds `$0` with `OFS`. On the BWK awk shipped with macOS the
bare version happens to be harmless here, because `$1` is `name`, `gsub` substitutes nothing,
and `$0` is left unrebuilt. The caveat is a warning rather than a claim about output, and the
failure it warns about is real on the most widely deployed awk, so I would **keep it as is** —
noted only so nobody "corrects" it after failing to reproduce on a Mac.

---

## Recommendations, in priority order

1. **H1 — close the `A`-tag gap** before merge, via the `let start = out.len()` hoist, plus one
   unit test on the `Character` arm. All four existing tag tests stay green unmodified. This is
   the same half-applied-rule argument the plan used to fold §2.8 in, and it now has a measured
   13-lines-for-3-records reproduction plus a Trim-Galore-produced malformed BAM behind it.
2. **M1 — fix the docs sentence** so it does not deny the corruption that justifies the scope.
3. **M3 — replace the README verify block**: drop the false stderr claim, and use a command that
   can fail.
4. **M2 + L1 — regenerate the four samtools fixtures through the documented pipe**, which makes
   the recipe reproducible and drops the absolute path in one step.
5. **L2 — reword the "three splitters" derivation** as a decision, so §2.7's correction survives
   in the code and not only in the plan.
6. L3–L6 at the maintainer's discretion; L6 is a "do not change this" note.

Everything else — predicate sets, call-site placement for QNAME, error-message accuracy under
invalid UTF-8, the `with_context` addition, fixture determinism for the byte-level pair, the
remediation recipe, and the whole test design — I checked by executing it, and it holds.
