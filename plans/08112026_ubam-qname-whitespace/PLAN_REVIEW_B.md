# PLAN_REVIEW_B — #415, plan v2 (reject whitespace in a uBAM QNAME)

**Reviewer:** B (independent, fresh context; no shared state with Reviewer A)
**Reviewed:** `plans/08112026_ubam-qname-whitespace/PLAN.md` v2, against `dev` @ `1fc00f6`
**Method:** read the cited code at `1fc00f6`, then **ran** the probes below. Everything marked
*(ran it)* was executed; everything marked *(read it)* is code-reading only. Nothing in
`/Users/fkrueger/Github/TrimGalore` was modified — tracked tree verified clean before and after
(`git status --porcelain` shows only the 12 pre-existing untracked entries). No `cargo build` or
`cargo test` in the main tree; behavioural probes used the already-built
`target/release/trim_galore` (stamped `3008b37`), which is behaviourally current for these paths:
`git diff 3008b37..HEAD -- src/` touches only `src/io.rs`, adding an env-gated test hook, and does
not touch `bam.rs`, `main.rs`, `adapter.rs` or `fastq.rs` *(ran it)*. Scratch lives in
`$TMPDIR/revB`.

**Verdict:** the plan's diagnosis is right, the chosen fix is right, and the two facts the
maintainer flagged as load-bearing both hold. But **A7 is false, not merely unverified**, and the
validation section cannot discover that, because every proposed integration case is intercepted
three entry points upstream of the code it means to test. Two Critical items, five Important.

---

## 0. What I confirmed of the plan's own claims

| Plan claim | Status |
|---|---|
| `bam_record_to_fastq` at `bam.rs:871`; empty-name check `:894-896`; UTF-8 at `:897`; id built `:898` | correct *(read it)* |
| `parse_name_and_data` `:707`, splits on `find(['\t',' '])` at `:709`; exactly one production caller, `BamWriter::write_record` `:589` | correct *(read it)* |
| `read_id_prefix` at `fastq.rs:208`, `split_ascii_whitespace` at `:210` | correct *(read it)* |
| `write_to` `fastq.rs:78-96` writes the id verbatim then `\n` | correct *(read it)* |
| §2.7's predicate table (VT excluded from `is_ascii_whitespace`; `split_ascii_whitespace` does not split on VT) | correct *(ran it,* `rustc` *)* |
| §2.1 row 1 — uBAM out: exit 0, QNAME `name`, `CB` gone, notice echoes `"with space\tCB:Z:AAACCC"` | reproduced exactly *(ran it)* |
| §2.1 row 2 — FASTQ out: `@name with space\tCB:Z:AAACCC` written verbatim | reproduced exactly *(ran it)* |
| A4 — samtools 1.21 builds the space fixture and reads it back intact | confirmed *(ran it)* |
| A5 — noodles refuses to write a whitespace QNAME | confirmed *(ran it —* `Error: failed to write BAM record / Caused by: invalid input parameter` *)* |
| §2.6 — after the fix, `find(['\t',' '])` sees the tab first, so the notice is FASTQ-input-only | holds *(read it)*. The id is `@{QNAME}` + `\t`-prefixed tags; with QNAME whitespace-free the first hit must be the first tab. `append_to_id` splices `:clip5:SEQ` (no whitespace) before that tab, and `--demux`/specialty modes never reach `parse_name_and_data` (rejected for uBAM output). A6 is sound |
| §2.4 — three production callers of `bam_record_to_fastq`, one check covers all BAM input | correct as a statement about call sites; **the labels are swapped and the entry points are four** — see C2/C3 |
| §4's signature is implementable as written | yes *(read it)*. `rec.name()` → `Option<&BStr>` (noodles-bam 0.73.0 per `Cargo.lock`, `record.rs:140`) and `BStr: Deref<Target=[u8]>` (bstr 1.12.3, `impls.rs:551`), so `&BStr` coerces at the call; no `.as_ref()` needed |
| Coverage of `--clump_only` and the specialty modes | genuine *(read it)*. `clump_only.rs:766` and `:923`, `specialty.rs:53/96/151/203/268/381` all enter through `format::open_sync_reader` or `BamReader::open_paired_interleaved_with_tags`, so all BAM input funnels through the three call sites |

---

## 1. Critical

### C1 — A7 is false. A refused run leaves a partial output file, and the partial uBAM is a valid, samtools-readable BAM. *(ran it)*

O2 asked whether the writer opens before the first record is read, and asked for it to be settled by
running rather than reasoning. I ran it, using the existing error at the **same call site**
(`bam.rs:874`'s aligned-BAM bail inside `bam_record_to_fastq`) so the propagation path is identical
to the proposed whitespace refusal. Fixture: record 1 clean unmapped, record 2 genuinely aligned
(needs an `@SQ` line — htslib silently forces `FUNMAP` when `RNAME` is `*`, which will bite anyone
building this fixture). Pre-passes suppressed with `-a AGATCGGAAGAGC --no_poly_g` (see C2 for why
that is necessary).

| Run | Exit | Left behind |
|---|---|---|
| SE, uBAM out | 1 | `alnlate_trimmed.bam`, 319 B — **valid BGZF, EOF block present, samtools reads it without complaint**, contains record 1 plus our `@PG` line |
| SE, FASTQ out | 1 | `alnlate_trimmed.fq`, 1 of 2 records |
| SE `--cores 2`, FASTQ out | 1 | 0-byte `alnlate_trimmed.fq` |
| `--paired` interleaved, FASTQ out | 1 | two 0-byte `ilv_val_{1,2}.fq` |
| `--paired` interleaved, uBAM out | 1 | `ilv_val.bam`, 293 B |

So the writer is opened before the read loop on every path, and O2's second branch — "the plan needs
a cleanup step or the refusal needs to move earlier" — is the live one. The partial BAM is the worst
shape: it is structurally complete (BGZF EOF terminator written on drop), so nothing downstream flags
it. The only signal is that the trimming reports are absent, which is easy to miss in a batch.

Note that "move the refusal earlier" does **not** fix it: a sanity-check-time peek sees only record 1
(`main.rs:38-46`). The durable fix is to remove the output paths on failure, or write to a temp name
and rename on success. This repo already carries a scar from partial/destructive output (#409), so I
would price the cleanup into this plan rather than leave it as an open non-blocking question.

Also worth knowing, from the same probes: an **LF** QNAME on the uBAM-output path already fails today
with `Error: failed to write BAM record / Caused by: invalid input parameter` and a 275-byte partial
BAM. The new check therefore also converts an opaque write-time failure into a named read-side one —
a benefit the plan does not claim.

### C2 — Every proposed integration case is intercepted by `sanity_check_any`, so V6 passes for the wrong reason and V7 is a false pass. *(ran it)*

There are **four** entry points into the three call sites, not one, and three of them run before the
trimming loop:

1. `main.rs:243` — `sanity_check_any(&cli.input[0])`, before all dispatch. For uBAM this is
   `BamReader::open` + `next_record()` (`main.rs:38-46`) → record 1 through the `:451` Direct caller.
2. `adapter.rs:115` — adapter auto-detection, up to 1M records, via `open_sync_reader` → `:451`.
3. `adapter.rs:268` — `detect_poly_g`, up to 1M records, when an explicit adapter was given and
   `--no_poly_g` was not → `:451`.
4. Only then the trimming loop (`:451` sync, `:380` threaded, `:261` interleaved producer).

Observed *(ran it)*:

- Bad **record 1**: dies before `Auto-detecting adapter type...`, output dir **empty**.
- Bad record 2 with `-a` but no `--no_poly_g`: dies during `Scanning for poly-G content...`, output
  dir **empty**.
- Bad record 2 with `-a` **and** `--no_poly_g`: reaches the loop — and that is the only configuration
  that produces C1's partial file.

Consequences for §9, all of them concrete:

- **V6 is a false pass.** With the single-record `ubam_ws_qname.bam` fixture the refusal happens at
  the sanity check, before any writer exists, so "nothing written" holds — and the plan would record
  A7 as confirmed while the hazard in C1 is real and untested. This is the "prove the check can fail
  before trusting that it passed" trap.
- **V7 is a false pass.** I tested the mechanism directly by putting the aligned record **first** in
  an interleaved BAM and running `--paired`: the error came back as `ilvfirst.bam: BAM record 1` +
  `Caused by:` — the sanity check / Direct shape (`bam.rs:451`) — not the `:261` producer's one-line
  wrapper `{path}: BAM record {n}: {msg}` (`:268`). The producer thread never starts, so the
  background-thread caveat the plan carefully adds to that test is moot and the path stays uncovered.
- **`:380` is covered by nothing.** Reaching it needs `--cores ≥2` *and* the offending record past
  both pre-passes.

Fix, and it is cheap: give both fixtures a **clean first record/pair** with the whitespace on a later
one, and run the cases with `-a` (plus `-a2` for paired) and `--no_poly_g`. Add a `--cores 2` case
for `:380`. Keep one record-1 fixture as the sanity-check arm — that is a real path, just not the
only one — and say in the test name which entry point each case pins. On a real >1M-read file the
pre-passes stop being a shield anyway, so the `-a`/`--no_poly_g` cases are the faithful ones.

### C3 — §2.4's labels for `:261` and `:380` are swapped. *(read it, confirmed by error shape when ran)*

`:261` sits in `open_paired_interleaved` (the `r1_queue`/`r2_queue`/`MAX_SLACK` de-interleaver, error
wrapper at `:268`); `:380` sits in `open_threaded_inner` (single stream, wrapper at `:394`). §2.4
calls `:261` the "threaded reader" and `:380` the "interleaved producer", and §5 step 4 then asks for
the paired test to exercise "the `:380` producer path". On its own this is a typo; combined with C2
it is what makes the paired test's design wrong, so fix it before the test is written.

---

## 2. Important

### I1 — "A newline-in-QNAME fixture is not constructible" is false. I built one. *(ran it)*

§5 step 3 accepts a stated coverage gap on the grounds that samtools cannot express an LF in SAM text
and noodles refuses the write. Both are true, and both are avoidable: BAM's QNAME is a
NUL-terminated byte field, so a ~25-line Python **stdlib** generator writes the record directly and
wraps it in a BGZF block (gzip with the `BC` extra field; `BSIZE` = block length − 1 — off-by-one
here is what produces htslib's `CRC32 checksum mismatch`). samtools 1.21 reads the result back at
exit 0 and prints the raw LF bytes.

With a three-record fixture (`readA`, `readB\nEVIL`, `readC`) on today's binary, FASTQ output is
**13 lines for 3 records** at exit 0 — record 2 emits a five-line block and every 4-line-block
consumer desynchronises from there on. So §2.2's symptom 2, the fact promoted at review to carry the
whole scope decision in §2.3, is verifiable end-to-end both before the fix (corruption) and after
(refusal). The plan's honest caveat — "§2.3's load-bearing fact is verified by reading `write_to`
plus a unit test, not by an end-to-end run" — can simply be removed.

Cost: the fixture is only reproducible from a committed script. That fits the existing convention —
`test_files/README.md:37-50` already carries samtools regen recipes — so add the generator there or
alongside it.

### I2 — §2.7's splitter union is incomplete: `append_to_id` is a fourth site, and its boundary set is a strict superset. *(ran it,* `rustc` *; read it for reachability)*

`FastqRecord::append_to_id` starts with `let id = self.id.trim_end();` (`src/fastq.rs:175`).
`str::trim_end` uses `char::is_whitespace` (Unicode White_Space), which strips **VT (U+000B), NBSP
(U+00A0) and U+3000** — I verified all three. Reachability from BAM input is real: `trimmer.rs:261-274`
calls it when `config.rename` is set and a `--clip_*` actually clipped, and #408 leaves
`--rename` admitted for BAM-in → FASTQ-out. So a QNAME with a **trailing** `\v` — accepted by
design — silently loses that byte from the output read name under `--rename`, which is precisely the
"silently repairing a name changes read identity" harm the edge-case table refuses elsewhere.

This does not overturn the settled predicate: rejecting VT read-side would over-reject on the same
grounds the plan gives, and trailing-VT-plus-`--rename` is vanishingly rare. But A2's claim that the
union is **exactly** `u8::is_ascii_whitespace` is not true as written, and V5 pins VT-acceptance as
an *anti-drift* test — so the next person to find `trim_end` will read the plan as wrong rather than
as knowingly scoped. Add the row with the "trailing only, `--rename` only, accepted collateral" note.
(`demux.rs:223/240`'s `record.id.replace(' ', "_")` is a fifth site, but space-only, so it does not
extend the union — worth a mention because it means demux currently *renames* such reads silently.)

### I3 — "`@` and non-ASCII UTF-8 … both harmless" is false on the uBAM-output path. *(ran it)*

All three classes §2.7 accepts as harmless fail today on uBAM output, exit 1, with
`Error: failed to write BAM record / Caused by: invalid input parameter` and a partial
`*_trimmed.bam` left behind; the clean control passes:

| QNAME | uBAM out | FASTQ out |
|---|---|---|
| `name\x0bvt` (VT) | exit 1, partial BAM | exit 0 |
| `na@me` | exit 1, partial BAM | exit 0 |
| `namée` (non-ASCII) | exit 1, partial BAM | exit 0 |
| `plainname` | exit 0 | exit 0 |

noodles enforces the SAM QNAME charset `[!-?A-~]{1,254}`, which excludes `@` (0x40), every byte
above `0x7E`, VT, and all whitespace. So the writer is a fourth validity constraint, and it is
neither a subset nor a superset of the reject set. Two consequences, neither a scope change:

1. The edge-case table should read "accepted read-side; on uBAM output noodles rejects the write
   today with an opaque message — pre-existing, out of scope" instead of "harmless". As written, an
   implementer could reasonably add an integration test asserting a `@`- or non-ASCII-named record
   round-trips to uBAM output, and be confused when it fails.
2. It mildly **strengthens** A3's cost argument: since VT/`@`/non-ASCII already fail on the uBAM
   path, the only class where this plan causes a genuinely new failure is space/tab/CR/FF on the
   FASTQ-output path — which is exactly the class §2.3 prices and accepts.

A follow-up issue for the opaque noodles message (no read name, partial output) seems worth filing;
it is the same defect shape as #415 with a different trigger.

### I4 — no CHANGELOG entry in §5 or §7. *(ran it)*

Every one of the last six commits on `dev` touched `CHANGELOG.md`, including #408's refusal
(`38deb1e`) and #398 (`11dfc38`), whose `### Unreleased` entry carries an explicit "**This refuses
some invocations that previously ran**" paragraph. This change is exactly that shape — the plan's own
header calls it "a currently-succeeding input class starts failing" — yet §5's six deliverables and
§7's file list ("`src/bam.rs` only, plus tests, one committed fixture, one docs sentence") both omit
it.

### I5 — the "no clean samtools one-liner" claim is falsifiable for the reported case. Shown running, as §5 step 1 asks. *(ran it)*

```bash
samtools view -h in.bam \
  | awk 'BEGIN{FS=OFS="\t"} /^@/{print; next} {gsub(/ /,"_",$1); print}' \
  | samtools view -b -o fixed.bam -
```

On the space fixture this produced `name_with_space` with `CB:Z:AAACCC` intact, and `fixed.bam` then
ran through `trim_galore --output-format ubam --preserve-tags CB` with the tag preserved. The whole
trick is `FS="\t"`, so `$1` is the entire QNAME; a bare `awk` splits on spaces and shifts every
field, which is presumably the trap behind the plan's claim.

It does not generalise, and I checked rather than assumed: I built tab-QNAME and LF-QNAME BAMs, and
`samtools view -h | samtools view -b` fails on both with `Parse error at line 3` — SAM text cannot
represent either field. So the accurate statement is "no one-liner covers tab or newline names",
which leaves the reported shape (a space) with a working recipe. Since the error text is what a user
acts on, I would offer it there with the caveat, or in the docs sentence from step 6.

---

## 3. Optional

- **O-a. Two severity facts worth folding into §2.1/§2.2** *(ran it)*. (i) `--clump_only
  --output-format ubam` on the space fixture also drops the tag silently at exit 0 — and
  `--clump_only` advertises byte-identical preservation, which is the exact ground on which
  `cli.rs:853-855` refuses `--rename`. A documented-lossless mode silently violating its contract is
  a stronger argument than anything currently in §2.2, and it upgrades §5 step 5 from a formality to
  a real guard. (ii) A multi-input SE run whose *second* input holds the bad record leaves input 1's
  outputs complete and input 2's partial — the sanity check only covers `cli.input[0]`. Worth one V
  row alongside C1's cleanup.
- **O-b. §5 step 6's "uBAM input section" does not exist by that name.** The input-focused section is
  `docs/src/content/docs/quickstart.md:55` "## Unaligned BAM (uBAM)";
  `docs/src/content/docs/guide/outputs.md:54` is output-only. Name the file.
- **O-c. §5 step 3's "the other five `test_files/ubam_*.bam` fixtures"** is five only if you count
  `ubam_out_{se,pe}_REFERENCE.bam`, which are TrimGalore outputs rather than samtools-built. The
  samtools-built inputs are three: `ubam_test.bam`, `ubam_test_with_tags.bam`,
  `ubam_paired_test.bam`.
- **O-d. V6's wording doesn't match the harness.** `tests/integration_ubam.rs::fresh_tmpdir` hands
  each test an empty directory and the fixture lives in `test_files/`, so "the output directory holds
  only the input" is not expressible. Assert that no `*_trimmed.*` / `*_val*` exists.
- **O-e. Escaping.** §5 step 1 asks for `\t` → `\\t` and `\n` → `\\n`. Add `\r` and `\f`, or the
  message for a CR-only name prints a name that looks unremarkable and overwrites itself in a
  terminal. `escape_debug` on the lossy string gives all four for free.

---

## 4. Efficiency

Nothing to add: `name.iter().any(u8::is_ascii_whitespace)` over ≤254 bytes against a per-record BGZF
decode is free, the path form avoids `redundant_closure`, and the escaping is error-path only. §6 is
correct.

One observation that belongs in §6 because it explains C1/C2: on a **default** run the predicate
executes 2–4 times per record — sanity check (record 1), adapter detection (≤1M), poly-G scan (≤1M,
when `-a` was given), then the loop. Still free, and it is the reason the refusal usually surfaces
before any writer opens, which is exactly what makes the untested residue (past 1M records, or with
`-a --no_poly_g`) the dangerous one.

---

## 5. Alternatives

1. **Add the cleanup, don't move the check** (addresses C1). Moving the refusal to sanity-check time
   only inspects record 1, so it cannot close the hole. Delete the created outputs on error, or write
   to a temp name and rename on success. The latter is the more general fix and would retire this
   whole class of partial-output bug (cf. #409), but it is a bigger change than #415 warrants —
   targeted cleanup on the refusal path, plus a V row, is the proportionate move.
2. **Keep the class split in reserve exactly as §2.3 says.** I found no evidence against A3, and one
   piece *for* A3b: `samtools import` puts only the first whitespace-delimited token in QNAME in all
   three modes I tried (default, `-i`, `-T '*'` — `@SRR123.1 1:N:0:ATCGATCG` → QNAME `SRR123.1`)
   *(ran it)*. The canonical FASTQ→uBAM converter therefore **cannot** produce a header-in-QNAME
   uBAM; such a file has to come from a hand-rolled converter. That is real support for "rare
   conversion mistake". I could not test Picard `FastqToSam` locally, so A3b stays a judgement — but
   it is now a judgement with one converter's behaviour behind it rather than none.
3. **For I1, commit the Python generator.** noodles refuses and samtools cannot express it, so a
   stdlib generator is the only route to the LF fixture — and it is 25 lines, which is cheaper than
   the documented coverage gap it removes.

---

## 6. Action items

**Critical**

1. **C1** — record A7 as **false**, with the evidence: a refused run leaves a partial output on every
   path, and the partial uBAM is valid, EOF-terminated and samtools-readable. Add a cleanup step (or
   temp-name-then-rename) to §5, and rewrite §11's risk 2 and §10's O2 accordingly.
2. **C2** — redesign the integration fixtures and invocations so the cases reach the code they name:
   clean first record/pair, whitespace later, `-a`/`-a2` + `--no_poly_g`, plus a `--cores 2` case for
   `:380`. Keep one record-1 case as the sanity-check arm. Note in §2.4 that four entry points
   (sanity check, adapter detection, poly-G scan, trimming loop) feed the three call sites.
3. **C3** — swap the `:261` / `:380` labels in §2.4 and fix the citation in §5 step 4.

**Important**

4. **I1** — drop the "not constructible" claim; add the stdlib BGZF generator recipe to
   `test_files/README.md` and an end-to-end LF case (corruption pre-fix, refusal post-fix), which
   retires §5 step 3's stated gap and §11's caveat about §2.3 being unit-tested only.
5. **I2** — add `append_to_id`'s `trim_end()` (`fastq.rs:175`, Unicode whitespace) to §2.7's table
   and soften A2's "exactly", with the trailing-VT-under-`--rename` exception stated as accepted
   collateral.
6. **I3** — replace "harmless" for `\v` / `@` / non-ASCII with the measured behaviour (uBAM output
   rejects them today with an opaque noodles error plus a partial file), and note the residual class
   as a follow-up rather than in scope.
7. **I4** — add a `CHANGELOG.md` entry to §5's deliverables and §7's file list.
8. **I5** — either offer the space-case remediation one-liner in the error text with its caveat, or
   narrow §5 step 1's claim to "no one-liner covers tab or newline names".

**Optional**

9. **O-a** — fold in the `--clump_only` lossless-contract violation and the multi-input partial-run
   shape; add a V row for the latter.
10. **O-b/O-c/O-d/O-e** — name `quickstart.md:55` as the docs target; correct the fixture count;
    reword V6 to match `fresh_tmpdir`; extend the escape set to `\r` and `\f`.

---

## Appendix — probe inventory

All under `$TMPDIR/revB`, all read-only with respect to the repo.

| Probe | Purpose | Result |
|---|---|---|
| `samtools view -b` from SAM, space QNAME | A4 | builds, reads back intact |
| `--output-format ubam --preserve-tags CB` | §2.1 row 1 | exit 0, QNAME `name`, tag gone, notice fires |
| `--preserve-tags CB` (FASTQ out) | §2.1 row 2 | `@name with space\tCB:Z:AAACCC` |
| aligned record 1, uBAM out | entry-point order | dies before adapter detection, nothing written |
| aligned record 2, `-a` only | entry-point order | dies in poly-G scan, nothing written |
| aligned record 2, `-a --no_poly_g` × {uBAM, FASTQ, `--cores 2`, paired × 2} | **A7** | partial output every time (C1) |
| aligned record 1, `--paired` interleaved | which caller reports it | sanity check, not `:261` (C2) |
| hand-rolled BGZF, LF QNAME (1 and 3 records) | §5 step 3 | constructible; 13 lines / 3 records, exit 0 (I1) |
| hand-rolled BGZF, LF QNAME, uBAM out | pre-fix behaviour | opaque noodles error + partial BAM |
| hand-rolled BGZF: `\v`, `@`, non-ASCII, control | §2.7 "harmless" | first three fail on uBAM out (I3) |
| hand-rolled BGZF, tab QNAME → SAM round-trip | I5 caveat | `Parse error at line 3` |
| `samtools view -h \| awk FS=OFS="\t" \| samtools view -b` | I5 | works; tag survives TrimGalore |
| `samtools import` × {default, `-i`, `-T '*'`} | A3b | QNAME is the first token only |
| `--clump_only --output-format ubam` | §5 step 5 | tag dropped silently, exit 0 (O-a) |
| two SE inputs, second bad | O-a | input 1 complete, input 2 partial |
| `rustc` probe of `is_ascii_whitespace` / `is_whitespace` / `trim_end` / `split_ascii_whitespace` | §2.7, I2 | table correct; `trim_end` strips VT/NBSP/U+3000 |
