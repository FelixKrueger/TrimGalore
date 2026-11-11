# Plan Coverage Report

**Mode:** B (code vs. plan — the design plan doubles as the implementation spec; there is no `IMPL.md`)
**Plan:** `plans/08132026_partial-output-on-refusal/PLAN.md` (v2, 2026-08-14)
**Code audited:** commit `0a6f3a3` on `fix/428-no-partial-output-on-refusal` (parent `285afba`)
**Date:** 2026-08-14
**Verdict:** INCOMPLETE — 6 items unresolved (5 substantive, 1 documentation nit), plus T9, an
unrun CI gate rather than a gap

## Working-tree drift after this audit (recorded 2026-08-14, late in the run)

**This report audits commit `0a6f3a3`.** While it was being written, five files were modified in
the working tree, uncommitted. Everything below describes the commit, not the current tree; the
uncommitted changes are **not** covered by any verification in this report, and the 701-passing
run and the clean fmt/clippy results were measured before the drift was complete.

`git diff --stat HEAD` at the time of writing: `docs/.../guide/outputs.md` (6), `src/io.rs` (30),
`src/fastq.rs` (20), `src/demux.rs` (5), `src/bam.rs` (4). No test file changed.

Read against the gaps below, that drift already closes three of them:

- **Gap 3** — `outputs.md:145`'s "a mid-stream refusal would leave a partial BAM" is rewritten to
  the `main.rs:399-403` wording.
- **Gap 6** — "a full disk" is dropped from the "Refused runs" paragraph, and the leftover-
  `.partial` sentence now covers the failed-rename case as well as the killed-process one.
- **Gap 1** — `partial_output_name` gains the exact note the plan asked for ("No suffix may be
  appended to a name that has been through this"), plus the injectivity statement the tripwire
  rests on.

It also carries one change that is **not** a documentation edit and that nothing in this report
has exercised: `partial_output_name` now builds the temporary as an `OsString` rather than via
`to_string_lossy`, so the mapping is byte-exact on non-UTF-8 filenames. That is what makes the
new injectivity claim true, and it wants a test run of its own.

Gaps 2, 4 and 5 are untouched by the drift.

## Summary

- Total items: 42
- DONE: 31
- PARTIAL: 9
- MISSING: 0
- DEVIATED: 2

Nothing the plan specified is absent wholesale, and there is no MISSING entry. The 9 PARTIALs
collapse to **6 distinct findings**, because §5 step 6's unkept promise surfaces three times in
the ledger (items 6, 27, 37) and T7 twice (items 18, 30):

1. a behavioural promise the code does not make — a failed commit does not name what already
   committed (§5 step 6, §4's edge-case table);
2. a **fifth** surviving statement of the reversed contract, in `docs/.../guide/outputs.md`;
3. one §2.4 arm (`--dont_gzip`) neither tested nor waived;
4. T7's confirmation was never performed;
5. a user-facing docs sentence claiming the stronger guarantee §4 item 4 explicitly declined;
6. the `partial_output_name` invariant note the plan asked for, in one sentence (nit).

Plus T9, which cannot be evaluated: the branch is unpushed, so CI has not run.

The 2 DEVIATED items carry §12's deviations 1, 2 and 3. Deviation 4 (T2's decoder choice) is a
substitution the plan itself offered, so item 13 is DONE rather than DEVIATED.

## How the ledger was built

- **§5 Implementation outline**, steps 1–11 → items 1–11.
- **§7 Validation**, T1–T10 → items 12–21.
- **§4 Behaviour** items 1–8 → items 22–29; the edge-case table's ten rows → items 30–39.
- **§8 Assumptions** → items 40–42 (only those that constrain the code: A6, A7, A11; the rest
  are statements of fact about the pre-change tree, not work items).

Every claim below was checked against the working tree at `0a6f3a3` and, where the plan
specified a test, against a full `cargo test --release` run. The tree was not modified: `git
status` shows no tracked changes.

## Coverage ledger

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `PendingOutput` type in `io.rs`: `create`/`commit`/`tmp_path`/`Drop`, `released` set **before** the rename, `partial_output_name` dot-prefixed and appended to the whole file name | §5 step 1 | DONE | `src/io.rs:669-751`. `commit` sets `self.released = true` at `:729`, the line *before* `fs::rename` at `:730`, so a failed rename leaves the temporary. `Drop` at `:745-751`. Sub-clause not literally present — see Gaps §1 |
| 2 | `FastqWriter` holds the pending output and gains `finish()`; `writer` declared before `pending` | §5 step 2 | DEVIATED | §12 dev. 1 (no `Option<Box<dyn Write>>`) and dev. 2 (`flush()` removed). Both documented and both justified — see Deviations |
| 3 | `BamWriter`: same shape; `try_finish()`, drop `inner`, then commit | §5 step 3 | DONE | `src/bam.rs:525-534`, `:654-666`. `try_finish` → `pending.take()` → `drop(self)` → `commit()` |
| 4 | The raw `File::create` record-output sites converted, and the commit-point invariant honoured at each | §5 step 4 | DONE | All four named sites converted; every in-run consumer verified to run *after* the commit (see Test verification → ordering) |
| 5 | Every driver calls `finish()` — no writer relies on `Drop` | §5 step 5 | DONE | 18 `FastqWriter` sites (6 `specialty.rs`, 2 `demux.rs`, 10 `main.rs`) + 7 `BamWriter` sites, each with a matching `finish()`. No `File::create` for record output survives in `src/` outside `#[cfg(test)]` |
| 6 | Multi-output commit ordering: fixed order, **nothing fallible between the first commit and the last**, and on mid-sequence failure "error out naming what committed and what did not" | §5 step 6 | PARTIAL | Order and the no-fallible-work window hold in **all seven** multi-output drivers. The error message does **not** name what committed — see Gaps §2 |
| 7 | `main.rs:395-405` — the `--rename`/uBAM refusal keeps its refusal but loses the clause asserting the defect | §5 step 7 | DONE | `src/main.rs:399-403` now reads "because a per-record decision would keep the annotation for some reads and silently drop it for others". The **docs twin of the same sentence was not updated** — counted under item 9 |
| 8 | The residue-asserting test is inverted | §5 step 8 | DONE | See T1 (item 12) |
| 9 | The four places that document the reversed contract | §5 step 9 | PARTIAL | All four updated (`src/cli.rs:294-296` `--help`, `docs/.../modes/passthrough.md:71`, `src/parallel.rs:90`, `src/parallel.rs:271`). A **fifth** statement survives — see Gaps §3 |
| 10 | Tripwire header: `.partial` in "Known-unplanned writes", injectivity argument, plus "a red `.partial` means a missed commit, not a missing candidate" | §5 step 10 | DONE | `tests/integration_preflight_tripwire.rs:41-49`, inside the "Known-unplanned writes" section (`:23`) and before "## Limits" (`:51`), with the exact clause in bold. `snapshot()` uses `read_dir`, which lists dotfiles, so the claim that a red `.partial` is detectable is sound |
| 11 | Docs sentence in `guide/outputs.md` + one `CHANGELOG.md` sentence | §5 step 11 | DONE | New "## Refused runs" section at `docs/src/content/docs/guide/outputs.md:122-126`; CHANGELOG entry is one sentence in the file's existing format. One clause in it overstates — see item 25 |
| 12 | **T1** — the pinned residue arm now asserts absence | §7 T1 | DONE | `tests/integration_ubam.rs:306` `ws_qname_refused_in_trimming_loop_ubam_out_leaves_no_output`. Keeps the "BAM record 2" stderr assertion; the `BGZF_EOF` block is replaced by absence of the final name **and** of any `.partial`. `ws_qname_refused_in_trimming_loop_fastq_out` (`:345`) also gained an absence assertion at `:366-370` |
| 13 | **T2** — the trailer is not lost to the rename, on three arms, at decoder level | §7 T2 | DONE | Three arms present and each asserts decoder-level integrity, not a bare count — see Test verification. §12 dev. 4 documented |
| 14 | **T3** — one case per §2.4 arm, asserting both halves (no final name, no temporary) | §7 T3 | PARTIAL | 11 of §2.4's 12 rows have a case; `--dont_gzip` has neither a case nor a waiver — see Gaps §4. Both halves are asserted everywhere, via `assert_nothing_published` |
| 15 | **T4** — a previous good output survives a failed re-run, **same output path**, reports asserted too | §7 T4 | DONE | `failed_rerun_leaves_the_previous_good_output_intact` copies the fixture to `subject.bam` for both runs, so both target `subject_trimmed.fq`, and compares the output **and** `subject.bam_trimming_report.txt` byte-for-byte |
| 16 | **T5** — `--cores 2`'s zero-byte residue is gone (absent, not empty) | §7 T5 | DONE | `refused_parallel_publishes_nothing_not_an_empty_file`; `assert_nothing_published` tests `!exists()`, and reports the byte count in the failure message |
| 17 | **T6** — the paired arm, and the >4096-pair large-plausible-residue regime | §7 T6 | DEVIATED | Both paired arms assert absence (`ubam_ws_qname_paired.bam`, FASTQ and uBAM out). The >4096-**pair** case became a 5000-**record** single-end case, §12 dev. 3 — see Deviations |
| 18 | **T7** — legitimate empty output still appears (confirm an existing all-reads-filtered test asserts the file's *presence*) | §7 T7 | PARTIAL | The confirmation was not recorded anywhere, and no all-reads-filtered test exists to confirm — the nearest coverage is **zero-record-input**, §4's other edge-case row — see Gaps §5 |
| 19 | **T8** — no suite regression; do not gate on a hard-coded number | §7 T8 | DONE | Verified, not trusted: `cargo test --release` → **701 passed, 0 failed** across 16 binaries (436 lib + 265 integration). Matches §12's 701 and the 684 + 17 arithmetic. T8's own platform arithmetic is off by one — Linux is 699, not 700 — see Test verification |
| 20 | **T9** — `Validate vs Perl TrimGalore` + `Validate uBAM input` green, by name | §7 T9 | PARTIAL | Both jobs exist by exactly those names (`.github/workflows/ci.yml:313`, `:957`), but the branch is **not pushed**, so neither has run against `0a6f3a3`. Not a code gap; an unrun gate |
| 21 | **T10** — the FastQC artefact is named from the *published* output | §7 T10 | DONE | `fastqc_report_is_named_from_the_published_output` asserts `ubam_test_with_tags_trimmed_fastqc.html` is present, and additionally that no surviving name starts with `.` or contains `.partial` |
| 22 | §4.1 — drivers ask by the **final** path; naming, pre-flight and report naming unchanged | §4 item 1 | DONE | `PendingOutput::create(final_path)` and both writers' `create` keep their existing signatures; no candidate list or namer changed. The tripwire's 55 cases pass unchanged |
| 23 | §4.2 — bytes go to `.{file_name}.partial`, created truncating | §4 item 2 | DONE | `src/io.rs:678-687` + `:715` (`File::create`). No test pins the *name shape*, so O2's glob-invisibility decision rests on code reading — noted, not a gap (the plan specified no test for it) |
| 24 | §4.3 — commit at the point where the driver closes the writer today | §4 item 3 | DONE | Verified against every in-run consumer, in both `clump_only.rs` and `main.rs` — see Test verification → ordering |
| 25 | §4.4 — the guarantee is "the final name appears only after the writer has been closed", **not** that the bytes are complete | §4 item 4 | PARTIAL | The code and `FastqWriter::finish`'s doc comment state exactly the weaker guarantee, and per-sink tracing shows only serial gzip's trailer can publish on a teardown failure. The shipped user-facing sentence states the stronger guarantee for one of its three examples — see Gaps §6 |
| 26 | §4.5 — on any error/early return/unwinding panic the temporary is removed, the final path is never created, and a pre-existing file there is untouched | §4 item 5 | DONE | `Drop for PendingOutput` (`src/io.rs:745-751`); asserted by 14 refusal cases plus T4 for the untouched-pre-existing half |
| 27 | §4.6 — multi-output drivers commit after `Ok`, fixed order, nothing fallible between | §4 item 6 | PARTIAL | Same item as ledger #6 |
| 28 | §4.7 — multi-input runs unchanged (failure on input 3 of 5 leaves 1–2 committed) | §4 item 7 | DONE | The per-input and "Pair N of M" loops are untouched; commits happen inside each per-input driver, so the semantics follow by construction. Multi-pair suite cases pass |
| 29 | §4.8 — the temporary is not user-facing; named in error text only if the rename fails | §4 item 8 | DONE | `commit`'s only error context names both paths: `"Failed to publish {final} (the trimmed data is intact at {tmp})"`. `tmp_path()` exists as the plan specified but is never called |
| 30 | Edge: zero records survive filtering → commit runs, a valid empty output appears | §4 table | PARTIAL | Same item as T7 (#18). The behaviour is right by construction (`finish()` runs unconditionally after the loop) and is exercised where the emptiness comes from routing (`integration_output_collision.rs:806`, which requires the empty files to exist). The two tests that reach it **by filtering** assert only `created ⊆ planned` — see Gaps §5 |
| 31 | Edge: zero-record input → empty `.fq` / header-only BAM appears | §4 table | DONE | Covered by existing, still-green tests: `src/parallel.rs:1449 test_parallel_empty_input_no_deadlock` opens the output through `FastqReader::open(&output_path)?` (which fails if the file is absent) and asserts zero records; `src/clump_only.rs:1399 test_clump_only_empty_input` and `:1531 test_empty_bam_input_produces_valid_bam` do the clump and BAM halves, the latter with an explicit `exists()`. Note this row is only reachable **in-process**: `src/format.rs:544` rejects an empty regular FASTQ with "is empty" before any writer exists, so the CLI never takes it |
| 32 | Edge: a user file already at the temporary path → overwritten | §4 table | DONE | By construction (`File::create` truncates). No test; none was specified |
| 33 | Edge: a stale temporary from a SIGKILLed run → overwritten by the next run, never renamed by anyone else | §4 table | DONE | By construction; documented for users at `guide/outputs.md:126` ("safe to delete; the next run at the same output path overwrites it") |
| 34 | Edge: abort during unwinding → no `Drop`, a `.partial` is left, the final name is still never created | §4 table | DONE | Rests on code reading — the final name is only ever created by `commit`, which is only reached from `finish`. A `ParCompress` teardown panic on the *ordinary* path is better than the plan predicted: it unwinds through `drop(self)` while the taken `pending` is still a live local, so its `Drop` removes the temporary. Not testable without inducing a double panic; the plan specified no test |
| 35 | Edge: `--output_dir` on a different mount → irrelevant | §4 table | DONE | `partial_output_name` joins the temporary onto `final_path.parent()`, so it is a sibling by construction; no `EXDEV` path exists |
| 36 | Edge: rename fails → hard error naming both paths, temporary **left in place** | §4 table | DONE | The `released`-before-rename ordering makes this work (`src/io.rs:729-730`); the message names both paths. No test (would need a read-only directory); none was specified |
| 37 | Edge: commit #2 of N fails after #1 succeeded → hard error naming what committed and what did not | §4 table | PARTIAL | Same item as ledger #6 — the naming half is not implemented |
| 38 | Edge: two concurrent runs; and B is better for a concurrent *reader* | §4 table | DONE | By construction — the previous inode survives until the rename. No behaviour claimed beyond today's |
| 39 | Edge: `--fastqc`, `--demux`, `clump_only`'s `output_bytes` consume the output by its **final** path in the same run | §4 table | DONE | This is the §5 step 4 invariant; verified at all ten consumer sites, and pinned by T10 and `refused_demux_publishes_nothing` |
| 40 | A6 — `FastqWriter` is only ever constructed from a path | §8.6 | DONE | `create` remains the only constructor; no `new`/stdout sink exists |
| 41 | A7 — `writer` declared before `pending` in both writers | §8.7 | DONE | `src/fastq.rs:574-581` and `src/bam.rs:525-534`, each with a one-line comment stating why |
| 42 | A11 — no new runtime dependency | §8.11 | DONE | `Cargo.toml` is not in the diff; `tempfile` stays a dev-dependency and is not used by the new code |

## Deviations (all four from §12), and whether the documentation justifies them

**Dev. 1 — `Option<Box<dyn Write + Send>>` was not needed** (ledger #2). Justified, and
verifiable in the code: `finish(mut self)` does `self.writer.flush()?` → `self.pending.take()`
→ `drop(self)` → `commit()`. Because `FastqWriter` has a `Drop` impl (pre-existing, flushing
the box) and its fields drop in declaration order, `drop(self)` closes the compressor and emits
the trailer *before* `commit()` renames. That is the same ordering guarantee step 2 asked for,
with no `expect`/panic path added to `write_record`. `BamWriter::finish` is the identical shape.
Documentation is accurate.

**Dev. 2 — `FastqWriter::flush()` removed, because `finish()` consuming `self` does *not*
enumerate the call sites** (ledger #2). This corrects a false claim in §5 step 5, and it is the
more load-bearing of the four: `writer.flush()?; drop(writer);` would have compiled and
silently published nothing. Removing `flush()` is what forced all 18 sites to surface. Verified
independently: no `FastqWriter::flush` remains, and every creation site has a matching
`finish()`. Well documented, and the residual risk is now genuinely closed rather than asserted
— the plan's §11 risk (b) ("~20 mechanical sites, where a missed one is a silent hole") is
retired by the compiler rather than by care.

**Dev. 3 — T6's >4096-pair paired case became a 5000-record single-end case** (ledger #17).
Documented with two reasons, and both check out: `examples/mk_bam_fixture.rs::record()` does
hardcode `FLAG_UNMAPPED` (`:20`, used unconditionally), so mate-adjacent paired generation
needs a generator change; and `BAM_READER_BATCH_SIZE = 4096` batching caps a paired residue at
one batch, so 5000 single-end records is the strictly larger plausible residue. What the
substitution does **not** deliver is the property T6 named it for by construction: no test now
exercises the paired de-interleaver past its first batch. The plan's §11 risk (d) anticipated
this case being dropped for runtime cost; it was relocated rather than dropped, and the
relocation preserves "large plausible residue" while losing "paired, multi-batch". Adequate,
with that one loss stated rather than hidden.

**Dev. 4 — T2's decoder-level integrity uses `FastqReader` rather than `gzip -t`** (ledger
#13). Justified and verified: `FastqReader::open` builds a `MultiGzDecoder`, and
`next_record()` propagates the read error through `?` (`src/fastq.rs`, the `ReaderSource::Direct`
arm), so a stream whose CRC/ISIZE trailer is missing fails the `.expect("the output must decode
cleanly to EOF")` in `read_back`. This is exactly the "`MultiGzDecoder` read-to-end that
propagates the trailer error" alternative T2 itself offered, and it avoids an external test
dependency.

## Gaps (detail)

### Gap 1 — the `partial_output_name` invariant note is not the note the plan asked for (ledger #1, low)

**Expected:** §5 step 1: *"Put the 'no suffix may be appended to a committed name' note above
that function."*
**Found:** `src/io.rs:671-677` documents the *rationale* — the marker is appended to the whole
file name and dot-prefixed, "so it neither extends a recognised extension nor matches a `*`
glob" — with the `publishDir` reasoning from O2. The constraint on future namers is not stated
as an invariant.
**Gap:** one sentence above `partial_output_name`. The plan's own wording for it is ambiguous,
and the substance (why the marker goes on the whole name) is present, so this is a
documentation nit rather than a missing guard.

### Gap 2 — a mid-sequence commit failure does not name what committed (ledger #6, #27, #37)

**Expected:** §5 step 6: *"On a mid-sequence failure, error out naming what committed and what
did not."* §4's edge-case table repeats it: *"Commit #2 of N fails after #1 succeeded | hard
error naming what committed and what did not."*
**Found:** every multi-output driver commits in a fixed order with nothing fallible in between
— `main.rs::run_paired` (R1, R2, passthrough, unpaired R1, unpaired R2, `:1629-1639`),
`main.rs::run_paired_ubam_single_file` (`:1895-1902`), `parallel.rs::run_paired_end_parallel`
(handles dropped, then `for guard in pending_outputs { guard.commit()? }`, `:310-314`),
`demux.rs::demultiplex` (barcode-sorted, then NoCode, `:256-261`),
`clump_only.rs::clump_only_paired` (`:520-521`), `specialty.rs::clock` (`:332-333`) and
`specialty.rs::implicon` (`:427-428`). Three of the seven carry the "nothing fallible in
between" comment. But the error a failed commit produces is `PendingOutput::commit`'s own
context — `"Failed to publish <final> (the trimmed data is intact at <tmp>)"` — which names the
file that failed and its temporary, and says nothing about the siblings that already committed.
No driver wraps the commit sequence to add that.
**Gap:** either driver-level context that names the already-committed set (e.g. commit through
a helper that accumulates the published paths and, on failure, reports them), or an explicit
plan amendment saying the per-file message is enough. As it stands the plan makes a promise the
binary does not keep, and it is the promise that covers the one non-atomic window the design
knowingly leaves open.

### Gap 3 — a fifth statement of the reversed contract survives, in the same file step 11 edited (ledger #7, #9)

**Expected:** §5 step 9 enumerates four locations; the audit brief asks that no *fifth*
statement of the old contract survive in `src/`, `tests/` or `docs/` (the archived
`docs/.../reference/changelog.md` excluded).
**Found:** `docs/src/content/docs/guide/outputs.md:145`, in the "Annotating read IDs" section:

> "…so the whole run is refused rather than decided per record — **a mid-stream refusal would
> leave a partial BAM**. uBAM input is accepted…"

This is the docs mirror of `main.rs:399-403`, the exact sentence step 7 rewrote in the source.
It sits **20 lines below** the new "## Refused runs" section added by step 11 in the same file,
which says a refused run leaves nothing — so the page now contradicts itself.
**Why the sweep missed it:** §12 records the old-claim sweep as repo-wide over "partial output".
This line says "partial **BAM**", so a `partial output` grep does not match it. The three
surviving matches §12 clears are all genuine (`integration_non_restartable_input.rs:177/206`
about a *hang*, `integration_paired_format_guard.rs:13` about the multi-pair pre-flight); the
sweep's term, not its scope, is what let this through.
**Gap:** rewrite that clause the way `main.rs:399-403` was rewritten — the substantive reason
(representability depends on the individual header) already precedes it in the same sentence.

### Gap 4 — §2.4's `--dont_gzip` row is neither tested nor waived (ledger #14)

**Expected:** T3: *"Every measured arm in §2.4 leaves nothing … one case per arm."*
**Found:** §2.4 has twelve rows. Eleven are covered:

| §2.4 row | Case |
|---|---|
| `--output-format ubam` | `refused_se_ubam_publishes_nothing` (+ T1 in `integration_ubam.rs`) |
| *(default FASTQ out)* | `refused_se_fastq_publishes_nothing` |
| `--dont_gzip` | **none** |
| `--cores 2` | `refused_parallel_publishes_nothing_not_an_empty_file` |
| `--clumpify --cores 2` | `refused_clumpify_publishes_nothing` |
| `--paired` (interleaved) | `refused_paired_fastq_publishes_nothing` |
| `--paired --output-format ubam` | `refused_paired_ubam_publishes_nothing` |
| `--hardtrim5 10` | `refused_hardtrim5_publishes_nothing` |
| `--hardtrim3 10` | `refused_hardtrim3_publishes_nothing` |
| `--hardtrim5 10 --output-format ubam` | `refused_hardtrim5_ubam_publishes_nothing` |
| `--clump_only --output-format ubam` | `refused_clump_only_publishes_nothing` |
| `--demux` | `refused_demux_publishes_nothing` |

`--dont_gzip` appears in neither T3's own how-column enumeration nor the "Consciously
uncovered" note (`--clock`, `--implicon`, `--passthrough`, `--retain_unpaired`), so it is the
one row that is neither tested nor waived. Note the plan is internally inconsistent here: T3's
*what* column says "every measured arm", its *how* column lists eight.
**Gap:** either one more case (`refuse(&dir, &["--dont_gzip"], OFFENDER)`, three lines) or an
added line to the "Consciously uncovered" note. Materially this arm is the same code path as
the covered default-FASTQ arm — for a `.bam` input the default output is already plain `.fq`,
so `--dont_gzip` reaches the same `gzip == false` branch of `FastqWriter::create` — but the
plan drew the boundary at "every measured arm", and this one is inside it.

### Gap 5 — T7's confirmation was not performed (ledger #18, #30)

**Expected:** T7: *"confirm at least one existing all-reads-filtered test asserts the file's
**presence**"* — the guard against the fix overshooting into "leaves no file when the output is
legitimately empty".
**Found:** §12 does not mention T7 at all: no confirmation is recorded in the plan, the commit
message, or the new test file. Searching for it independently: the nearest existing coverage is
`src/parallel.rs::test_parallel_empty_input_no_deadlock` (`:1449`), which is a **zero-record
input** — §4's *other* edge-case row — and does assert presence, since
`FastqReader::open(&output_path)?` fails if the file is missing. `clump_only.rs`'s
`test_empty_bam_input_produces_valid_bam` covers the BAM half of the same row.

**Two tests do reach the all-reads-filtered state, and neither asserts presence** —
`tests/integration_preflight_tripwire.rs:660` `paired_trim_retain_unpaired` and `:675`
`paired_trim_retain_unpaired_cores`. Their fixture helper `write_fastq_short` (`:134`) writes
40 records of a 5 bp read for R2 against an 80 bp R1, so under the default `--length 20` every
pair is discarded and `s_R1_val_1.fq` / `s_R2_val_2.fq` are zero-record files. But `Tripwire`
asserts `created ⊆ planned` plus `!created.is_empty()`, and the comment at `:436-438` says the
reverse relation is false by design — so if `s_R2_val_2.fq` vanished entirely, both tests would
still pass. This is the precise shape T7 was written to guard against, sitting in the suite
already, one assertion short.

One existing test **does** require a zero-record output's presence, though the zero comes from
routing rather than filtering: `tests/integration_output_collision.rs:806`
`demux_writes_exactly_the_planned_paths` compares planned paths against a `read_dir` listing with
`assert_eq!` in **both** directions (`:843`), and its barcodes are chosen so that
`m_trimmed_sample2.fq` and `m_trimmed_NoCode.fq` are legitimately empty. That does exercise
"commit a legitimately empty output" on the widest commit set in the change — `demux.rs`'s
per-barcode fan-out, where every writer is created and finished unconditionally.

The overshoot guard is therefore covered in substance on two paths (zero-record input, and
demux's empty-by-routing outputs); what is missing is a filter-driven case and the confirmation
T7 asked for.
**Gap:** either point at an all-reads-filtered test that does assert presence (and record the
confirmation), or accept the zero-record-input tests as discharging T7 and say so, or add the
case. The substance of the overshoot guard is partly covered — an empty output *is* asserted
present on the zero-record path, and that path exercises the same commit-on-`Ok` code — but the
specific confirmation T7 asked for was not made.

### Gap 6 — the shipped docs sentence claims the stronger guarantee for disk-full (ledger #25)

**Expected:** §4 item 4 states the guarantee deliberately narrowly — "the final name appears
only after the writer has been closed. **Not** 'the bytes are complete'" — because every
teardown swallows its error, so a disk-full at trailer time still publishes a trailerless `.gz`
or an EOF-less BAM with exit 0. §9 flags exactly that as the one candidate clearing the filing
bar, and O5 defers the fix.
**Found:** `docs/src/content/docs/guide/outputs.md:124`: "a refusal — a malformed read name a
million records in, a truncated input, **a full disk** — leaves nothing behind".

Traced per sink, the claim is *nearly* true, and the exception is exactly the one the plan
carved out:

- **plain FASTQ** — `finish()`'s `self.writer.flush()?` propagates, so a full disk during or at
  the end of writing publishes nothing. Covered.
- **uBAM** — `BamWriter::finish` calls `try_finish()` explicitly and propagates before the
  commit. Covered.
- **`--cores N` gzip** — `gzp::ParCompress::drop`'s `.unwrap()` panics on a teardown failure.
  That panic unwinds through `drop(self)` inside `finish()`, at which point the taken `pending`
  is a live local, so its `Drop` removes the temporary: nothing is published. Covered, by
  panicking.
- **serial gzip** — `BufWriter::flush` reaches `GzEncoder::flush`, which does not write the
  trailer. So a disk-full in the last handful of bytes (deflate end + the 8-byte CRC/ISIZE) is
  swallowed by `GzEncoder::drop`'s `let _ = try_finish()`, `finish()` returns `Ok`, and a
  trailerless `.gz` is published with exit 0. **Not** covered.

**Gap:** narrow the example or drop "a full disk", so the page claims the guarantee §4 item 4
delivers ("closed", not "complete"). The window is one sink and a few bytes wide, but it is the
precise case §4 item 4 and §9's filing candidate were written to describe, and the docs are the
only place the stronger claim is made — the code and `FastqWriter::finish`'s doc comment both
state it correctly.

## Test verification

`cargo test --release` at `0a6f3a3`: **701 passed, 0 failed, 0 ignored**, 16 test binaries
(436 lib + 265 integration + 0 bin + 0 doc). Matches §12 and the plan's 684 + 17.

**The platform arithmetic in T8 is wrong, in the plan and in the brief given to this audit, and
the correct figure is worth recording** since a future reader will use it as the baseline. There
are **three** `#[cfg(target_os = "macos")]` tests, not two — `src/format.rs:689`, `:711` **and
`tests/integration_non_restartable_input.rs:252`** (`passthrough_dev_stdin_is_rejected`) — plus
the one `#[cfg(target_os = "linux")]` test at `src/cli.rs:1950`. So Linux reports **699**, two
fewer than macOS, not one. This does not affect T8's verdict: the plan's instruction was not to
gate on a hard-coded number, and the run is green either way.

Gates: `cargo fmt --all -- --check` clean; `cargo clippy --all-targets --release -- -D warnings`
reports no diagnostics (a cache hit at this source fingerprint, i.e. a replayed success for
these exact sources, not a fresh compile).

### New and changed tests

| Test | File | Plan item | Status |
|---|---|---|---|
| `refused_se_fastq_publishes_nothing` | `tests/integration_partial_output.rs:100` | T3 | PASS |
| `refused_se_ubam_publishes_nothing` | `:107` | T3 | PASS |
| `refused_parallel_publishes_nothing_not_an_empty_file` | `:116` | T3, T5 | PASS |
| `refused_clumpify_publishes_nothing` | `:125` | T3 (§2.4 tenth arm) | PASS |
| `refused_hardtrim5_publishes_nothing` | `:132` | T3 | PASS |
| `refused_hardtrim3_publishes_nothing` | `:139` | T3 | PASS |
| `refused_hardtrim5_ubam_publishes_nothing` | `:146` | T3 | PASS |
| `refused_clump_only_publishes_nothing` | `:157` | T3 | PASS |
| `refused_demux_publishes_nothing` | `:180` | T3 (the fan-out arm) | PASS |
| `refused_paired_fastq_publishes_nothing` | `:195` | T6 | PASS |
| `refused_paired_ubam_publishes_nothing` | `:203` | T6 | PASS |
| `refused_after_thousands_of_good_records_publishes_nothing` | `:217` | T6 (deviated: SE, 5000 records) | PASS |
| `failed_rerun_leaves_the_previous_good_output_intact` | `:239` | T4 | PASS |
| `published_serial_gzip_output_decodes_to_eof` | `:309` | T2(a) | PASS |
| `published_parallel_gzip_output_decodes_to_eof` | `:328` | T2(c) | PASS |
| `published_ubam_output_carries_the_bgzf_eof_marker` | `:349` | T2(b) | PASS |
| `fastqc_report_is_named_from_the_published_output` | `:382` | T10 | PASS |
| `ws_qname_refused_in_trimming_loop_ubam_out_leaves_no_output` | `tests/integration_ubam.rs:306` | T1 (inverted) | PASS |
| `ws_qname_refused_in_trimming_loop_fastq_out` | `tests/integration_ubam.rs:345` | T1 (absence assertion added) | PASS |
| all-reads-filtered presence guard | — | T7 | NOT FOUND (Gap 5) |
| `--dont_gzip` refusal | — | T3 | NOT FOUND (Gap 4) |

**One case does not assert *why* the run failed**, and was checked by hand rather than left in
doubt: `refused_clump_only_publishes_nothing` bypasses the `refuse()` helper and asserts only
`!status.success()`, so a CLI-validation rejection would satisfy it vacuously. Running the
binary directly reproduces the intended path — `Error: test_files/ubam_ws_qname_late.bam: BAM
record 2 … BAM read name contains whitespace` — with an empty output directory afterwards. The
case is sound; T3 did not require the stderr assertion, so this is an observation, not a gap.

**T2's three arms are genuinely three mechanisms**, and arm (c) really does reach
`gzp::ParCompress`: `main.rs:474-482` passes `cli.cores` into `specialty::hardtrim5`, which
forwards it to `FastqWriter::create`, whose `gzip && cores > 1` branch builds
`ParCompressBuilder::<Gzip>` (`src/fastq.rs:609-623`). The fixture is a `.gz` input, so `gzip`
is true and the asserted output is `BS-seq_10K_R1.30bp_5prime.fq.gz`. Integrity in (a) and (c)
is decoder-level via `read_back`, and in (b) is a byte-suffix check for the 28-byte BGZF EOF
marker.

### Commit-point ordering (§5 step 4's invariant), verified site by site

| Consumer | Commit | Consumer runs at |
|---|---|---|
| `clump_only.rs` SE FASTQ → `fastqc::run` | `pending_output.commit()?` `:349` | `:375` |
| `clump_only.rs` PE FASTQ → `fastqc::run` | `:520-521` | `:555-556` |
| `clump_only.rs` SE uBAM → `metadata`, then `fastqc::run` | `writer.finish()?` `:823` | `:828`, `:848` |
| `clump_only.rs` PE uBAM → `metadata`, then `fastqc::run` | `:1069` | `:1072`, `:1096` |
| `main.rs` SE → `fastqc::run`, `demux::demultiplex` | `writer.finish()?` `:1334` | `:1459`, `:1474` |
| `main.rs` PE serial → `fastqc::run` ×3 | `:1629-1639` | `:1775`, `:1781`, `:1793` |
| `main.rs` SE uBAM → `fastqc::run` | `:2154` | `:2240` |
| `main.rs` PE uBAM (interleaved in) → `fastqc::run` | `:2318` | `:2397` |
| `main.rs` PE uBAM (FASTQ pair in) → `fastqc::run` | `:2447` | `:2520` |
| `parallel.rs` SE/PE → the caller's `fastqc`/`demux` | inside `thread::scope`, before return (`:991`, `:310-314`) | after the call returns |

Every commit precedes every consumer of the same final path. The `?`-propagating
`metadata(&output_path)` calls at `clump_only.rs:828` and `:1072` — the ones a late commit
would have failed a *good* run on — both sit after their `finish()`.

### The four "consciously uncovered" arms

The plan waives `--clock`, `--implicon`, `--passthrough` and `--retain_unpaired` on the grounds
that each is the same code shape as a covered arm and the compiler is the guard. Checked, because
the risk this change carries on those arms is not a residue but the opposite — a missed
`finish()` publishing *nothing* on a **successful** run. Three of the four already have an
existing success-path test that asserts output presence and would have failed:
`specialty.rs:589-590` (`test_clock_happy_path_extracts_umi_and_clips`), `specialty.rs:720-721`
(`test_implicon_happy_path_umi_8`) and `tests/integration_passthrough.rs:135-145`
(`passthrough_cli_cores_1_smoke`, all three outputs). `--retain_unpaired`'s only local coverage
is in the tripwire, whose relation is `created ⊆ planned` plus a non-empty liveness check — a
single unpublished output would pass it — so that one rests on the compiler enumeration alone,
plus the CI job `Validate --retain_unpaired × multi-pair produces per-pair unpaired files`. The
waiver is sound as written; the suite is green.

### §2.5's site inventory, re-walked

All 34 non-test record-output creation sites listed in §2.5 now route through `PendingOutput`,
directly or via a writer that holds one: `main.rs` (10 `FastqWriter` + 3 `BamWriter`),
`specialty.rs` (6 + 2), `parallel.rs` (1 + 5 raw handles), `demux.rs` (2), `clump_only.rs`
(3 raw + 2 `BamWriter`). `grep File::create src/` leaves only report writers, `io.rs`'s own
`PendingOutput::create`, and `#[cfg(test)]` sites — matching the plan's "reports are out of
scope".

### The mutation-test claim in §12

Not re-run: it required editing `src/io.rs`, and this audit is read-only on a tree two other
reviewers are reading. The claim's arithmetic has a prose slip — "all 13 refusal cases, T4, and
both inverted `integration_ubam` cases" double-counts T4 (12 refusal cases + T4 + 2 inverted =
15, the number stated). The 4-that-stayed-green figure matches the four success-path publish
cases exactly, so the accounting is self-consistent once T4 is counted once.

## Verdict

**INCOMPLETE — 6 items unresolved (5 substantive, 1 nit).** Nothing was silently dropped: all
eleven implementation steps and all ten validations are present in some form, the 42-item ledger
has no MISSING entry, the 34-site inventory is fully converted, every commit precedes every
in-run consumer of its final path, and each of the four §12 deviations is documented with a
justification that holds up on inspection. What remains is five concrete pieces of work, one
documentation nit, and one gate that has not run.

To close it, in the order a reviewer would want them:

1. **Gap 3 — `docs/src/content/docs/guide/outputs.md:145`.** Rewrite "a mid-stream refusal
   would leave a partial BAM" the way `main.rs:399-403` was rewritten. This is the one item
   that ships a user-visible contradiction: the same page asserts both contracts, 20 lines
   apart, and the false half is the one step 7 was written to remove.
2. **Gap 2 — §5 step 6's error text.** Either make a failed commit name what already
   committed, or amend the plan to say `PendingOutput::commit`'s per-file message is
   sufficient. The ordering half of step 6 is done in all seven multi-output drivers; only the
   promise about the failure message is unkept.
3. **Gap 6 — `guide/outputs.md:124`.** Drop or narrow "a full disk" so the docs claim the
   guarantee §4 item 4 actually delivers ("closed", not "complete"). Three of the four sinks do
   leave nothing behind on a full disk; serial gzip's trailer is the exception.
4. **Gap 5 — T7.** Cheapest close, now that the all-filtered fixture is located: add a presence
   assertion to `tests/integration_preflight_tripwire.rs:660` / `:675`, which already drive every
   pair below `--length 20` and produce zero-record `_val_1` / `_val_2` files that nothing
   currently requires to exist. Failing that, record `integration_output_collision.rs:806`
   (`demux_writes_exactly_the_planned_paths`, which *does* require legitimately-empty outputs to
   be present, in both directions) plus the two zero-record-input tests as discharging the
   overshoot guard, and say so in the plan. Any new wording needs to not read as contradicting
   `integration_partial_output.rs:114-116`, whose whole point is that a valid empty FASTQ must
   **not** appear on a refusal.
5. **Gap 4 — T3's `--dont_gzip` row.** Three lines of test, or one line added to
   "Consciously uncovered".
6. **Gap 1 — the `partial_output_name` note.** One sentence above `src/io.rs:678`. A
   documentation nit; it does not need to gate anything.
7. **T9 (item 20).** Not a gap — `Validate vs Perl TrimGalore` and `Validate uBAM input` both
   exist by name (`.github/workflows/ci.yml:313`, `:957`), but the branch is unpushed, so
   neither has run against `0a6f3a3`. It stays unverified until CI runs.

Nothing in this report is a request to change the design. Gaps 2, 3, 5 and 6 are places where
the plan and the code disagree about what was promised; Gap 4 is a test the plan's own scope
table calls for. Each can be closed without touching `PendingOutput` or the commit points.
