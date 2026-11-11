# PLAN — #428: a mid-stream refusal must leave no output file

**Issue:** [#428](https://github.com/FelixKrueger/TrimGalore/issues/428) (`bug`)
**Branch base:** `dev` @ `285afba`
**Test baseline:** 684 on macOS / **683 on Linux** — see §7 T8
**Version:** v2 (revised after dual plan review; see §12)

---

## 1. Goal

A run that is refused after the first output record has been written must leave **no file at
the output path**. Today it leaves one, and on every output encoding that residue is
structurally valid — a BGZF-EOF-terminated BAM, a well-formed truncated FASTQ, or a
zero-byte FASTQ — so nothing downstream flags it. The only signal is the absence of the
trimming reports, which are written after the loop and therefore never appear on a refused
run.

The fix is a write-then-publish discipline: record output goes to a hidden sibling
temporary and is renamed to its final name only once the driver has closed the writer. A
refused, panicking or erroring run therefore never creates the final name at all.

Secondary goal, and the reason this is worth more than an `unlink()`: a **failed re-run must
not destroy the previous run's good output**. Reviewer A reproduced the harm — a complete
1479-byte output was destroyed by a re-run that then failed, and replaced with a well-formed
103-byte file. `File::create` truncates before the first record is read, so the old result is
gone the moment the new run starts. That is the #409 class of harm, reached with no collision
involved.

---

## 2. Context

### 2.1 The mechanism, measured

`main.rs` opens the reader, then the writer. Every writer constructor creates its file
immediately (`BamWriter::create` at `src/bam.rs:565` does `File::create` + `write_header`;
`FastqWriter::create` at `src/fastq.rs:597` does `File::create`). The read loop then bails,
`?` propagates, and the writer is dropped.

Non-emptiness of the residue does **not** come from the eager header write. It comes from
`Drop`, and all three teardowns were read from the vendored sources by both reviewers:

| Sink | `Drop` behaviour |
|---|---|
| `flate2::write::GzEncoder` | `if inner.is_present() { let _ = self.try_finish(); }` — writes the trailer, **discards the error** |
| `noodles_bgzf::Writer` | `if inner.is_some() { let _ = self.try_finish(); }` — same shape |
| `gzp::ParCompress` | `self.finish().unwrap()` — **panics** on error |

Plain FASTQ needs no trailer at all: a record boundary is a line boundary, so a truncated
file is a well-formed one. The residue is therefore indistinguishable from a short-but-complete
run on **all three** encodings. This is not a uBAM-specific problem.

### 2.2 Why moving the check earlier cannot fix it

Three entry points run before the trimming loop — `sanity_check_any`, adapter
auto-detection (≤1 M reads) and the poly-G scan (≤1 M reads) — and each sees a prefix of the
file, never the whole thing. A fixture whose *first* record is the offender is caught before
any writer exists (`ws_qname_refused_at_sanity_check_entry_point` pins that). A fixture whose
offender is at record 2 with `-a … --no_poly_g` reaches the loop, and so does any file longer
than the 1 M-read scan window. **`-a … --no_poly_g` is not a corner case: it is the arm every
production-scale file takes.**

The aux-tag-value class has no benign regime at all. `sanity_check_any` (`main.rs:39`) opens
the reader **without a tag list**, so a malformed tag value is never looked at during the
peek — it always reaches the trimming loop, and therefore always leaves a residue.

### 2.3 The precedent this closes

`main.rs:401-403`, in the `--rename` + uBAM refusal added by #408 (quote verified exact):

> "The whole run is refused rather than decided per record, because a mid-stream refusal
> would leave a partial BAM behind."

#408 declined a per-record decision *because* of this defect, so fixing it retires a standing
constraint on the design of future refusals. #409 is the precedent for treating partial or
destroyed output on this path as data integrity rather than cosmetics.

### 2.4 Measured scope — which arms leave a residue

All runs: `test_files/ubam_ws_qname_late.bam` (3 records: `readA`, `readB with space`,
`readC`) with `--preserve-tags CB -a AGATCGGAAGAGC --no_poly_g`, binary built at `285afba`.
Every one exits 1.

| Arm | Residue | Shape |
|---|---|---|
| `--output-format ubam` | `…_trimmed.bam`, ≈360 B | valid BAM, BGZF EOF marker present |
| *(default FASTQ out)* | `…_trimmed.fq`, **103 B** | one well-formed record |
| `--dont_gzip` | `…_trimmed.fq`, **103 B** | same |
| `--cores 2` | `…_trimmed.fq`, **0 B** | valid empty FASTQ — reads as "no reads survived" |
| **`--clumpify --cores 2`** | `…_trimmed.fq`, **0 B** | the tenth arm; see the note below |
| `--paired` (interleaved) | `…_val_1.fq` + `…_val_2.fq`, **0 B each** | ditto, ×2 |
| `--paired --output-format ubam` | `…_val.bam`, ≈325 B | valid BAM |
| `--hardtrim5 10` | `….10bp_5prime.fq`, **31 B** | one record |
| **`--hardtrim3 10`** | `….10bp_3prime.fq`, **31 B** | one record |
| `--hardtrim5 10 --output-format ubam` | `….10bp_5prime.bam`, ≈365 B | valid BAM |
| `--clump_only --output-format ubam` | `…_clumped.bam`, ≈305 B | valid BAM |
| `--demux` | `…_trimmed.fq`, **103 B** | bails at the trimming stage, before any barcode file |

**BAM byte counts are approximate on purpose.** The `@PG CL:` line embeds the full command
line, so the same fixture and flags give 359 B or 364 B depending on the length of the `-o`
path. Only the FASTQ counts are stable. Nothing in this plan depends on the numbers — T1/T3
assert absence.

Two corrections to v1's reading of this table, both from Reviewer A:

- **`--clumpify` is reachable and is a live arm.** The refusal at `cli.rs:612-616` sits
  *inside* the `output_format == ubam` block, so `--clumpify --cores 2` with FASTQ output on
  a uBAM input is not refused anywhere. v1 dismissed it on two probes that each failed for a
  different wrong reason (`--clumpify` alone → "requires `--cores >= 2`"; with
  `--output-format ubam` → "not applicable"). Measured above.
- **The de-interleaver *is* the per-record bail site for paired-BAM input.** `bam.rs:261`
  calls `bam_record_to_fastq` on the producer thread; the trimming loop never calls it for
  this input shape at all — it pulls already-converted records off a channel. So the paired
  rows above are *not* "refused for the wrong reason": they are the paired arm's genuine bail
  site. The 0-byte residue comes from batching, not from the de-interleaver as such:
  `BAM_READER_BATCH_SIZE = 4096` (`bam.rs:38`), and the error path returns without flushing
  the partial batch (`bam.rs:264-275`). **Consequence for testing: any paired fixture with
  fewer than 4096 good pairs before the offender lands deterministically in the 0-byte
  regime, which is the one regime where the residue already looks suspicious.** The dangerous
  paired shape — thousands of good pairs written, then a bail, leaving a large plausible
  output — needs >4096 pairs. See §7 T6.

The five per-record bails in `bam_record_to_fastq` sit on the **read** side, so they fire on
every output path. The scope is "every record-output file", not "the uBAM writer".

### 2.5 Output-creation site inventory

**Both reviewers independently verified this table as complete with correct line numbers**,
each enumerating every `File::create` / `FastqWriter::create` / `BamWriter::create` in `src/`
and classifying it against its module's `#[cfg(test)]` boundary. Every omitted site is
inside `mod tests`.

| Module | Sites |
|---|---|
| `src/main.rs` | `:1331` SE serial; `:1586,1587` PE serial R1/R2; `:1597` passthrough; `:1603,1604` unpaired R1/R2; `:1873,1874` PE-from-interleaved-BAM R1/R2; `:1877,1878` its unpaired; `:2148`, `:2307`, `:2440` the three uBAM drivers |
| `src/specialty.rs` | `:54` hardtrim5; `:97` hardtrim3; `:152`, `:204` their uBAM twins; `:270,271` clock; `:383,384` implicon |
| `src/parallel.rs` | `:956` SE parallel; `:212-216` PE parallel R1/R2/passthrough/unpaired×2 (raw `File::create`) |
| `src/demux.rs` | `:197` one per barcode; `:203` nocode |
| `src/clump_only.rs` | `:298` SE; `:440,444` PE; `:767`, `:1005` uBAM |

`parallel.rs`'s output handles are created and written **on the main thread** inside the
`thread::scope` closure (`:212-216`, written at `:238-253`) and never moved into workers, so
holding an RAII guard across them is straightforward.

Reports are **out of scope**: `main.rs:1431/1450/2215/2232/2598/2623`, `demux.rs:283`,
`clump_only.rs:590` all run after the loop returns `Ok`, so a mid-stream bail never produces
a partial report. Their absence is today's only failure signal, and it stays.

---

## 3. Design decision — priced, with a recommendation

### Option A — unlink the output on the error path

Each driver catches the `Err` from its loop, drops the writer, `fs::remove_file`s each output
it created, and re-raises.

- **Diff:** ~15 driver-level sites. No writer API change.
- **Covers:** the ten arms and every other `?` between create and return.
- **Does not cover:** panics, SIGINT/SIGKILL.
- **Does not fix the truncation hole.** `File::create` has already destroyed the previous
  run's output before the first record is read. A failed re-run ends with *nothing* where a
  good file used to be. Unlinking makes the residue honest; it does not give the old file
  back. Worse in one sense: it converts "silently wrong file" into "no file at all" — more
  honest and more destructive at once.
- **Leaves #408's constraint standing**, because a panic still leaves a partial BAM.
- **In A's favour, and worth stating plainly:** A introduces no ordering invariant, no
  drop-before-rename hazard, no new residue class and no tripwire blind spot. Every Critical
  finding in the dual review is a cost of B specifically. A's risk profile is genuinely
  flatter, and that is a correctness argument, not only a diff-size one.

### Option B — write to a hidden temporary, rename on success *(recommended)*

Output goes to a dot-prefixed sibling, `.{file_name}.partial`. One RAII type owns the pair; a
`commit()` renames, and `Drop` removes the temporary if `commit()` never ran.

- **Diff:** one new type in `io.rs`, a `finish()` on `FastqWriter`, ~20 creation sites.
- **Covers:** every early return and every panic that unwinds (one exception below), and
  leaves a *previous good output intact* on a failed re-run — the property A cannot have at
  any diff size.
- **Recoverable where A is destructive.** A 500 GB run that fails at 99 % leaves a
  `.partial` the user can rename by hand and salvage; under A the same run leaves nothing.
  For a trimmer whose performance story is Buckberry-scale files, this is a stronger argument
  than the truncation one.
- **Better for a concurrent reader.** Today `File::create` truncates the previous output in
  place, so anything reading it watches it shrink to zero; under B the old inode survives
  until the atomic rename.
- **Retires the class, not the trigger:** disk-full, a mid-input CRC error, a future
  per-record refusal — none of them leave a valid-looking file, without being enumerated.
  That is what lets #408's constraint at `main.rs:402` be lifted later.
- **Same filesystem by construction**, the temporary being a sibling, so `rename` cannot fail
  `EXDEV` even when `--output_dir` is on another mount.
- **The one panic it does not cover:** `gzp::ParCompress::drop` is `self.finish().unwrap()`.
  A teardown failure panics; if that happens while already unwinding, Rust aborts, `Drop for
  PendingOutput` never runs, and a `.partial` is left. The *final* name is still never
  created, so the guarantee holds — but "covers every panic that unwinds" is one case
  narrower than v1 claimed. Reachable on `cores > 1 && gzip`: `specialty.rs:54/97/270/271/383/384`,
  `demux.rs:197/203`.
- **Cost on the error path:** dropping the inner writer now deflates and writes a full
  trailer — and joins compressor threads under `--cores N` — into a file about to be
  unlinked. A refusal returns measurably slower than today. Correct trade, worth stating so
  nobody files it.

### Ruled out

- **Create the output lazily, on the first record write.** The bail is at record 2, after the
  first write. Same defect as "move the check earlier".
- **Buffer the whole output and write once at the end.** Infeasible for a streaming trimmer.
- **`unlink`-on-`Drop` (an RAII form of A).** Gets A's small diff plus panic coverage, still
  cannot preserve a previous good output, and a SIGKILL leaves the misleading *final* name.
  Strictly dominated by B.
- **`O_TMPFILE` + `linkat(AT_EMPTY_PATH)`.** The only design with genuinely zero SIGKILL
  residue, since the temporary is never linked into the directory. Linux-only — the CI matrix
  includes `macos-latest` — and `linkat(AT_EMPTY_PATH)` needs `CAP_DAC_READ_SEARCH` on many
  kernels.

### Recommendation

**Option B.** Both reviewers independently reached the same verdict: B, conditionally. The
conditions are §5's corrections — above all the commit-point invariant in step 4, without
which B fails *good* runs (§5, C1). B is right and it needs those corrections before it is
safe to implement.

---

## 4. Behaviour

1. A driver asks for an output file by its **final** path, as today. Output naming, the
   collision pre-flight and report naming are unchanged.
2. The bytes go to `.{file_name}.partial` in the same directory, created with `File::create`
   (truncating). The dot prefix keeps a leftover out of `*` globs — see §10 O2.
3. On success the driver commits **at the point where it closes the writer today**. The
   temporary is closed, then renamed onto the final path.
4. What is guaranteed from that instant: **the final name appears only after the writer has
   been closed.** Not "the bytes are complete" — every teardown swallows its error
   (`let _ = try_finish()`, `BufWriter::drop`'s flush), so a disk-full at trailer time still
   publishes a trailerless `.gz` or an EOF-less BAM under the final name. That is *not worse
   than today* (the same bytes already reach the final name today, with exit 0), and making it
   propagable is out of scope — see §10 O5.
5. On any error, early return, or unwinding panic before the commit, the temporary is removed
   and the final path is never created. If a file already existed at the final path, **it is
   untouched**.
6. Multi-output drivers commit each output after the loop returns `Ok`, in a fixed order, with
   **no fallible operation between the first commit and the last** (§5 step 6).
7. Multi-input runs are unchanged: a failure on input 3 of 5 leaves inputs 1–2 committed.
   That matches the existing per-file and "Pair N of M" semantics; cross-file rollback is
   **not** in scope.
8. The temporary is not a user-facing artefact. It appears in error text only if the rename
   itself fails, in which case both paths are named.

### Edge cases

| Case | Behaviour |
|---|---|
| Zero records survive filtering (legitimate empty output) | commit runs; a valid empty output is created, as today |
| Zero-record input | same — writer created and committed, so an empty `.fq`/BAM-with-header appears, as today |
| A user file already sits at the temporary path | overwritten. Acceptable: `<final>` is overwritten unconditionally today |
| A stale temporary from a SIGKILLed run | overwritten by the next run at the same output path; never renamed by anyone else |
| **Abort during unwinding** (`ParCompress::drop`'s `.unwrap()` panicking inside another panic) | no `Drop` runs, so a `.partial` is left. The final name is still never created |
| `--output_dir` on a different mount from CWD or `$TMPDIR` | irrelevant — the temporary is a sibling of the final path |
| Rename fails (permissions, read-only dir) | hard error naming both paths; the temporary is **left in place** so the data is not lost (§5 step 1 makes this work) |
| Commit #2 of N fails after #1 succeeded | hard error naming the failing output and its temporary, which holds its data. The outputs already published stay; **later outputs in the order lose their temporaries**, because propagating the error drops their guards. The error does not enumerate what committed — `?` releases the un-iterated guards before any handler could, so meeting that promise needs a `release_all` on the collection, which is more machinery than the fault deserves (Reviewer B, F6). Perfect N-file atomicity is not achievable without a directory transaction and is not claimed. Still better than today, where the equivalent fault gives a complete `_val_1` and a *silently truncated* `_val_2` at its final name |
| Two runs writing the same output path concurrently | corrupt each other today and under B. B is strictly better for a concurrent *reader* of the previous output (item 4 above) |
| `--fastqc`, `--demux`, `clump_only`'s `output_bytes` | all consume the output by its **final** path inside the same run. This is the invariant in §5 step 4, not an "unaffected" row |

---

## 5. Implementation outline

**Step 1 — `io.rs`: the type.**

```rust
/// A file whose final name appears only once the writer has been closed.
///
/// Writes land on a dot-prefixed sibling; [`commit`] renames onto `final`. If
/// the value is dropped without committing, the temporary is removed and
/// `final` is never created.
pub struct PendingOutput {
    tmp: PathBuf,
    final_path: PathBuf,
    released: bool,
}

impl PendingOutput {
    pub fn create(final_path: &Path) -> Result<(Self, File)>;
    pub fn commit(mut self) -> Result<()>;
    pub fn tmp_path(&self) -> &Path;
}

impl Drop for PendingOutput { /* remove tmp when !released; ignore errors */ }
```

> **`released` is set *before* the `rename` is attempted, not after.** `PendingOutput`
> implements `Drop`, so `commit` cannot destructure out of `self`; the flag is the only
> mechanism. If it were set only on success, a failed rename would return `Err`, `self` would
> drop with `!released`, and `Drop` would **delete the temporary** — the exact opposite of
> §4's promise that a failed rename leaves the data in place. The field name states the
> meaning: "no longer owned by `Drop`".

`partial_output_name(final_path)` returns `parent.join(format!(".{}.partial", file_name))` —
dot-prefixed, and appended to the *whole* file name so it never collides with an
extension-stripping namer. Put the "no suffix may be appended to a committed name" note above
that function.

**Step 2 — `FastqWriter`: hold the pending output, gain `finish()`.**

Change `writer: Box<dyn Write + Send>` to `Option<Box<dyn Write + Send>>` and add
`pending: Option<PendingOutput>`. Declare `writer` **before** `pending` so declaration-order
field drop closes the file before the temporary is removed.

```rust
/// Flush, close the compressor, and publish the output under its final name.
/// Consumes `self` so a caller cannot forget. Mirrors `BamWriter::finish`.
pub fn finish(mut self) -> Result<()>
```

> **Two distinct hazards live here, and v1 stated only the first.**
>
> **(a) Rename before drop publishes a truncated file.** The gzip trailer and `gzp`'s final
> blocks are written by the inner value's `Drop`, so `finish()` must `Option::take()` the box
> and drop it *before* renaming. `take()` needs only `&mut self`, which is what makes it work
> inside a type that implements `Drop`.
>
> **(b) Drop-without-propagating publishes a truncated file *and reports success*.** A trait
> object exposes only `Write`, and `Write::flush()` on `BufWriter<GzEncoder<File>>` does not
> write the gzip trailer — only `try_finish()` does, and it is unreachable through
> `Box<dyn Write>`. So on a trailer-write failure `finish()` returns `Ok(())` and commits a
> trailerless `.gz`. This is **not a regression** — today the same bytes reach the same final
> name with exit 0 — and it is why §4 item 4 states the guarantee as "closed", not "complete".
> The fix is the sink enum in §10 O5, deliberately **out of scope** for this PR.
>
> Hazard (a) is caught by T2. Hazard (b) is structurally invisible to any success-path test.

`BamWriter` is the model for (a) and already correct: it holds the concrete type, so its
`finish()` calls `try_finish()` and propagates (`bam.rs:644-649`).

**Step 3 — `BamWriter`: same shape.** Add `pending: Option<PendingOutput>`; in `finish`,
after `try_finish()`, drop `self.inner`, then commit. Same drop-before-rename rule.

**Step 4 — the raw `File::create` record-output sites, and the commit-point invariant.**

> **The invariant, which is what makes Option B safe:** FastQC, `--demux`,
> `write_clump_only_report` and `clump_only`'s `output_bytes` all consume the output **by its
> final path, within the same run**. So the commit must happen exactly where the writer is
> closed today — never later, and never at the end of the enclosing function.

The consumers, all verified at `285afba`:

| Site | What it does with the final path |
|---|---|
| `main.rs:1459`, `:2241` | `fastqc::run(&output_path, …)` |
| `main.rs:1474` | `demux::demultiplex(&output_path, …)` — reads the trimmed output, fans out per-barcode files |
| `clump_only.rs:374`, `:847` | `fastqc::run(&output_path, …)` |
| **`clump_only.rs:827`, `:1071`** | **`stats.output_bytes = std::fs::metadata(&output_path)?.len();`** |

The last is the sharp one: a `?`-propagating stat of the final path on the **success** path. A
late commit fails a **good** run with `ENOENT`. The writer-level design in steps 2–3 is safe
only because the commit lives inside `finish()`, which is called where `drop` is today.

So for the raw-`File` sites, commit at today's teardown point, named explicitly:

- `clump_only.rs:298` (SE FASTQ) → commit at `out.flush()?; drop(out);` (`:347-348`), which
  is *before* `fastqc::run` at `:374`.
- `clump_only.rs:440,444` (PE FASTQ) → commit at `out_r1.flush()?; out_r2.flush()?;`
  (`:517-518`).
- `clump_only.rs:767`, `:1005` (uBAM) → inherit step 3 via `writer.finish()?` (`:822`), which
  is *before* the `metadata` stat at `:827`.
- `parallel.rs:212-216`, `:956` → bare `File` handles, no user-space buffer and no trailer, so
  there is nothing to order; commit after the concatenation loop returns `Ok`.

**Do not remove the existing `flush()` calls** on the grounds that `PendingOutput` now owns
the file. They are what marks the commit point.

**Step 5 — every driver calls `finish()`.** ~15 sites. `main.rs:1331-1334`'s
`writer.flush()?; drop(writer);` becomes `writer.finish()?`. Because `finish` consumes
`self`, the compiler enumerates every site that still relies on `Drop`.

**Step 6 — multi-output commit ordering.** Commit N outputs in a fixed order (R1, R2,
passthrough, unpaired R1, unpaired R2) after the loop, with **no fallible operation between
the first commit and the last** — no report write, no stats computation — so the non-atomic
window stays minimal. A mid-sequence failure propagates `PendingOutput::commit`'s own error,
which names the failing output and the temporary holding its data; it does not enumerate what
already committed (see §4's edge-case row for why that promise was withdrawn). Do not
un-commit: deleting a successfully published file on an I/O fault would be a second data loss.

**Step 7 — `main.rs:395-405`.** Leave the `--rename` refusal in place; this PR does not relax
it. Update only the clause stating the reason, so it stops asserting a defect that no longer
exists. One sentence, current state only.

**Step 8 — the residue-asserting test.** §7 T1.

**Step 9 — the four places that document the contract this fix reverses.** v1 called a single
docs sentence "the only user-visible surface", which was wrong:

| Location | Text that becomes false |
|---|---|
| **`src/cli.rs:295`** | `--passthrough` help: "On a mid-stream reader error … partial output files may remain on disk — re-run after fixing the input." This is **`--help` output** |
| `docs/src/content/docs/modes/passthrough.md:71` | "Mid-stream errors leave partial output on disk … files already-written to disk are not rolled back" |
| `src/parallel.rs:90` | "partial outputs may still exist on disk on mid-stream error — v1 contract" |
| `src/parallel.rs:263` | "Partial output files may remain on disk (v1 contract; Cli --help documents this)" |

Leave the archived changelog entry at `docs/.../reference/changelog.md:829` alone — it
describes a shipped version, and the synced changelog page is not to be edited.

**Step 10 — the tripwire header.** Add `.partial` to `tests/integration_preflight_tripwire.rs`'s
"Known-unplanned writes" section with the injectivity argument (§11), plus one clause: a red
`.partial` means a **missed commit**, not a missing candidate. Without it, the tripwire's own
remediation message points at the wrong fix.

**Step 11 — docs + `CHANGELOG.md`.** One sentence in `docs/.../guide/outputs.md`: a refused
run leaves no output file. One changelog sentence, existing format, no migration text.

---

## 6. Efficiency

One extra `rename(2)` per output file per run — a directory metadata operation, unmeasurable
against a trimming run. No extra copy: the temporary is a sibling, so the rename never falls
back to copy-and-delete. No extra memory, no change to the compression path.

The **error** path gets slower: dropping the inner writer deflates and writes a full trailer,
and joins compressor threads under `--cores N`, into a file about to be unlinked. Pointless
work, but skipping it would require a non-`Drop` teardown, i.e. §10 O5's enum. Worth one
sentence in the PR so nobody files the slower refusal as a regression.

**Peak disk usage on an in-place re-run doubles.** `File::create` used to truncate the previous
output before the first record was read, so a re-run's peak footprint was one copy. The previous
output now survives until the rename, so a re-run over an existing result needs room for both.
That is the same mechanism as the recovery property §3 sells B on, but at Buckberry scale a
re-run that used to fit can now hit `ENOSPC` — and per §4 item 4 that is the failure mode which
publishes a trailerless file with exit 0 on the serial-gzip path. Both reviewers raised it
independently (A-M2, B-F7). PR clause, not a docs change.

---

## 7. Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| **T1** | The residue is gone on the arm that pinned it | invert `tests/integration_ubam.rs::ws_qname_refused_in_trimming_loop_ubam_out_leaves_partial` → `…_leaves_no_output`; keep the `BAM record 2` stderr assertion, replace the `BGZF_EOF` block with an assertion that neither `…_trimmed.bam` nor any temporary exists | passes only with the fix |
| **T2** | **The trailer is not lost to the rename** — three arms, not one | (a) serial gzip output; (b) uBAM read back with `assert_ubam_eq`; (c) **`--cores 2` + gzip on a specialty or `--demux` path**, which is the only way to exercise `gzp::ParCompress`'s teardown. Assert **decoder-level integrity** (`gzip -t`, or a `MultiGzDecoder` read-to-end that propagates the trailer error) **in addition to** the record count | full record count *and* a valid trailer. A record count alone passes on a trailerless gzip — the exact failure T2 exists to catch |
| **T3** | Every measured arm in §2.4 leaves nothing | one case per arm: FASTQ serial, `--cores 2`, **`--clumpify --cores 2`**, `--hardtrim5` FASTQ + uBAM, **`--hardtrim3`**, `--clump_only --output-format ubam`, **`--demux`** | no output file **and** no leftover temporary — assert both halves |
| **T4** | **A previous good output survives a failed re-run** (B only; the discriminating property) | run to success, record the output's bytes **and its report pair**; re-run at the **same output path** (same input filename — copy the offending fixture onto the successful run's stem) with the offending fixture; compare | byte-identical output and reports. Written with two different input names it compares unrelated files and passes under Option A too |
| **T5** | `--cores 2`'s zero-byte residue is gone | the `--cores 2` case in T3 | file absent, not merely empty |
| **T6** | The paired arm, in the regime where the harm is real | `test_files/ubam_ws_qname_paired.bam` **already exists** and already drives `ws_qname_refused_via_interleaved_deinterleaver` — assert absence there. **Then add the case that matters:** >4096 pairs before the offender, so the residue is large and plausible rather than 0 B | neither `_val_1/2.fq` nor `_val.bam` exists, in both the small and the >4096-pair regime |
| **T7** | Legitimate empty output still appears | confirm at least one existing all-reads-filtered test asserts the file's **presence** | unchanged — an empty output file is created. This is the guard against the fix overshooting into "leaves no file when the output is legitimately empty" |
| **T8** | No regression in the suite | `cargo test --release`. **The baseline is platform-dependent** — `src/cli.rs:1950` is `#[cfg(target_os = "linux")]` and `format.rs` has two `#[cfg(target_os = "macos")]` tests, so macOS reports 684 and Linux 683. Do not gate on a single hard-coded number | zero failures on both platforms |
| **T9** | CI's byte-identity invariant holds | `Validate vs Perl TrimGalore` + `Validate uBAM input` by **name**, not job count | green |
| **T10** | `--fastqc` still produces a correctly *named* artefact | the `--fastqc`-active arms of the tripwire cover this via T1; assert it rather than resting on a table row | `x_trimmed_fastqc.html`, not `x_trimmed.fq_fastqc.html`. This is the concrete user-visible symptom of a commit-ordering slip |

**Consciously uncovered**, so the coverage audit does not have to rediscover them: `--clock`,
`--implicon`, `--passthrough` and `--retain_unpaired`. Each is the same code shape as a
covered arm, and `finish()`-consumes-`self` plus the compiler is the guard. `--demux` is
*not* in this list — it is the only stage that reads a committed output and then writes more,
so it gets a real case.

`cargo fmt --all -- --check` and `cargo clippy --all-targets --release -- -D warnings` are
gates, not validations.

---

## 8. Assumptions

1. **Fixed:** the five bails are read-side and reach every output path — including the
   paired arm, where the bail site is the de-interleaver's producer thread (§2.4).
2. **Fixed:** reports are never partial **from a mid-stream bail**, because they are written
   after the loop. A disk-full *during* `report::write_report_header` does leave a truncated
   report at its final name; reports stay out of scope (a footerless text report is
   self-evidently truncated, and it is not the data).
3. **Fixed:** `rename` is atomic within a directory on both CI platforms — `ci.yml:28` is
   `os: [ubuntu-latest, macos-latest]`, no Windows — so unlink-while-open and
   rename-while-open are both safe.
4. **Fixed:** the temporary is a sibling of the final path, never in `$TMPDIR`. This is what
   makes the `EXDEV` claim sound.
5. **Fixed, as an invariant rather than a non-existence claim:** four consumers read the
   output by its final path *inside the same run*, so **the commit happens exactly where the
   writer is closed today, never later**. v1 asserted no such consumers existed; that was
   false, and `clump_only.rs:827`/`:1071`'s `?`-propagating `metadata` call means getting it
   wrong fails a good run. See §5 step 4.
6. **Fixed:** `FastqWriter` is only ever constructed from a path — there is no `new`/stdout
   constructor. Option B depends on it: a writer over a non-file sink could not hold a
   `PendingOutput`. Stated because it is a constraint on future work.
7. **Fixed:** field declaration order (`writer` before `pending`) closes the file before the
   temporary is removed. On POSIX the reverse is harmless — writes land on the unlinked inode
   — so this is a free invariant, not a correctness risk.
8. **Configurable → deliberately not:** the temporary name is a compile-time derivation. No
   flag.
9. **Assumed:** overwriting a pre-existing temporary is acceptable, since `<final>` is
   overwritten unconditionally today.
10. **Assumed:** cross-input rollback is out of scope (§4.7).
11. **Fixed:** no new runtime dependency. `tempfile` is a **dev**-dependency
    (`Cargo.toml:50`/`:52`); the temporary name is derived, so nothing needs promoting.

---

## 9. Notes from planning and review — not filed

Per the filing bar (data-integrity or user-visible defects only; one batch if several
accumulate):

- **`--cores 2` leaves a zero-byte output rather than a short one** — same defect, different
  face. T5 covers it.
- **§2.4's BAM byte counts are not reproducible** — the `@PG CL:` line embeds the command
  line, so the count varies with the `-o` path length. Presentation fixed in v2.
- **`ws_qname_refused_in_trimming_loop_fastq_out` asserts only the refusal**, never the
  residue, so the FASTQ arm was unpinned in both directions. T3 closes it. A test gap.
- **`--gzip` (hidden, `cli.rs:425`) is consulted only at `cli.rs:1102`**, its deprecation
  notice. A documented no-op behaving as intended.
- **`the_site_inventory_is_complete`** in the tripwire pins the number of
  `preflight_output_collisions` call sites. Option B adds none, so it is unaffected.
- **Behaviour changes at the output path, all inherent to publish-by-rename, none previously
  written down** (A-L4, B-F8, both verified empirically): a **symlink** at the output path is
  replaced by a regular file rather than written through, so the link's former target keeps the
  old bytes; a **write-protected (0444) existing output is now silently replaced** where
  `O_TRUNC` used to return `EACCES` and refuse the run; and any non-default **mode** on a
  pre-existing output is reset. The hard-link case changes in the better direction — a
  previously hard-linked copy keeps its content instead of silently seeing the new bytes. PR
  clause.
- **A 9-character band of legal output names now fails with `ENAMETOOLONG`** (B-F2, measured by
  direct syscall). The temporary is nine bytes longer than the final name, so where the final
  name is 247–255 bytes it exceeds `NAME_MAX` and the run dies at output-creation time. Only
  reachable where no report is written — with reports on, `{input}_trimming_report.txt` was
  already too long — so `--hardtrim5/3`, `--clock`, `--implicon` and `--no_report_file`. It
  fails loudly, before any work, with no data loss; no derived-name scheme avoids it, and
  shortening the base name is what the pre-flight's injectivity argument depends on *not*
  happening. Accepted, with a PR clause.
- **Removing `pub fn flush` is a breaking change to the published library surface** (B-F9).
  `trim-galore` publishes a lib; practical impact is near zero for a
  `command-line-utilities` crate. PR clause, deliberately **not** a CHANGELOG entry — that file
  is for user-facing behaviour, one sentence each.
- **Duplicate barcodes in a `--demux` samplesheet now lose a planned file outright** (A-L6).
  `read_barcode_file` validates alphabet and length but not uniqueness, and `writers` is keyed
  by barcode, so a second entry with the same barcode drops the first writer — previously an
  empty file at the first sample's final name, now nothing. Pre-existing missing validation
  whose symptom this change alters; arguably the more honest symptom. Below the bar alone.
- **Two concurrent runs at the same output path share one temporary** (B-F12), so the second's
  rename fails `ENOENT` while the error says the data is intact at a path that no longer
  exists. Both runs already corrupt each other today via the shared final path, so not a
  regression — but it is the one case where the new message can lie.
- **The one candidate that *does* clear the bar, and is not part of this PR:** every teardown
  swallows its error, so a disk-full or EIO at trailer time makes `trim_galore` **exit 0 while
  leaving a silently truncated `.gz`** (`flate2`: `let _ = try_finish()`) or an EOF-less BAM,
  and on `--cores N` gzip paths `gzp::ParCompress::drop` **panics** instead. That is a
  user-visible data-integrity symptom, it exists today independently of #428, and it is what
  §10 O5's sink enum would fix. Both plan reviewers *and* both code reviewers found it
  independently. File it as **one** issue covering all three sinks once #428's PR lands — not
  three, and not inside #428. Two further bullets belong on that same issue rather than on new
  ones: the **SE parallel path** only `eprintln!`s a worker error and then publishes a short
  output with **exit 0** (A-N1, B-F11) — the paired path fixed exactly this by putting
  `Result<…>` in the channel, and the single-end twin never got the same treatment; and a
  **worker panic in the paired parallel path** commits before `thread::scope` re-raises, giving
  a truncated output at the final name with exit 101 (A-N2). Both are pre-existing — the same
  bytes reached the same final name before — and both are the same family: an error discarded
  at teardown or in a worker, producing a truncated file under the final name.

---

## 10. Questions

- **O1 — answered by both reviewers: Option B, conditional on §5's corrections.** The
  conditions are the commit-point invariant (step 4), the `released`-before-rename rule
  (step 1), and the scope corrections in §2.4. Felix's sign-off is still wanted on B over A,
  since §3 now states A's flatter risk profile as a real argument.
- **O2 — answered: a dot-prefixed sibling, `.{file_name}.partial`.** Decided on glob
  visibility, not legibility. Measured: `*.fq*` matches `x_trimmed.fq.gz.partial`, and a
  dot-prefixed name matches neither `*.fq*` nor `*.fq.gz`. TrimGalore's dominant consumer is
  nf-core, where outputs are collected by glob and `publishDir` patterns, so a leftover under
  a `*`-visible name can silently enter a pipeline's output set. The honest trade-off: a
  dot-prefixed leftover is invisible to the *user* too, so a SIGKILLed run leaves silent junk
  rather than a breadcrumb. The glob risk is sharper — hidden junk only wastes disk.
- **O3 — Open, assumption taken: no.** The `--rename`/uBAM refusal at `main.rs:395-405` is
  not relaxed here; its substantive reason (representability depends on the individual header)
  stands on its own. Step 7 edits only the clause.
- **O4 — dissolved, no decision needed.** `test_files/ubam_ws_qname_paired.bam` already
  exists and is already the fixture v1 proposed, already driven by
  `ws_qname_refused_via_interleaved_deinterleaver`. (v1 also miscounted the fixture set:
  `test_files/` holds 13 `.bam` fixtures, not eight.) What T6 still needs is the >4096-pair
  case.
- **O5 — Open, assumption taken: out of scope.** Replacing `Box<dyn Write + Send>` with an
  enum over `Plain`/`Gz`/`ParGz` is the only way to make teardown errors propagable, and it
  is where the two reviewers disagreed: one wanted it inside this PR as the thing that makes
  B's guarantee real; the other wanted it deferred, with §4's claim weakened instead. Taken:
  **defer.** It touches `write_record` on the hot path so it needs `scripts/benchmark.sh`
  numbers rather than an assertion, its failure mode needs a size-limited filesystem to test,
  and it is orthogonal to the residue question. §4 item 4 states the weaker, true guarantee
  instead. Tracked by the §9 issue.

---

## 11. Self-review

**Efficiency.** One `rename` per output file; nothing in the hot loop changes. The error path
is slower by one trailer write and a thread join — stated in §6 rather than discovered later.

**Logic.** Traced create → write → error and create → write → commit for all five writer
kinds. Two things v1 got wrong and v2 fixes: the commit point is an *invariant* set by four
in-run consumers, not a free choice (§5 step 4), and the `released` flag must precede the
rename or `Drop` destroys the data §4 promises to keep (§5 step 1).

**Edge cases.** Added the abort-during-unwinding row, the concurrent-writer/reader row, and
the N-file partial-commit row. Kept the zero-record and all-filtered rows: "leaves no file"
must not become "leaves no file when the output is legitimately empty" (T7).

**Integration.** No change to output naming, the collision pre-flight, report naming, or
candidate order — which matters because #423/#424 land next and #423's constraint is that
candidate *order* is load-bearing. The temporary names stay out of the pre-flight because
`final → temporary` is **injective**, and injective under `collision_key`'s case-folding too:
if no two finals collide, no two temporaries can, so adding them buys zero detection power.
v1 gave test churn as the reason, which invites a future maintainer to "fix" the exclusion by
adding candidates and re-pinning. (v1 also named the wrong test: the pinned one is
`tests/integration_preflight_tripwire.rs`, which consumes `DUMP_PLANNED_ENV` and asserts
`created ⊆ planned`, not `integration_output_collision.rs`.) Reviewer A verified Option B does
not break the tripwire: every case asserts exit 0 and snapshots after the process exits, so a
committed temporary never appears in `created`.

**Remaining risks.** (a) Hazard (b) in step 2 — a swallowed trailer error still publishes a
truncated file with exit 0. Not a regression, explicitly deferred to O5, and the reason §4
says "closed" not "complete". (b) ~20 mechanical sites, where a missed one is a silent hole;
mitigated by `finish()` consuming `self`. (c) `--demux`'s per-barcode fan-out is the widest
commit set and now gets a real T3 case rather than an admission. (d) The >4096-pair paired
case (T6) is the only test that exercises a large plausible residue; if it is dropped for
runtime cost, the plan covers file *creation* everywhere and partial *content* nowhere.

---

## 12. Implementation notes (2026-08-14)

Implemented on `fix/428-no-partial-output-on-refusal` off `dev` @ `285afba`.
**707 tests pass** (684 → +23), `cargo fmt --all -- --check` clean, `cargo clippy
--all-targets --release -- -D warnings` clean.

### Post-review round (dual code review + coverage audit)

`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md` and `COVERAGE.md`. Both reviewers verified the
commit-point invariant by control flow rather than line order, confirmed no output byte changes
on any success path, and judged all four deviations below correct. Neither found a Critical.
The coverage audit returned INCOMPLETE — 4 of 42 — with no item MISSING.

**Deviation 2's safety argument was wrong, and both reviewers proved it by mutation.** §12 first
claimed the compiler enumerates the ~20 mechanical sites. It does not: replacing a `finish()?`
with `drop(w)` compiles, and where no test asserts the file exists, *nothing* protects it.
Reviewer A and Reviewer B independently mutated `--retain_unpaired`'s two writers to `drop(w)`
and the whole suite stayed green. What guards those sites is **tests**, not the compiler.
Deleting `flush()` still earns its place — it forces each site to call *something* — but the
argument stops there. Under this change a missed commit is also *worse* than before: the file is
absent rather than present-and-empty.

Fixes applied in this round:

- **The docs claimed a full disk leaves nothing behind.** It is the one case in that list this
  change does not cover, and it fails with exit 0. Reviewer A traced it per sink: plain FASTQ
  propagates via `flush()?`, uBAM via `try_finish()`, `--cores N` gzip is covered *by panicking*
  (`ParCompress::drop`'s `.unwrap()` unwinds while the taken `pending` is a live local, so its
  `Drop` removes the temporary) — only serial gzip loses the trailer silently. Clause deleted.
- **A fifth statement of the reversed contract survived**, at `docs/.../guide/outputs.md:145`,
  twenty lines below the new "Refused runs" section in the same file, so the page contradicted
  itself. The first sweep grepped `partial output`; this said "partial **BAM**". Rewritten as
  `main.rs:399-403` was. The lesson — sweep the varying noun, not the phrase you wrote, and
  grep the file you just edited first — is memorised.
- **The "safe to delete" advice could have cost data.** A failed rename deliberately leaves the
  temporary and the error says so, which the docs contradicted (B-F1). Now both are stated.
- **Both ordering rules had zero coverage.** Two new tests, each proved load-bearing by
  re-running the reviewers' own mutations: `io::tests::a_failed_commit_leaves_the_temporary_on_disk`
  (a non-empty directory at the final path makes `rename` fail; goes red when `released` is set
  after the rename) and `retain_unpaired_publishes_both_rescued_outputs` (goes red on `drop(w)`).
- **Two absence assertions could have passed vacuously.** `refuse()` now pins the record
  *ordinal*, not just the refusal reason — if the whitespace check ever moved to
  `sanity_check_any`, every absence assertion in the file would have passed with no writer ever
  created. `refused_clump_only_publishes_nothing`, the one case that bypasses `refuse()`, now
  asserts why the run failed. `refuse()`'s unused `-> String` is gone.
- **`partial_output_name` was lossy** (A-L1): `to_string_lossy` maps distinct invalid-UTF-8
  names onto one temporary, so the injectivity claim the tripwire header now rests on held only
  for UTF-8. Rebuilt on `OsString` — byte-exact, unconditional. Two unit tests pin the naming.
- **`create`'s error named the final path** while the failing syscall targeted the nine-bytes-longer
  temporary (A-L3, and the same nine bytes as B-F2's `ENAMETOOLONG` band). Both are named now.
- **T7, `--dont_gzip` and the coverage audit's remaining gaps closed** — see the ledger.
- **Three comments exceeded the house rule** and both reviewers flagged the same three
  (`fastq.rs::finish`'s doc, `demux.rs`'s four-liner, `io.rs`'s glob rationale). Trimmed.
  `fastq.rs::finish` also *overstated its mechanism*: Reviewer B showed that committing before
  the drop keeps the suite green, because on POSIX `rename` does not invalidate an open
  descriptor — the trailer still lands in the same inode. The current order is right precisely
  because it survives a SIGKILL in that window, and the comment now says only that.
- **`PendingOutput::tmp_path()` was dead code** (A-L2) — kept, because the new `commit`-failure
  test uses it.

Withdrawn rather than implemented: §5 step 6's "error naming what committed and what did not".
Reviewer B showed the natural implementation cannot work — `?` drops the un-iterated guards
before any handler runs, so a mid-sequence failure additionally *removes* the later outputs'
temporaries. Meeting the promise needs a `release_all` on the collection; the plan now describes
what ships instead.

### Deviations from the plan

1. **`Option<Box<dyn Write + Send>>` was not needed.** Step 2 called for wrapping the inner
   writer in an `Option` so `finish()` could take it out before renaming. Only `pending`
   needs to be an `Option`: `finish(mut self)` takes `pending`, then `drop(self)` — which
   closes the compressor and writes the trailer — then commits. Same ordering guarantee,
   no `expect`/panic path on `write_record`, smaller diff. `BamWriter::finish` uses the
   identical shape.
2. **`FastqWriter::flush()` was removed, and this is what made the fix safe to land.** The
   plan claimed `finish()` consuming `self` would make the compiler enumerate every site
   still relying on `Drop`. **That was false:** `writer.flush()?; drop(writer);` compiles
   fine and would silently publish nothing. Deleting the public `flush()` turned the
   enumeration into a build error and listed all 18 sites (8 in the lib, 10 in the
   binary). Nothing needed a mid-stream flush, so nothing was lost.
3. **T6's >4096-pair paired case became a 5000-record single-end case.** `record()` in
   `examples/mk_bam_fixture.rs` hardcodes `FLAG_UNMAPPED`, so mate-adjacent paired records
   need a generator change; and paired batching (`BAM_READER_BATCH_SIZE`) only makes a
   paired residue *smaller* than the single-end equivalent. New fixture
   `test_files/ubam_ws_qname_bulk.bam` (14,860 B, byte-reproducible from
   `cargo run --release --example mk_bam_fixture -- ws_qname_bulk`) puts 5000 clean records
   ahead of the offender, so `refused_after_thousands_of_good_records_publishes_nothing`
   covers the large-plausible-residue shape. The paired arms still get absence cases via
   the existing `ubam_ws_qname_paired.bam`.
4. **T2's decoder-level integrity uses `trim_galore::fastq::FastqReader`**, not `gzip -t` —
   the same `MultiGzDecoder` the tool reads with, and no external test dependency.

### Verification beyond the suite

**Mutation-tested, because a clean run proves nothing.** Planting the pre-#428 behaviour —
`Drop for PendingOutput` renaming instead of removing — turned **15 tests red**: all 13
refusal cases, T4, and both inverted `integration_ubam` cases. The 4 that stayed green are
exactly the success-path publish cases, which do not touch the abandon path. `src/io.rs`
was restored from a byte-compared backup afterwards (md5 verified, zero `MUTANT`
occurrences remaining).

**Clippy coverage of the new test target was proved, not assumed** — a planted `len_zero`
violation in `tests/integration_partial_output.rs` failed the run; the clean run that
followed is therefore meaningful rather than a cache hit.

**Commit ordering verified against every in-run consumer.** `clump_only.rs`: commit `:349`
< `fastqc::run` `:375`; `:521-522` < `:555-556`; `writer.finish()` `:823` < `metadata`
`:828` < `fastqc` `:848`; `:1069` < `:1072` < `:1096`. `main.rs`: `:1333` < `demux` /
`fastqc` `:1458`; `:2153` < `:2239`; `:2317` < `:2396`; `:2446` < `:2519`.

**The old-claim sweep was run repo-wide**, not just on the sections the review pointed at.
Three surviving matches for "partial output" are about different claims — a *hang* leaving
output (`integration_non_restartable_input.rs:177/206`) and the multi-pair pre-flight
(`integration_paired_format_guard.rs:13`) — both still true.

### Not done

The sink enum (§10 O5) stays out of scope, so a teardown error at trailer time still
publishes a trailerless file with exit 0. §4 item 4 states the weaker guarantee that is
actually delivered. The §9 issue covering all three sinks is still to be filed.

## 13. Revision history

**v2 (2026-08-14)** — revised after dual independent plan review
(`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B_gen2.md`). Both reviewers verified the §2.5 site
inventory as complete with correct line numbers, and both answered O1 as "B, conditionally".
Substantive changes: Assumption 5 replaced by the commit-point invariant after A found four
in-run consumers of the final path, one of them a `?`-propagating `metadata` on the success
path (§5 step 4); `--clumpify` and `--hardtrim3` added to §2.4 after A showed `--clumpify` is
reachable with FASTQ output on uBAM input; the paired-arm diagnosis corrected — the
de-interleaver *is* the bail site, and the 0-byte residue comes from `BAM_READER_BATCH_SIZE`
batching, which makes the >4096-pair case the only one testing a large plausible residue;
O4 dissolved (the fixture exists); `released`-before-rename rule added (§5 step 1); the
"complete" claim in §4 weakened to what is actually delivered, with the sink enum deferred to
O5; T2 given a third arm and decoder-level integrity; step 9 extended from one docs sentence
to the four places that document the reversed contract, one of them `--help`; §11's pre-flight
reasoning corrected to injectivity and the right test name; O2 answered on glob visibility
with a dot-prefixed temporary; T8's baseline noted as platform-dependent (684 macOS / 683
Linux).

**v1 (2026-08-13)** — initial plan.
