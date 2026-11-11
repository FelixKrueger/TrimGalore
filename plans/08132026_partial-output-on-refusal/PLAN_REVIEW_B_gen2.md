# PLAN REVIEW — Reviewer B (generation 2)

**Plan:** `plans/08132026_partial-output-on-refusal/PLAN.md` (#428)
**Tree:** `dev` @ `285afba`, read-only. Binary `target/release/trim_galore` verified to report
`285afba — macos/aarch64 — built 2026-08-13T20:16:22Z`.
**Method:** every load-bearing claim re-derived from source at `285afba` or measured with that
binary. Probes under `$TMPDIR/revb2/`. No build, no test run, no edit.

**Verdict on O1: Option B — but not as specified.** Two of B's three headline properties do not
hold as written. The drop-to-finalise mechanism in §5 step 2/3 silently discards the error from
the very write that completes the file, and then renames the truncated result into place (C1).
The compile-time enforcement §11 relies on does not exist (C2). Both are cheap to fix inside B,
and B's core advantage — a previous good output surviving a failed re-run — is real and I
measured it. Fix C1 and C2 and B is the right call.

---

## 0. What I verified, and what came back different

| Plan claim | Verdict |
|---|---|
| `FastqWriter::create` does `File::create` at `fastq.rs:597` | ✅ exact |
| `BamWriter::create` does `File::create` + `write_header` at `bam.rs:565` | ✅ exact |
| `main.rs:402` quote | ✅ verbatim (lines 401–402) |
| `tempfile` is a dev-dependency at `Cargo.toml:52` | ✅ `[dev-dependencies]` opens at line 50 |
| No Windows in CI | ✅ `os: [ubuntu-latest, macos-latest]` (ci.yml:28); all 7 other jobs `ubuntu-latest` |
| §2.5 site inventory | ✅ complete and line-accurate for production sites (see §1.4) |
| Reports are written after the loop, never partial | ✅ all 8 report sites post-loop |
| `flate2::GzEncoder::drop` finalises | ✅ `flate2-1.1.9/src/gz/write.rs:161` → `try_finish()` |
| `gzp` parallel blocks come from `Drop` | ✅ `gzp-2.0.2/src/par/compress.rs:391` → `self.finish().unwrap()` |
| `rename` cannot hit `EXDEV` (sibling) | ✅ for a plain path — ❌ for a symlinked output path (I4) |
| §2.4's nine arms leave a residue | ✅ all nine reproduce (BAM byte counts differ, I7) |
| A failed re-run destroys the previous good output | ✅ **measured**: 1639 B → 103 B |
| `--clumpify` is unreachable from a BAM input | ❌ **false** — reachable, and it is a tenth arm (C3) |
| The `--paired` rows were refused by mate-adjacency | ❌ actual error is invalid flag bits (I3) |
| T6 needs a new ninth BAM fixture (O4) | ❌ the fixture already exists (I3) |
| §2.1's three residue encodings are established | ⚠️ gzip was asserted, not measured — I measured it, and it is the worst case (I1) |
| Step 9 is "the only user-visible surface" | ❌ `--help` and a docs page already promise the opposite (I2) |
| §6: no efficiency cost beyond one `rename` | ❌ peak disk doubles on a re-run (I8) |
| §11: the compiler enumerates the ~20 sites | ❌ it does not (C2) |

---

## 1. Logic review

### 1.1 The mechanism as described is right, and its strongest evidence is missing

§2.1's causal chain is correct in every particular: writers create eagerly, the bail propagates,
and non-emptiness comes from `Drop` rather than the header write. I confirmed all three
finalisers.

But §2's *evidence* is nine arms of a 3-record BAM fixture, and the largest residue in the table
is 359 bytes. That understates the defect by five orders of magnitude, because a BAM input can
never produce gzipped FASTQ output — `io::is_gzipped` is extension-based (`io.rs:33`), `x.bam`
has no `.gz`, so the output mode is always plain. The gzip arm is therefore unreachable from the
five bails, which is why it is absent from the table, and §2.1's `flate2` bullet is the one claim
in that section with no measurement behind it.

It is reachable from a *FASTQ* input, with a fixture already in the repo:

```
$ trim_galore -a AGATCGGAAGAGC --no_poly_g -o D test_files/truncated.fq.gz
Error: unexpected end of file            (exit 1)
$ ls -l D/truncated_trimmed.fq.gz        →  6,124,349 bytes
$ gunzip -t D/truncated_trimmed.fq.gz    →  clean
$ gunzip -c D/truncated_trimmed.fq.gz | wc -l   →  188,330 records
$ tail -c8 …                             →  intact CRC32/ISIZE trailer
```

`--cores 2` gives the same shape (6,000,099 B, 184,286 records). A downstream consumer receives a
valid, fully-terminated, 188,330-read trimmed FASTQ with nothing at all to distinguish it from a
complete run. That is the residue class at production scale, on a code path with no BAM anywhere
in it, and it belongs in §2 as the headline rather than in a footnote — it is the only arm that
substantiates the plan's own claim that the scope is "every record-output file, not the uBAM
writer".

### 1.2 §2.4's `--clumpify` exclusion is wrong

§2.4 note 2 reasons: `--clump_only` needs `--output-format ubam` on uBAM input, and `--clumpify`
is refused *with* uBAM output, therefore "`--clumpify` is unreachable from a BAM input and
contributes no per-record-bail arm."

The second premise is about output format and does not constrain the input. Measured:

```
$ trim_galore --clumpify --cores 2 --preserve-tags CB -a AGATCGGAAGAGC --no_poly_g \
      -o D test_files/ubam_ws_qname_late.bam
Caused by: BAM read name contains whitespace: "readB with space"     (exit 1)
$ ls -l D  →  ubam_ws_qname_late_trimmed.fq   0 bytes
```

Tenth arm, per-record bail, zero-byte residue. (`--clump_only` on uBAM input without
`--output-format ubam` *is* refused, as the plan says — that half is right.)

It shares `parallel.rs:956` with `--cores 2`, so Option B covers it with no extra code; the
consequence is only that §2.4 and T3 are incomplete, and that a claim of unreachability is on
record which a future reader will rely on.

### 1.3 §2.4's paired-arm note misreads its own measurement

The note says the two `--paired` rows were refused by `open_paired_interleaved` because
"this fixture's records are not mate-adjacent with matching names". The actual message is:

```
Error: Paired uBAM 'test_files/ubam_ws_qname_late.bam': record 1 has invalid flag bits
(FREAD1=false, FREAD2=false, raw=0x0004).
```

An SE fixture used as paired input, rejected at record 1 — i.e. those two rows demonstrate only
that the writers are created before the first record is read. Real residue, real table rows,
wrong explanation, and the wrong conclusion drawn from them (see I3).

### 1.4 Site inventory: complete

`grep -n "FastqWriter::create\|BamWriter::create\|File::create" src/*.rs` returns exactly the
production sites §2.5 lists, at exactly those lines. Nothing missed. Two notes:

- The **test-only** construction sites are under-listed. §2.5 names `adapter.rs:762` and
  `fastq.rs:794`; there are also `parallel.rs:1308`, `:2024–2026`, `:2109–2111`,
  `specialty.rs:541`, `clump_only.rs:1113/1544/1606`, and seven in `bam.rs`. Every one of them
  will stop producing a file the moment `create` starts returning a pending output, so they are
  part of the diff, not a footnote.
- **`parallel.rs` does not move output handles into workers.** `File::create` at `:212–216` and
  `:956` runs on the main thread inside `std::thread::scope`; workers pre-compress into
  `Vec<u8>` and the main thread concatenates. Holding a `PendingOutput` beside each handle for
  the scope's lifetime is straightforward. Step 4 is feasible as written.

### 1.5 Interaction with the pre-flight and the namers — the plan is right, for a partly wrong reason

Keeping `.partial` out of the pre-flight is correct. §11's stated reason (candidate *order* is
load-bearing, and #423 lands next) checks out — `integration_output_collision.rs:1148` carries
"No filename assertion: candidate order decides which report is flagged". But the stronger reason
is that `.partial` paths are derived 1:1 from final paths, so they carry no collision information
the final paths do not: any two runs that collide on `<final>` collide on `<final>.partial` too.
Say that, and the exclusion stops looking like a test-driven compromise.

`the_site_inventory_is_complete` (tripwire:1140) asserts `main.rs` has exactly 11
`preflight_output_collisions(` call sites. Untouched by this plan. ✅

Naming: `.partial` appended to the whole file name is safe against every namer I checked.
`is_gzipped` reads only the final extension and is never asked about a temporary;
`strip_fastq_extensions` (`io.rs:641`) is only ever applied to *inputs*. §5 step 1's rule
("append to the whole file name") is the right one and the invariant is placed correctly.

**The tripwire trap.** `tests/integration_preflight_tripwire.rs` asserts `created ⊆ planned` in
one direction only (`compare`, T1, line ~426 — "over-planning is legitimate … `planned ⊆ created`
is false by design"). Two consequences:

- A `.partial` that survives a *successful* arm turns the tripwire red with
  `created but never planned: [...fq.partial]`, and the failure text instructs the reader to
  *"Add each to the candidate list at {site}"* — precisely the fix §11 wants to avoid. Give the
  tripwire a `.partial` clause so that failure says "a driver forgot to commit" instead.
- More seriously, the tripwire is **blind to a forgotten commit**: a missing output shrinks
  `created`, and a subset relation still holds. See C2.

### 1.6 Step 6's multi-file commit: stated behaviour is right

R1 commits, R2's rename fails → R1 published, R2 at its temporary, hard error naming both. Worse
than today only for a fault that requires the directory to become unwritable between the last
record write and the rename. No data is lost either way, and refusing to un-commit is the right
call. One scenario the plan does not mention comes out *better*: a SIGKILL between commit 1 and
commit 2 leaves a complete R1 and no R2, which a paired-end consumer fails on loudly — today the
same kill leaves two valid, differently-truncated FASTQs, which is the silent case.

### 1.7 §11's other three claims

- *"Nothing reads an output file while the run writing it is still going."* Holds, with a
  constraint the plan states for `--fastqc` but not for `--demux`: `main.rs:1474` calls
  `demux::demultiplex(&output_path, …)`, which re-opens the primary output with
  `FastqReader::open` (`demux.rs:211`). So the commit-before-consumer ordering has **two**
  post-loop consumers of the final name. In the SE serial driver both sit well after
  `main.rs:1333`, so the natural implementation is correct — but §4's edge-case table should name
  demux alongside FastQC, because it is the one that would break silently rather than with
  "file not found".
- *"`--clumpify` writes no intermediate."* ✅ `clump_only.rs:298–348` bins in memory, writes
  through one `BufWriter<File>`, and only ever stats the **input** afterwards.
- *"Reports are never partial."* ✅ for a mid-stream bail. Not for a write fault: the txt report
  at `main.rs:1431–1436` has no `flush()` and relies on `BufWriter::drop`, so an error writing it
  is discarded and the run still exits 0 (the JSON report at `:1450` does flush). Pre-existing,
  out of scope, note only.

---

## 2. Critical findings

### C1 — `Option<Box<dyn Write + Send>>` fixes the ordering and locks in silent error-swallowing at the publish moment

This is the answer to the question §5 step 2 raises and does not finish. The `Option` is
sufficient to get the *sequence* right. It is not sufficient to get the *outcome* right, because
`Box<dyn Write>` erases the type, and every fallible finaliser lives on the concrete type:

| Concrete writer | Fallible finaliser that exists | What `drop(box)` does instead |
|---|---|---|
| `BufWriter<File>` | `into_inner() -> Result<File, IntoInnerError>` | flush error discarded |
| `BufWriter<GzEncoder<File>>` | `into_inner()?` then `GzEncoder::finish() -> io::Result<File>` (write.rs:125) | trailer write error discarded (write.rs:163: `let _ = self.try_finish()`) |
| `ParCompress<Gzip, File>` | `ZWriter::finish() -> Result<W, GzpError>` (compress.rs:377) | `self.finish().unwrap()` — **panics** (compress.rs:391) |

So step 2 as written produces this: the disk fills while the gzip trailer is being written, `Drop`
swallows the error, `commit()` renames, and a **truncated file appears under the final name**.
That is the harm the step's own warning box exists to prevent, reached by a different route — and
T2 cannot see it, because T2 is a happy-path record count.

The same hole is already open in `BamWriter::finish` today and step 3 inherits it verbatim.
`inner: bam::io::Writer<bgzf::Writer<BufWriter<File>>>` (`bam.rs:526`): `try_finish()` propagates
the BGZF EOF write into the `BufWriter`, and the `BufWriter → File` flush happens in `Drop` after
`finish` has already returned `Ok(())`. A disk-full on the last uBAM block is exit 0 with a
truncated BAM, today. Step 3's "drop `self.inner` and then commit" preserves that exactly.

**Fix, inside B.** Make the finalisation fallible before the rename. Either hold a small private
enum over the three concrete writers, or keep one boxed trait object and give the trait a
consuming fallible finaliser:

```rust
trait FinishWrite: Write + Send {
    fn finish_boxed(self: Box<Self>) -> Result<()>;
}
```

with three impls that call the table's middle column. Then `FastqWriter::finish` is
`finalise (fallible) → commit → return`, no `Drop` in the success path at all, and `Drop` is left
holding only the abort path. Check the `'scope` lifetime on `ParCompress` when you do it — that is
the one place the shape may need an enum rather than a trait object.

This also makes gzp's `unwrap()` reachable only on the abort path, where a panic is harmless.

### C2 — the compile-time enforcement §11 relies on does not exist

§11 residual risk (b): "~20 mechanical sites, where a missed one is a silent hole — mitigated by
`finish()` consuming `self`, so the compiler enumerates them."

It does not. Adding `finish(self)` alongside the existing `pub fn flush(&mut self)`
(`fastq.rs:635`) leaves every current site compiling unchanged: `writer.flush()?; drop(writer);`
type-checks exactly as before, `PendingOutput::drop` then removes the temporary, and the run
**exits 0 having written no output file**. On a successful run. That is a worse failure than the
bug being fixed.

What is actually true, which is better than the plan realises but needs saying differently:

- Every production `FastqWriter` site does currently call `flush()` or sits in a container —
  `main.rs:1333/1627/1897`, `specialty.rs:69/113/334/431`, `demux.rs:254/256`,
  `clump_only.rs:348/519`. So **removing or privatising `FastqWriter::flush`** does turn all nine
  driver sites into compile errors. That is the lever, and step 2 must say so; adding `finish`
  without removing `flush` buys nothing.
- The raw-`File` sites (`parallel.rs:212–216`, `:956`, `clump_only.rs:298/440/444`) are converted
  by hand in step 4 anyway.
- **`demux` is the one site no API change reaches.** `writers: HashMap<String, FastqWriter>`
  (`demux.rs:191`) followed by `drop(writers)` (`:257`) compiles forever. It is also the widest
  commit set, exactly as §11(c) says.
- Renaming the constructor (`FastqWriter::create` → e.g. `create_pending`, same for `BamWriter`)
  is the only thing that forces a deliberate visit to **all** ~30 construction sites including
  the test-only ones. Cheap, since the plan touches them all regardless, and it is a one-time
  guarantee rather than a standing one.
- For regression, the honest mitigation is the existing suite, not the compiler: virtually every
  integration test and every arm of the CI validation matrix asserts on output *content*, so a
  forgotten commit goes red immediately. Say that instead. But note that the one test built
  specifically to catch "a writer fell out of the covered set" — the #414 tripwire — is
  structurally blind to it, because T1 only checks `created ⊆ planned`.

**Net effect on O1.** B does not remove the per-driver discipline the plan holds against A; it
replaces "remember to unlink on error" with "remember to commit on success". B still wins, for
three reasons the plan should state in these terms: forgetting is caught by the existing suite
under B and by nothing under A; B's failure mode is a missing file rather than a plausible one;
and B fixes the truncation hole, which A cannot.

### C3 — `--clumpify` is a tenth arm, and the plan says it is unreachable

§1.2 above. Correct the note, add the arm to §2.4, add a T3 case. No code consequence.

---

## 3. Important findings

### I1 — add the gzipped arm; it is the strongest evidence and the cheapest test

§1.1. `test_files/truncated.fq.gz` is committed, needs no fixture work, exercises the serial
`GzEncoder` path that step 2's hazard is about, and produces a valid 188,330-record residue.
Add it to §2.4 (and to §2's framing) and as a T3 arm — plus its `--cores 2` variant, which
covers `ParCompress`. Without it, no test in §7 exercises a gzipped residue at all, because
every §2.4 arm yields plain output.

### I2 — three existing statements become false and are outside step 9's scope

Step 9 calls `docs/.../guide/outputs.md` "the only user-visible surface of the change". These
already promise the opposite:

- `src/cli.rs:295` — `--passthrough` help text: *"On a mid-stream reader error (truncated or
  desynced passthrough), partial output files may remain on disk — re-run after fixing the
  input."* This is `--help`, i.e. user-visible.
- `docs/src/content/docs/modes/passthrough.md:71` — *"**Mid-stream errors leave partial output on
  disk.** … files already-written to disk are not rolled back."*
- `src/parallel.rs:90` and `:263` — the "v1 contract" doc comment and the in-loop comment, both
  asserting partial outputs may remain.

`docs/.../reference/changelog.md:829` says the same thing about multi-pair runs; that one is a
historical entry and must be left alone.

### I3 — O4 dissolves: the T6 fixture already exists

`test_files/ubam_ws_qname_paired.bam` (184 B) is a mate-adjacent 4-record fixture with the
offender as the first read of pair 2, already used by
`integration_ubam.rs::ws_qname_refused_via_interleaved_deinterleaver`, which pins "BAM record 3".
Measured, three consecutive runs, identical each time:

```
--paired                     → _val_1.fq 0 B + _val_2.fq 0 B, "BAM record 3: BAM read name
                               contains whitespace: \"pairB with space\""
--paired --output-format ubam → _val.bam 240 B
```

So T6 needs no new fixture and O4 needs no answer. One correction to T6's wording: on this arm
the bail *is* raised inside the de-interleaver (it converts records, so the whitespace check
fires there) — what the fixture proves is that the refusal is the per-record whitespace bail at
record 3 rather than a structural rejection at record 1. Assert the whitespace message and
`BAM record 3`, as the existing test already does, and add the residue-absence assertion.

### I4 — a symlinked output path changes behaviour; measured

Today's `File::create` writes *through* a symlink. Measured at `285afba`: with
`x_trimmed.fq` → `real/target.fq`, a successful SE run leaves the symlink intact and 18,931 bytes
in `real/target.fq`. Under B, `rename("x_trimmed.fq.partial", "x_trimmed.fq")` **replaces the
symlink with a regular file**, and the data lands in the source directory — which may be the small
volume the symlink existed to avoid, and is the one shape in which the "sibling, so never `EXDEV`"
argument does not hold in spirit.

Same class, same cause: a pre-existing **hard link** to the output is broken by the rename but
preserved by truncate-in-place; and a pre-existing output's **mode and ownership** are preserved
today and replaced by `0644 & ~umask` plus the invoking user under B, which matters in shared
group-writable pipeline directories.

Cheap fix that preserves today's semantics: if the final path is a symlink, resolve it once and
derive both the temporary and the rename target from the resolved path. Two lines, keeps the
sibling property, keeps cross-mount staging working. At minimum this needs a row in §4 and a
stated decision.

### I5 — `--demux` cannot commit "in a fixed order" as specified

§5 step 6 requires a fixed commit order. `demux.rs:191` is a `HashMap<String, FastqWriter>`, whose
iteration order is randomised per process, so `for (_, w) in writers.drain()` gives a
nondeterministic order and therefore a nondeterministic partial-commit error message. Sort by
barcode (or by output path) before committing. §11(c) already flags demux as the widest and least
exercised set; this is the concrete edge inside it, and it is worth one T3 case rather than the
acknowledged gap — `--demux` is in the CI validation matrix, so a forgotten commit there would go
red anyway, but the *ordering* would not.

### I6 — step 1 drops the parent-directory creation

`FastqWriter::create` currently does `create_dir_all` on the parent when it does not exist
(`fastq.rs:589–595`) before `File::create`. `PendingOutput::create(final_path)` as specified does
not, and the raw-`File` sites converted in step 4 never had it. Decide where it lives — inside
`PendingOutput::create` is the natural home now that it is the single choke point — or `-o
newdir/` regresses.

### I7 — §2.4's BAM byte counts are argv-dependent

I reproduced all nine arms; the four FASTQ counts match to the byte (103, 103, 0, 0+0), the BAM
counts do not (285 vs 359, 256 vs 321, 274 vs 365, 230 vs 303). Cause: `@PG CL:` embeds the full
command line (`bam.rs:662` `build_output_header`), so a BAM's size depends on how long the
invocation string was. Not an error in the plan, but no test may assert those numbers, and the
table should say they are illustrative.

### I8 — §6 is wrong about the cost: peak disk usage doubles on a re-run

§6 says "No extra copy … No extra memory", which is true, and stops there. It misses the one real
resource change: during a re-run over an existing output, the old final file and the new
`.partial` coexist until the rename, so peak disk usage for that output is **2×** where today the
`File::create` truncation frees the old blocks at open time. For a WGBS run writing a
200 GB `_val_1.fq.gz` on a volume with 300 GB free, today succeeds and B fails with `ENOSPC`.

This is inherent to write-then-publish and is not a reason to prefer A — but it is the honest
price of B at Buckberry scale, it belongs in §6, and "a previous good output survives a failed
re-run" and "peak footprint doubles on a re-run" are the same sentence read two ways. Felix should
see both halves when answering O1.

---

## 4. Assumptions

§8's nine are individually checkable and I checked all of them. Results:

1. Five bails are read-side, reach every output path — ✅, and now ten arms not nine (C3).
2. Reports never partial — ✅ for a bail; not for a write fault (§1.7).
3. `rename` atomic within a directory on both CI platforms; no Windows — ✅.
4. Temporary is a sibling of the final path — ✅ as a design statement, ⚠️ the symlink case makes
   "sibling of the final path" and "same filesystem as the data" two different things (I4).
5. Nothing reads an output while the run writes it — ✅, with `--demux` as an unlisted second
   post-loop consumer (§1.7).
6. Suffix is a compile-time literal — ✅ appropriate.
7. Overwriting a pre-existing `<final>.partial` acceptable — ✅. Residual, note-only: an *input*
   literally named `*.partial` sitting at a colliding path is outside the pre-flight's #383
   input-overwrite check, because `.partial` candidates are (rightly) not pre-flighted.
8. Cross-input rollback out of scope — ✅, and consistent with the existing per-file semantics.
9. No new runtime dependency; `tempfile` is dev-only — ✅ verified, `[dev-dependencies]` at
   `Cargo.toml:50`, `tempfile = "3"` at `:52`.

**Unstated assumptions worth promoting:**

- **Field declaration order in `FastqWriter`/`BamWriter` is load-bearing** for the abort path: the
  writer must drop before `pending`. Benign on POSIX either way (unlinking an open file just
  sends the trailer to a detached inode, and no final name appears), so it is a clarity
  invariant rather than a correctness one — one line above the struct, not in `finish`.
- **`panic = "abort"` is not set** (`[profile.release]`: `opt-level`, `lto`, `codegen-units`,
  `strip` only), so §4.4's "unwinding panic" coverage is real. Worth stating, because setting
  `panic = "abort"` later would silently downgrade the guarantee to "leaves a `.partial`".
- **A `.partial` on a successful run is always a bug.** That is the property that makes the
  tripwire clause in §1.5 self-policing rather than an allowance.

---

## 5. Efficiency

One `rename(2)` per output per run: agreed, unmeasurable. No per-record cost, nothing in the hot
loop, no allocation change, compression path untouched — all correct.

The unlink on the abort path is O(1) metadata regardless of residue size (I unlinked a 6.1 MB
residue; a 200 GB one is the same call).

The miss is I8: peak disk doubles on a re-run over an existing output.

One more, minor: `Drop for PendingOutput` ignoring `remove_file` errors is right (nothing useful
to do while unwinding), but it means a failed cleanup is invisible. Since the failure mode is a
stray `.partial` — the documented SIGKILL residue shape — that is the correct trade, and §5 step 1
already specifies it.

---

## 6. Alternatives

**Option A (unlink on error) is priced fairly and is the weaker call.** The plan's decisive claim
— A "structurally cannot" preserve a previous good output — is correct and I measured the harm it
cannot fix (1639 B → 103 B on a failed re-run). Any in-place scheme destroys the old file at open
time; the only escape is copying it first, which is worse than B at every scale. A remains the
right fallback if the diff is judged too wide, and the plan says so.

**One argument against B that §3 does not make, and should.** Preserving the previous good output
converts a loud failure into a *silent staleness* for any caller that ignores the exit code: today
such a caller consumes a truncated file, under A it finds no file at all (loudest), under B it
finds a complete-looking file from the previous run. This does not flip the recommendation —
today's behaviour is a truncated output *plus* stale reports, which is strictly worse than a
consistent set of stale outputs *and* stale reports — but it is the one place where A is arguably
better than B, and O1 should be decided with it visible.

**One argument for B that §3 does not make.** `ParCompress::drop` calls `self.finish().unwrap()`,
so a finalisation fault in the parallel-gzip path already panics. Under B that unwind removes the
temporary and no final name is ever created; under A the panic skips the unlink and leaves a
valid-looking output. B degrades correctly where A degrades into the reported bug. Same for
SIGKILL: `.partial` versus a plausible final name.

**Two options neither the plan nor I would take, recorded so they are not re-proposed:**

- *`O_TMPFILE` + `linkat`* (Linux) leaves no name at all until publication, which is strictly
  better than `.partial`, but it is Linux-only and would need the `.partial` path as a macOS
  fallback — two mechanisms for one guarantee. Not worth it.
- *`fsync` before rename.* Would upgrade the guarantee from "no torn file after a crash of the
  process" to "no torn file after a crash of the machine". Costs a real sync per output on a
  streaming trimmer. Out of scope, and the plan is right not to claim crash-consistency.

The two the plan rules out (lazy create, buffer-everything) are correctly ruled out.

---

## 7. Validation sufficiency

**T2 does not cover the hazard it is assigned to.** It covers the *ordering* slip — rename before
drop truncates the output, and a record count over a decompressed stream catches that. It cannot
cover C1, because C1 is a fault-path defect and T2 is a happy path. The mitigation for C1 is
structural (a fallible finaliser, with no `let _ =` on it) plus one cheap test that is actually
writable: create a writer, write a record, remove the output directory, call `finish()`, assert
the error names both paths. That is a unit test in `fastq.rs`, and it exercises the rename-failure
row of §4's edge-case table, which nothing in §7 currently touches.

**Arms untested after §7.** T1 covers SE uBAM-out. T3 covers FASTQ serial, `--cores 2`,
`--hardtrim5` ×2, `--clump_only` uBAM. T6 covers both paired arms. That leaves:

- `--dont_gzip` — identical bytes to FASTQ serial on this input; genuinely redundant. Fine.
- `--clumpify` — missing because the plan believes it unreachable (C3). Add.
- **any gzipped output** — missing structurally, not by oversight (I1). Add
  `truncated.fq.gz` serial + `--cores 2`.
- `--demux` — §11(c) admits it. With the HashMap ordering issue (I5) it is worth one case, not
  just an admission.
- `--passthrough` — not in §2.4 at all, yet it is the arm whose `--help` text this change
  invalidates (I2), and it has the widest simultaneous commit set after demux. One case.

**Would any test pass for the wrong reason?**

- **T1**: asserting "neither `…_trimmed.bam` nor any `*.partial` exists" passes vacuously if the
  run fails earlier than intended. Keep the `BAM record 2` stderr assertion, as the plan says —
  that is what makes it non-vacuous. ✅ as specified.
- **T3/T5**: same shape, and `--cores 2` is the one where "absent" versus "0 bytes" is the whole
  point, so assert `!path.exists()` and not `metadata().len() == 0`. The plan says this; hold to
  it.
- **T4** is the discriminating test and is specified loosely. Make it: run to success, hash the
  output, plant that file at the path the offending fixture will target, re-run, assert
  byte-identity. Without the planting step the two runs have different output names and the test
  proves nothing. (My probe did exactly this and is the measurement behind §1's verdict.)
- **T7**: "existing tests where every read is filtered out" — verify at least one such test
  actually asserts the file *exists*, rather than only asserting its record count is zero, or
  T7 is a no-op.
- **T8**: `cargo test --release` also means `debug_assert!` is compiled out, so any
  debug-assertion-based guard against a forgotten commit will not fire in CI. Do not rely on one.
- **T9**: naming the jobs rather than counting them is right. Both exist:
  `Validate vs Perl TrimGalore` (ci.yml:313) and `Validate uBAM input` (ci.yml:957).

**One gap in §7 with no case at all:** the rename-failure and partial-commit rows of §4's
edge-case table. Both are stated behaviours; neither is validated. The unit test above covers the
first cheaply; the second (commit #2 of N fails) needs fault injection and is fair to leave
unvalidated *if* the plan says so out loud.

---

## 8. Action items

### Critical

1. **C1 — replace drop-to-finalise with a fallible finaliser** in `FastqWriter::finish` and
   `BamWriter::finish`. `Option<Box<dyn Write + Send>>` gets the ordering right and makes the
   error from the file's final write unreachable; a disk-full at the trailer then publishes a
   truncated file under the final name. Use a consuming fallible finaliser per concrete writer
   (`BufWriter::into_inner`, `GzEncoder::finish`, `ZWriter::finish`) behind one trait or enum, and
   flush `BamWriter`'s `BufWriter<File>` explicitly before committing rather than in `Drop`.
2. **C2 — fix §11(b) and get the enforcement you thought you had.** The compiler does not
   enumerate the sites: `finish()` added *alongside* `flush()` leaves every current site
   compiling and silently producing no output on a successful run. Remove or privatise
   `FastqWriter::flush` (that does break all nine driver sites — verified), rename both
   constructors to force a visit to the ~30 construction sites including test-only ones, and
   restate the regression mitigation as "the existing content-asserting suite and the CI
   validation matrix", noting that the #414 tripwire is blind to a missing file. Call out
   `demux`'s `HashMap` + `drop(writers)` as the one site no API change reaches.
3. **C3 — correct §2.4's `--clumpify` claim** (measured reachable from uBAM input with FASTQ
   output; 0 B residue; tenth arm) and add the T3 case.

### Important

4. **I1 — add the gzipped arm** (`test_files/truncated.fq.gz`, serial and `--cores 2`) to §2.4,
   §2's framing, and T3. It is the only measurement of the gzip encoding, the only non-BAM arm,
   and by five orders of magnitude the most consequential face of the bug: a valid, `gunzip -t`
   clean, 188,330-record 6.1 MB output.
5. **I2 — widen step 9.** `src/cli.rs:295` (`--help`), `docs/.../modes/passthrough.md:71`, and the
   source contracts at `src/parallel.rs:90` and `:263` all promise partial output survives. Leave
   `changelog.md:829` alone.
6. **I3 — close O4 with no work.** `test_files/ubam_ws_qname_paired.bam` already is T6's fixture;
   measured to leave the residue deterministically. Also fix §2.4's paired-arm note (the actual
   error is invalid flag bits at record 1, not mate-adjacency).
7. **I4 — decide the symlink case** and put it in §4. Measured: today writes through the symlink;
   B replaces it. Same class for hard links and for a pre-existing output's mode/ownership.
   Resolving a symlinked final path once preserves today's semantics in two lines.
8. **I8 — correct §6.** Peak disk usage doubles on a re-run over an existing output; at WGBS scale
   that can turn a run that succeeds today into `ENOSPC`. It is the true price of B and belongs
   next to the property it is the flip side of.
9. **I5/I6 — two mechanical corrections to §5.** Commit `--demux`'s `HashMap` in a sorted order;
   decide where `create_dir_all` lives now that `PendingOutput::create` is the choke point.
10. **Validation — tighten T4 (plant the previous output at the colliding path, or it proves
    nothing), add the rename-failure unit test, and add `--clumpify`, `--passthrough` and
    `--demux` T3 arms.** State explicitly that the partial-commit row is not validated.

### Optional

11. Give the tripwire a `.partial` clause so a surviving temporary does not get "fixed" by adding
    it to a candidate list — the failure text currently instructs exactly that.
12. Add the stronger reason for keeping `.partial` out of the pre-flight: derived 1:1 from final
    paths, so it carries no collision information the final paths do not.
13. Record in §3 both arguments the plan omits: against B, that preserving the previous output
    turns a loud failure into silent staleness for an exit-code-ignoring caller; for B, that
    `ParCompress::drop`'s `unwrap()` and SIGKILL both degrade to `.partial` under B and to a
    plausible final name under A.
14. Name `--demux` beside `--fastqc` in §4's ordering row — it re-reads the primary output at
    `main.rs:1474`.
15. Note in §8 that field declaration order is the abort path's invariant and that
    `panic = "abort"` is not set (so §4.4's panic coverage is real, and setting it later would
    silently downgrade the guarantee).
16. Note-only, no issues per the filing bar: the txt trimming report at `main.rs:1431–1436` has no
    `flush()` and discards write errors while exiting 0; an input literally named `*.partial`
    is outside the #383 input-overwrite check; §2.4's BAM byte counts are argv-dependent via
    `@PG CL:` and must not be asserted. If any of these is worth acting on, they are one batch,
    not three.

---

## Appendix — probe log

All under `$TMPDIR/revb2/`, binary `target/release/trim_galore` @ `285afba`.

| Dir | Command shape | Result |
|---|---|---|
| `p1_*` | §2.4 arms 1–4 on `ubam_ws_qname_late.bam` | 285 B BAM / 103 B / 103 B / 0 B |
| `p2_*` | §2.4 arms 5–9 | 0+0 B / 256 B / 31 B / 274 B / 230 B; paired arms refused at record 1 for flag bits |
| `p3_paired_{1,2,3}` | `--paired` on `ubam_ws_qname_paired.bam`, ×3 | 0+0 B, "BAM record 3 … whitespace", identical each run |
| `p3_paired_ubam` | same, `--output-format ubam` | 240 B `_val.bam`, same message |
| `p4_rerun` | success run, plant output at the offender's path, re-run | 1639 B → **103 B** |
| `p5_clumpify` | `--clumpify --cores 2` on uBAM input | exit 1, whitespace bail, 0 B residue |
| `p5_clumponly_fq` | `--clump_only` on uBAM input, FASTQ out | refused at CLI (plan correct) |
| `p6_symlink` | SE run onto a symlinked output path | symlink intact, 18,931 B in the target |
| `p7_trunc_serial` | `truncated.fq.gz`, serial | 6,124,349 B, `gunzip -t` clean, 188,330 records |
| `p7_trunc_cores2` | same, `--cores 2` | 6,000,099 B, valid, 184,286 records |
