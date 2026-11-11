# PLAN_REVIEW_A — #428 "a mid-stream refusal must leave no output file"

**Reviewer:** A (independent; Reviewer B reviewing in parallel, no shared state)
**Plan:** `plans/08132026_partial-output-on-refusal/PLAN.md` v1
**Verified against:** main working tree on `dev` @ `285afba` (source read directly; my worktree was
pinned at `c6528e5`, ~13 k lines behind, so **no** line number below comes from the worktree)
**Binary used for probes:** `/Users/fkrueger/Github/TrimGalore/target/release/trim_galore` @ `285afba`

**Verdict on O1: Option B — but conditionally.** B is the right call and the reasons the plan gives
for it are true (I reproduced the truncation harm). But the plan has **mis-priced B by omitting its
real cost**: B introduces an ordering invariant that four existing in-run consumers of the output
already depend on, and the plan's Assumption 5 asserts that those consumers do not exist. Fix
Assumption 5 and §5 step 4's wording before implementing; do not implement B as written.

---

## 0. What I verified, and what checked out

Confirmed correct, by reading the source at `285afba`:

| Claim | Status |
|---|---|
| `tempfile` is a **dev**-dependency only (`Cargo.toml:52`, under `[dev-dependencies]` at `:50`) | ✅ |
| CI OS matrix is `os: [ubuntu-latest, macos-latest]`, no Windows (`ci.yml:28`) | ✅ |
| `main.rs:401-403` quote is exact ("The whole run is refused rather than decided per record, because a mid-stream refusal would leave a partial BAM behind.") | ✅ |
| `BamWriter::create` does `File::create` + `write_header` (`bam.rs:565`, `:570`) | ✅ |
| `FastqWriter::create` does `File::create` (`fastq.rs:597`) | ✅ |
| SE serial site reads `writer.flush()?; drop(writer);` (`main.rs:1333-1334`) | ✅ exact |
| `flate2::write::GzEncoder::drop` → `let _ = self.try_finish()` (`flate2-1.1.9/src/gz/write.rs:161-167`) | ✅ |
| `noodles_bgzf::Writer::drop` → `let _ = self.try_finish()` (`noodles-bgzf-0.47.0/src/io/writer.rs:231-240`) | ✅ |
| `BamWriter::finish(mut self)` already exists and calls `try_finish()` (`bam.rs:644-649`) | ✅ |
| **Site inventory (§2.5) is complete and every line number is right.** I enumerated every `File::create` / `FastqWriter::create` / `BamWriter::create` in `src/` and cross-checked each against its module's `#[cfg(test)]` boundary. Every production record-output site the plan lists exists at the stated line; every site I found that the plan omits is inside `mod tests` (`adapter.rs:762` > `:492`; `fastq.rs:794` > `:648`; `parallel.rs:1308/2024-2026/2109-2111` > `:1218`; `specialty.rs:541` > `:528`; `clump_only.rs:1113/1544/1606` > `:1104`; `format.rs:443-527` > `:428`; `bam.rs:1614+` > `:1139`; `demux.rs:337/378` > `:308`) | ✅ no missed site |
| Report sites are all post-loop (`main.rs:1431/1450/2215/2232/2598/2623`, `demux.rs:283`, `clump_only.rs:590`) | ✅ |
| `parallel.rs`'s output handles are created and written **on the main thread** inside the `thread::scope` closure (`:212-216`, written at `:238-253`), never moved into workers | ✅ holding an RAII guard across them is straightforward — the concern raised against this is unfounded |
| `rename` cannot fail `EXDEV`: a sibling is in the same directory, hence the same filesystem | ✅ (the §4 table already handles other rename failures) |
| §11's conclusion that `.partial` must stay out of the pre-flight | ✅ correct — see I6 for why the stated *reason* is not |

Reproduced empirically (my command lines differ from the plan's, so BAM byte counts differ — N1):

| Arm | Plan | Measured |
|---|---|---|
| SE FASTQ out | 103 B | **103 B** ✅ |
| SE uBAM out | 359 B | 364 B (N1) |
| `--cores 2` | 0 B | **0 B** ✅ |
| `--hardtrim5 10` | 31 B | **31 B** ✅ |
| `--clump_only --output-format ubam` | 303 B | 308 B (N1) |
| `--paired` interleaved, FASTQ | 0 B ×2 | **0 B ×2** ✅ |

**And I reproduced the harm that justifies Option B over Option A**, which is the plan's load-bearing
claim and worth restating as measured rather than argued:

```
good run  : ubam_test_trimmed.fq  = 1479 bytes  (complete)
failed re-run at the same output stem:
            ubam_test_trimmed.fq  =  103 bytes  (plausible short output)
```

A complete 1479-byte output was destroyed by a re-run that then failed, and replaced with a
well-formed 103-byte file. Option A cannot recover it — `File::create` truncates before the first
record is read. §1's secondary goal and T4 are correct, and this is the strongest single argument in
the plan.

---

## 1. Logic review

### C1 — Critical. Assumption 5 is false, and §5 step 4's wording walks straight into it.

Assumption 5 says:

> **Fixed:** no code reads an output file while the run that writes it is still going. FastQC runs
> after the writer is finished; `--clumpify` bins in memory and writes no intermediate.

Five production sites read a record-output file **by its final path, inside the same run**:

| Site | What it does with the final path |
|---|---|
| `main.rs:1459` | `fastqc::run(&output_path, …)` |
| `main.rs:2241` | `fastqc::run(&output_path, …)` (uBAM SE driver) |
| `main.rs:1474-1484` | `demux::demultiplex(&output_path, …)` — reads the trimmed output and fans out per-barcode files |
| `clump_only.rs:374`, `:847` | `fastqc::run(&output_path, …)` |
| `clump_only.rs:827`, `:1071` | **`stats.output_bytes = std::fs::metadata(&output_path)?.len();`** |

The last one is the sharpest: it is a `?`-propagating stat of the final path on the **success** path.
Under Option B, if the commit has not happened by then, the run fails with `ENOENT` **on a good
input**.

This matters because §5 step 4 says, of the raw-`File` sites:

> hold the `PendingOutput` for the function's lifetime, and commit after the concatenation loop
> returns `Ok`.

"Hold it for the function's lifetime" and "commit after the loop" are two different instructions, and
an implementer who follows the first will commit at the end of the function. For `clump_only.rs`'s SE
FASTQ path that is **after** `fastqc::run` at `:374`; for the two BAM paths it is after
`metadata(&output_path)` at `:827`/`:1071`. Today's code marks the correct point unambiguously —
`out.flush()?; drop(out);` at `clump_only.rs:347-348`, `writer.finish()?` at `clump_only.rs:822` —
and both sit *before* those consumers.

The plan's writer-level design happens to be right for exactly this reason: because the commit lives
inside `finish()`, and `finish()` is called where `drop` is today, the consumers keep seeing the
final name. But the plan never states that this is *why* it works, so nothing stops the raw-`File`
conversion from getting it wrong. The failure would be build-clean.

**Action:** rewrite Assumption 5 as an invariant rather than a non-existence claim — *"FastQC,
`--demux`, `write_clump_only_report` and `clump_only`'s `output_bytes` all consume the output by its
final path within the same run, so the commit must happen exactly where the writer is closed today,
never later"* — and restate step 4 as "commit at the `flush()`/`drop()`/`finish()` point that exists
today", naming `clump_only.rs:347-348` and `:822` explicitly.

### C2 — Critical (scope). §2.4 wrongly excludes `--clumpify`; there is a tenth arm.

§2.4 note 2 says:

> `--clumpify` is refused with uBAM output, so `--clumpify` is unreachable from a BAM input and
> contributes no per-record-bail arm.

The premise is right (`cli.rs:612-616` refuses `--clumpify` under `--output-format ubam`) but the
conclusion does not follow: `--clumpify` with **FASTQ** output on a **uBAM input** is not refused
anywhere. Measured:

```
trim_galore --clumpify --cores 2 -a AGATCGGAAGAGC --no_poly_g --preserve-tags CB \
            -o DIR test_files/ubam_ws_qname_late.bam
→ exit 1, "BAM record 2: BAM read name contains whitespace"
→ residue: ubam_ws_qname_late_trimmed.fq, 0 bytes
```

`--clumpify` requires `--cores >= 2`, which is why a naive probe reports a CLI error instead and
makes the arm easy to dismiss.

This is a per-record-bail arm with a residue, missing from the scope table that T3 keys off. As
written, the coverage audit will tick "every arm in §2.4 covered" with an arm missing. `--hardtrim3`
is in the same position, though less consequentially — I measured a 31 B residue, and §2.5 lists
`specialty.rs:97` so the fix covers it; only §2.4/T3 omit it.

The fix does cover `--clumpify` mechanically (`parallel.rs:956` / the clumpy PE path, both in step
4). The defect is in the scope statement, not the fix.

### C3 — Critical (test design). T6 and O4 rest on a misdiagnosis of the paired arm.

§2.4 note 1 says the two `--paired` rows "were refused by the **de-interleaver**, not by a per-record
bail", and infers that "pinning the five bails on a paired arm needs a new fixture".

The de-interleaver **is** the per-record bail site for paired-BAM input. `bam.rs:261` calls
`bam_record_to_fastq` inside the producer thread; the trimming loop never calls it for paired-BAM
input at all — it pulls already-converted `FastqRecord`s off a channel. So T6's stated objective, "A
paired arm bails per-record, not in the de-interleaver", is **unreachable by construction**.

Three consequences:

1. **O4 dissolves.** `test_files/ubam_ws_qname_paired.bam` already exists and is already exactly what
   T6 proposes: 4 records, mate-adjacent, names matching per pair, offender at record 3 (= pair 2 R1).
   `ws_qname_refused_via_interleaved_deinterleaver` (`integration_ubam.rs:386`) already drives it. I
   ran it both ways: FASTQ out → `_val_1.fq` + `_val_2.fq` at **0 B each**; uBAM out → `_val.bam` at
   333 B. No new fixture is needed. (O4 also miscounts — `test_files/` already holds 13 `.bam`
   fixtures, not eight.)

2. **The plan has the wrong cause for the 0-byte paired residue.** It is not the de-interleaver as
   such — it is batching. `BAM_READER_BATCH_SIZE = 4096` (`bam.rs:38`): the producer accumulates
   pairs in `batch_r1`/`batch_r2` and only sends at `len() >= 4096` (`bam.rs:287-301`); on the error
   path it `return`s **without flushing the partial batch** (`bam.rs:264-275`). So for any paired
   fixture with fewer than 4096 good pairs before the offender the residue is deterministically 0
   bytes — not a race, which is good news for test stability, but it also means:

3. **The dangerous paired shape is untested by any case in §7.** A paired run that has already
   written thousands of good pairs and *then* bails leaves a large, plausible, complete-looking
   output — the actual #428 harm. Every paired case in the plan lands in the 0-byte regime, where the
   residue is already suspicious. Covering the real shape needs > 4096 pairs before the offender.
   That is the one paired test worth having; T6 as specified tests file *creation*, not partial
   *content*.

T6 will still go red before the fix and green after (today the files exist at 0 B; the assertion is
"neither exists"), so it is not a broken test — its rationale and its fixture requirement are wrong.

### I1 — Important. `PendingOutput::commit` as sketched contradicts the §4 edge-case table.

§4's table says:

> Rename fails (permissions, read-only dir after the fact) | hard error naming both paths; **the
> temporary is left in place so the data is not lost**

The §5 sketch is `commit(mut self) -> Result<()>` with `Drop` removing the temporary "when
`!committed`". If `commit` sets `committed` only after a successful `rename`, a failed rename returns
`Err`, `self` drops, `!committed` holds, and **`Drop` deletes the temporary** — the data is lost,
which is the opposite of the stated behaviour. `PendingOutput` implements `Drop`, so `commit` cannot
destructure out of it; the flag is the only mechanism.

**Action:** state that the flag means "no longer owned by `Drop`" and is set **before** the `rename`
is attempted. Consider naming it `released` so the invariant is legible at the field.

### I2 — Important. The drop-order hazard is real and correctly identified; the completeness claim built on top of it is not earned.

The mechanism checks out at every layer, from the vendored sources:

- `flate2-1.1.9/src/gz/write.rs:161-167` — `Drop for GzEncoder` → `let _ = self.try_finish()`.
- `gzp-2.0.2/src/par/compress.rs:391-402` — `Drop for ParCompress` → `self.finish().unwrap()`.
- `noodles-bgzf-0.47.0/src/io/writer.rs:231-240` — `Drop for Writer` → `let _ = self.try_finish()`.

So yes: rename-then-drop publishes a truncated file. The `Option<Box<dyn Write + Send>>` field is
sufficient for the *ordering* — `Option::take()` needs only `&mut self`, which is what makes it work
inside a type that implements `Drop`, and Rust's declaration-order field drop means a
`writer`-before-`pending` declaration also gets the abandon path right. Nothing further is needed for
order to be correct.

Two things the plan should say and does not:

- **`ParCompress::drop` calls `.unwrap()`.** It is reachable in production — `demux.rs:197/203` and
  `specialty.rs:54/97/270/271/383/384` pass `cores` through to `FastqWriter::create`, so `--cores N`
  + gzip on those paths builds a `ParCompress`. A teardown failure therefore panics rather than
  erroring; on an unwinding abandon path that is a double panic → abort → `Drop for PendingOutput`
  never runs → a `.partial` is left. That still satisfies the plan's guarantee (the *final* name is
  never created), so it is not a defect in the fix — but §3's "covers every panic that unwinds" is
  one case narrower than stated, and the `.partial`-after-abort case belongs in the §4 table next to
  the SIGKILL row.

- **Every teardown swallows its errors** (`let _ = try_finish()`, and `BufWriter::drop`'s flush). So
  §4 step 3's "From this instant the output exists, complete" is not something the design delivers: a
  disk-full at trailer time still publishes a trailerless `.gz` or an EOF-less BAM under the final
  name. This is *not worse than today* (today the same bytes are already at the final name), so it is
  not a regression — but the plan should either drop the word "complete" or take the opportunity,
  since B is touching exactly this code. See A2.

### I3 — Important. Four existing places document the contract this fix reverses; step 9 updates none of them.

Step 9 is "One sentence in `docs/.../guide/outputs.md`" and calls that "the only user-visible surface
of the change". It is not:

| Location | Text that becomes false |
|---|---|
| **`src/cli.rs:295`** | `--passthrough` help: "On a mid-stream reader error (truncated or desynced passthrough), partial output files may remain on disk — re-run after fixing the input." — this is **`--help` output**, i.e. user-visible |
| `docs/src/content/docs/modes/passthrough.md:71` | "**Mid-stream errors leave partial output on disk.** … files already-written to disk are not rolled back" |
| `src/parallel.rs:90` | "partial outputs may still exist on disk on mid-stream error — v1 contract" |
| `src/parallel.rs:263` | "Partial output files may remain on disk (v1 contract; Cli --help documents this)" |

The last is a pointer to the first, so the code itself tells you the `--help` text exists.

Leaving the archived changelog entry at `docs/src/content/docs/reference/changelog.md:829` ("pair K
may have a partial output file") alone is correct — it describes a shipped version and the synced
changelog page is not to be edited.

This is the "corrective copy contains the defect it corrects" shape: the fix removes a documented
behaviour, the plan documents the removal in one new place, and four statements of the old behaviour
survive — one of them in `--help`.

### I4 — Important. Step 6 needs one more constraint: nothing fallible between the first and last commit.

Step 6's stated failure behaviour is **right**, and better than today's. If R2's rename fails after
R1's succeeded you get a complete `_val_1.fq`, a `_val_2.fq.partial`, and a hard error naming both.
Today the equivalent fault (a swallowed flush error at close) gives you a complete `_val_1.fq` and a
*silently truncated* `_val_2.fq` at its final name. Loud-and-half-published beats
silent-and-looks-complete, and refusing to un-commit is correct — deleting a successfully-published
file on an I/O fault would be a second data loss. I could not construct a scenario where step 6
leaves a state worse than today.

What the plan does not say is the invariant that keeps the window minimal: **no fallible operation
may sit between the first commit and the last.** "Commit them in a fixed order after the loop"
permits an implementation that interleaves a report write or a stats computation between two commits,
widening the non-atomic window for no reason.

### I5 — Important (O2). `.partial` is visible to `*` globs and to workflow output collectors.

O2 weighs `.partial` against `.tg-partial` on legibility. Neither addresses the operational risk,
which is that a leftover temporary sits in the output directory under a name that `*` matches.
TrimGalore's dominant consumer is nf-core pipelines, where the output directory is a Nextflow task
workdir and outputs are collected by glob. A stale `sample_trimmed.fq.gz.partial` can be picked up by
a permissive `path("*.fq.gz*")`-style declaration or by a subsequent `trim_galore out/*`, and
`strip_fastq_extensions` will fall through to `file_stem()` on it and name a second-order output from
it.

A **dot-prefixed sibling** — `parent.join(format!(".{}.partial", file_name))` — costs nothing, is the
conventional spelling for exactly this, and is invisible to `*` globs and to `publishDir` patterns.
The trade-off is worth stating honestly: a dot-prefixed leftover is also invisible to the *user*, so
a SIGKILLed run leaves silent junk rather than an obvious breadcrumb. My reading is that the glob
risk is sharper, because it can silently enter a pipeline's output set, whereas hidden junk only
wastes disk.

Appending to the **whole** file name rather than inserting before the extension is right either way,
and the plan's reason for it — never colliding with an extension-stripping namer — is correct. I
checked `io.rs`: `is_gzipped` (`:33`) is extension-based, and every namer
(`single_end_output_name` `:278` … `clumped_paired_bam_output_name` `:548`) derives from the **input**
path, so nothing re-derives from a `.partial`; the writers take `gzip` as an explicit parameter
rather than sniffing the output name. No namer misclassification exists in either spelling.

### I6 — Important. §11's pre-flight reasoning reaches the right conclusion for the wrong reason, and the tripwire needs one line.

§11 says `.partial` names are kept out of the pre-flight "because adding candidates would reorder the
lists that `tests/integration_output_collision.rs` pins".

The conclusion is right. The reason is not, in two ways:

1. **The substantive reason is redundancy, not test churn.** `final → final + ".partial"` is
   injective, and injective under `collision_key`'s case-folding too, so if no two finals collide no
   two `.partial`s can. Adding them buys zero detection power. Stated as "it would churn a pinned
   test", the exclusion looks like a workaround, and a future maintainer may "fix" it by adding the
   candidates and re-pinning — churn for nothing. Say the real reason. (This also matters for the
   #423 hand-off, whose constraint is that candidate *order* is load-bearing: "we chose not to add
   candidates because they add no information" survives a reordering; "we chose not to because a test
   pins the order" invites revisiting it the moment that test is re-pinned.)

2. **The pinned test is `tests/integration_preflight_tripwire.rs`, not
   `integration_output_collision.rs`.** The tripwire (#414, `1fc00f6`) is the one that consumes
   `DUMP_PLANNED_ENV` / `DUMP_PLANNED_PREFIX` (`:63`, `:175`, `:190`) and asserts `created ⊆ planned`
   (`:426-433`, T1). `integration_output_collision.rs` uses neither symbol.

I checked whether Option B breaks the tripwire, since that was the real question: **it does not.**
Every case asserts exit 0 (`:285-289`) and snapshots after the process exits (`:291`), so on a
successful run the `.partial` has been renamed away and never appears in `created`. The `unchanged()`
hash check (`:315-325`) is likewise untouched, since a create-then-remove does not modify a
pre-existing file.

But there is a second-order point worth one line in the tripwire's own header. Option B introduces
~20 file-creation paths that are deliberately unplanned and that the tripwire is *structurally blind
to*, because they exist only between create and commit. If the fix ever fails to commit on a success
path, the tripwire goes red — correctly — with a message that says *"The fix is to add the missing
path to that arm's candidate list in `src/main.rs`"*, which for a `.partial` is precisely the wrong
fix, and the header already warns that widening the allowance is "how this test dies".

**Action:** add `.partial` to the tripwire's "Known-unplanned writes" section with the injectivity
argument and one clause saying a red `.partial` means a missed commit, not a missing candidate.

---

## 2. Assumptions

Walking §8 item by item.

| # | Plan's claim | Assessment |
|---|---|---|
| 1 | The five bails are read-side and reach every output path | **Holds**, and I measured a tenth arm the plan excludes (C2) |
| 2 | Reports are written only after the loop succeeds, so they are never partial | **Holds for the bail class**; slightly overstated in general — a disk-full *during* `report::write_report_header` leaves a truncated report at its final name (`main.rs:1431` onward). Reports being out of scope is defensible (a footerless text report is self-evidently truncated, and it is not the data), but say "never partial *from a mid-stream bail*" rather than "never partial" |
| 3 | `rename` atomic within a directory on both CI platforms; no Windows target | **Verified** (`ci.yml:28`). Unlink-while-open and rename-while-open are indeed both safe on the matrix |
| 4 | The temporary is a sibling, never in `$TMPDIR` | **Holds**, and this is what makes the `EXDEV` claim sound |
| 5 | No code reads an output file while the run that writes it is still going | **FALSE — see C1.** This is the one assumption that changes the implementation |
| 6 | Suffix is a compile-time literal, no flag | Agreed. See I5 on which literal |
| 7 | Overwriting a pre-existing `<final>.partial` is acceptable | **Holds.** `<final>` is overwritten unconditionally today, so this is strictly less surprising |
| 8 | Cross-input rollback out of scope | Agreed, and consistent with the existing per-pair semantics documented in the changelog |
| 9 | No new runtime dependency; `tempfile` is dev-only | **Verified** (`Cargo.toml:50`/`:52`). The derived-name approach means nothing needs promoting |

### Implicit assumptions the plan does not surface

- **`FastqWriter` is only ever constructed from a path.** True today — `FastqWriter { writer }` is
  built only inside `create`, and there is no `new`/stdout constructor. Option B's design depends on
  it: a writer over a non-file sink could not hold a `PendingOutput`. Worth one line, because it is a
  constraint on future work rather than an observation about current code.
- **Field drop order.** The abandon path relies on the inner writer dropping before `PendingOutput`.
  On POSIX the reverse order is harmless (writes land on the unlinked inode), so this is not a
  correctness risk — but it is a free invariant to state at the struct.
- **Concurrency.** Two runs writing the same output path concurrently corrupt each other today
  (interleaved writes to one final name) and under B (interleaved writes to one `.partial`, then a
  racing rename). B is *strictly better* for a concurrent **reader**: today `File::create` truncates
  the file in place, so a reader watching the old output sees it shrink to zero; under B the old
  inode survives until the atomic rename. Worth one line in §4, because it is a genuine second
  benefit of B that the plan does not claim.

---

## 3. Efficiency

§6 is correct and I have nothing to add on cost: one `rename(2)` per output file per run is a
directory metadata operation, invisible next to trimming; the sibling placement means no
copy-and-delete fallback; nothing in the per-record path changes; no new allocation.

One efficiency *opportunity* the plan does not mention, in the opposite direction: on the abandon
path, dropping the inner writer now deflates and writes a full gzip/BGZF trailer into a file that is
about to be unlinked. For `ParCompress` that also means joining compressor threads. It is pointless
work on an error path and it is the right trade (skipping it would require a non-`Drop` teardown
path, i.e. exactly the enum refactor in A2), but on a large `--cores N` run the refusal will now take
noticeably longer to return than it does today. Not a defect; worth a sentence so nobody files it as
one.

---

## 4. Alternatives

### A1 — Option A, honestly re-priced

The plan prices A as "~40 lines smaller, fixes only the reported symptom". That undersells A in one
respect and oversells it in another.

**In A's favour, which the plan does not say:** C1, I1, I2, I5 and I6 all *disappear* under Option A.
A introduces no ordering invariant, no drop-before-rename hazard, no new residue class, no
tripwire blind spot, and no writer API change. A's risk profile is genuinely much flatter, and the
plan's framing ("the honest fallback if the diff is judged too wide") undersells that as a
*correctness* argument rather than a diff-size one.

**Against A, decisively:** A does not fix the truncation hole, which I measured. Worse, A converts
"silently wrong file" into "no file at all" — the previous good 1479-byte output is deleted and
nothing replaces it. That is more honest and more destructive at the same time. And A leaves #408's
constraint at `main.rs:402` standing, because a panic still leaves a partial BAM, so the stated
secondary benefit of the whole exercise is not delivered.

There is a second argument for B over A that the plan never makes and that I think is stronger than
the one it leads with: **B is recoverable, A is not.** A 500 GB run that fails at 99% leaves, under
B, a `.partial` the user can rename by hand and salvage; under A the same run leaves nothing. For a
trimmer whose whole performance story is Buckberry-scale files, that asymmetry is worth stating in
§3.

**Conclusion: B.** But note the shape of the risk — every one of my Critical findings is a cost of B
specifically. B is right *and* it needs the corrections in C1/C2/C3 before it is safe to implement.

### A2 — Replace `Box<dyn Write + Send>` with an enum, so `finish()` can propagate teardown errors

This is a variant of step 2 rather than an alternative to B, and it is the only way to make §4 step
3's "complete" claim true. Today the trailer errors are unreachable through `Box<dyn Write>`:
`GzEncoder::try_finish` and `ParCompress::finish` both exist as inherent methods, but the trait
object erases them, so `finish()` can only drop the box and hope.

```rust
enum Sink {
    Plain(BufWriter<File>),
    Gz(BufWriter<GzEncoder<File>>),
    ParGz(ParCompress<'static, Gzip, File>),
}
```

`finish()` then matches and calls the concrete teardown, propagating the error — and as a bonus it
sidesteps `ParCompress::drop`'s `.unwrap()` panic entirely on the success path. Cost: an `impl Write
for Sink` with a three-arm match per `write_record`, which replaces a vtable indirection with a
branch — likely neutral or slightly faster, but it does touch the hot path, so it needs the
`scripts/benchmark.sh` numbers rather than an assertion.

My recommendation: **do not fold this into #428.** It is a real improvement and it is orthogonal to
the residue question. The right move for this PR is to delete the word "complete" from §4 step 3 and
say what is actually guaranteed — *the final name appears only after the writer has been closed* —
then leave the enum for a follow-up if Felix wants error-checked teardown. If it stays in scope,
T2 must be extended to the disk-full case, which needs a size-limited filesystem in the test and is
a meaningful extra cost.

### A3 — Ruled out, but worth naming so nobody proposes it

`O_TMPFILE` + `linkat(AT_EMPTY_PATH)` is the only design that gets SIGKILL residue to genuinely
zero: the temporary is never linked into the directory at all, so a killed process leaves nothing,
not even a `.partial`. It is unusable here — Linux-only (no macOS `O_TMPFILE`, and the CI matrix
includes `macos-latest`), and `linkat(AT_EMPTY_PATH)` needs `CAP_DAC_READ_SEARCH` on many kernels.
The plan's "Ruled out" list has the two weak candidates and not this one; adding it in a clause
forecloses the obvious "why not make the residue actually zero?" question.

### A4 — De-risking split, if one PR is judged too wide

Not "A instead of B", but **B in two PRs: the writers first, the raw-`File` sites second.**
`FastqWriter` and `BamWriter` are where the commit point already coincides with today's teardown, so
C1 cannot bite there. `parallel.rs:212-216/956` and `clump_only.rs:298/440/444` are exactly where the
commit point is ambiguous and where the four in-run consumers live. Splitting leaves `--cores N` PE
and `--clump_only` FASTQ still leaking residue in the interim, which undercuts the "retires the
class" framing — so I would not recommend it, but it is the right split if a split is wanted.

---

## 5. Validation sufficiency (§7)

### Is T2 sufficient against the drop-order hazard?

**Almost, and the gap is nameable.** T2 as written is "a success-path run whose output is gzipped,
decompressed and record-counted end to end; plus one uBAM run read back with `assert_ubam_eq`".

That catches the *single-threaded* gzip trailer (`GzEncoder`) and the BGZF EOF marker. It does **not**
cover `gzp`'s `ParCompress`, which is a different teardown with a different failure mode (it joins
compressor threads and flushes their blocks, and it panics rather than erroring). `ParCompress` is
only reachable via `FastqWriter::create` with `gzip && cores > 1`, i.e. `--demux` and the specialty
modes — precisely the paths §11 already admits are least exercised.

**Action:** T2 needs a third arm — a `--hardtrim5` (or `--clock`) run with `--cores 2` and gzip
output, decompressed and record-counted. Without it, a drop-order slip in the parallel-compressor
case ships truncated output on the exact paths the plan says are thinnest.

A second, subtler sufficiency point: T2 counts records. A gzip stream that is missing only its
8-byte CRC/ISIZE trailer will still decompress all its records under most decoders, and `flate2`'s
reader may or may not error depending on how the stream ends. **Record count alone can pass on a
trailerless file.** T2 should assert on `gzip -t` (or an explicit `MultiGzDecoder` read-to-end that
propagates the trailer error) in addition to the count, otherwise it can pass for the wrong reason —
which is the specific risk given that a lost trailer is exactly the failure it exists to catch.

### Which of the nine (ten) measured arms end up untested?

| §2.4 arm | Covered by |
|---|---|
| `--output-format ubam` | T1 ✅ |
| default FASTQ out | T3 ✅ |
| `--dont_gzip` | — (identical bytes to the default arm; acceptable) |
| `--cores 2` | T3 / T5 ✅ |
| `--paired` interleaved | T6 ✅ (see C3 on its rationale) |
| `--paired --output-format ubam` | T6 ✅ |
| `--hardtrim5` FASTQ | T3 ✅ |
| `--hardtrim5` uBAM | T3 ✅ |
| `--clump_only --output-format ubam` | T3 ✅ |
| **`--clumpify --cores 2` (the tenth arm, C2)** | **— untested** |
| `--hardtrim3` (measured by me, absent from §2.4) | — untested |
| `--demux` per-barcode fan-out | — untested, admitted in §11 |
| `--clock` / `--implicon` (§2.5 sites `270/271`, `383/384`) | — untested |
| `--passthrough` / `--retain_unpaired` (§2.5 sites `1597`, `1603/1604`) | — untested |

So §11's "`--demux` is the one gap" is **not** the only gap: `--clumpify`, `--hardtrim3`, clock,
implicon, passthrough and unpaired are all uncovered. Most are the same code shape as a covered arm
and I would not add a case for each — but the plan should say which ones it is consciously leaving to
`finish()`-consumes-`self` plus the compiler, rather than implying the list is `--demux` only.

`--demux` deserves better than "least exercised", though, because it is qualitatively different: it
is the only stage that **reads a committed output and then produces more outputs** (`main.rs:1474`
→ `demux.rs:197/203`). That makes it the arm where a commit-ordering slip is most likely and hardest
to reason about. My measurement: `--demux` on a uBAM input reaches the record-2 bail and leaves a
103 B `_trimmed.fq` residue at the trimming stage, so it is a live arm, not a theoretical one. It
should get one T3-style case.

### Would any test pass for the wrong reason?

- **T2** — yes, on the trailer, as above. Fix by asserting decoder-level integrity, not record count.
- **T3** — assert `!path.exists()` for the final name **and** that no `*.partial` remains. The plan
  says this; the trap is asserting only the first, which today's `--cores 2` arm would already fail
  and therefore looks like a working test. Both halves matter.
- **T4** — the discriminating test. It must re-run at the **same output path** (my probe K did this
  by copying the offending fixture to the successful run's stem, since the output name derives from
  the input name). Written naively — two different input names — it compares two unrelated files and
  passes vacuously against Option A as well. Worth spelling out in the plan, because it is the one
  test whose whole purpose is to distinguish B from A.
- **T7** — "existing tests where every read is filtered out". Confirm at least one such test exists
  and asserts the file's *presence*, not just its content; otherwise T7 is a claim about coverage
  rather than a validation. This is the guard against the fix over-shooting into "leaves no file when
  the output is legitimately empty", which is the one way this change could break a green pipeline.
- **T8** — the "684" baseline: I did not verify it (running the suite was outside my budget and would
  have risked the timeout that cost a previous review round three of four agents). Treat it as
  unconfirmed by this review.
- **T9** — naming the two jobs rather than counting them is right, and matches the lesson from the
  earlier CI-verification incidents.

### One missing validation

There is no test that `--fastqc` still produces a correctly-**named** artefact after the change. This
is the concrete user-visible symptom of a commit-ordering slip (C1): FastQC on
`x_trimmed.fq.gz.partial` derives its report name from the path it is given, so you would get
`x_trimmed.fq_fastqc.html` or a hard error instead of `x_trimmed_fastqc.html`. The
`--fastqc`-active arms of `integration_preflight_tripwire.rs` cover this incidentally via T1, so
T8 catches it — but the plan should say so rather than resting on the §4 table's "`--fastqc` |
unaffected" row, which is an assertion, not a check.

---

## 6. Notes (not for filing — no user-visible symptom, per the filing bar)

- **N1 — §2.4's BAM byte counts are not reproducible.** The `@PG CL:` line embeds the full command
  line, so the same fixture and flags give 359 B or 364 B depending on the length of the `-o` path. I
  got 364/308/333 where the plan reports 359/303/321. Nothing in the plan depends on the numbers
  (T1/T3 assert absence), but "**359 B**" in bold reads as a stable measurement. Say "≈360 B (varies
  with the command line, which the `@PG CL:` line embeds)".
- **N2 — `--hardtrim3` is absent from §2.4.** Measured 31 B residue. §2.5 lists `specialty.rs:97`, so
  the fix covers it.
- **N3 — O4 miscounts the fixture set.** `test_files/` holds 13 `.bam` fixtures, not eight, so "a
  ninth BAM fixture" understates the existing set — and per C3 no new fixture is needed at all.
- **N4 — `the_site_inventory_is_complete`** in `integration_preflight_tripwire.rs` pins the number of
  `preflight_output_collisions` call sites (the header says "catches a 12th call site"). Option B adds
  no pre-flight call sites, so it is unaffected — noting it only so the implementer does not go
  looking when the tripwire file comes up.
- **N5 — §9's `--gzip` no-op check was worth doing and its conclusion is right.** `cli.rs:425` →
  `cli.rs:1102` deprecation notice only; output compression comes from `is_gzipped` on the input plus
  `--dont_gzip`. Nothing for the fix to do there.

---

## 7. Action items

### Critical — resolve before implementation

1. **C1: rewrite Assumption 5 and §5 step 4.** Assumption 5's non-existence claim is false. Five
   sites consume the output by its final path inside the same run: `fastqc::run`
   (`main.rs:1459`, `:2241`, `clump_only.rs:374`, `:847`), `demux::demultiplex` (`main.rs:1474`),
   `write_clump_only_report`, and `std::fs::metadata(&output_path)?` (`clump_only.rs:827`, `:1071`).
   Replace it with the invariant *"the commit happens exactly where the writer is closed today,
   never later"*, and restate step 4 to name `clump_only.rs:347-348` and `:822` as the commit points
   instead of "for the function's lifetime". The `metadata` sites are `?`-propagating on the success
   path, so getting this wrong fails a good run with `ENOENT`.
2. **C2: add the `--clumpify` arm to §2.4 and T3.** `--clumpify --cores 2` on uBAM input with FASTQ
   output is *not* refused, reaches the record-2 bail, and leaves a 0-byte `_trimmed.fq` (measured).
   §2.4 note 2's "unreachable from a BAM input" is wrong. Add `--hardtrim3` to the table while there.
3. **C3: rewrite T6 and close O4.** The de-interleaver *is* the per-record bail site for paired-BAM
   input (`bam.rs:261`), so T6's stated objective is unreachable. `ubam_ws_qname_paired.bam` already
   is the fixture T6 proposes — **O4 needs no decision, only deletion**. Separately, add the one
   paired case that actually matters: > 4096 pairs before the offender
   (`BAM_READER_BATCH_SIZE`, `bam.rs:38`), which is the only way to get a paired residue that is
   large and plausible rather than 0 bytes.

### Important — fold into the plan

4. **I1:** state that `PendingOutput`'s flag is set **before** the `rename` is attempted, so a failed
   rename leaves the temporary as §4's table promises. Consider `released` over `committed`.
5. **I3:** add the four contract-reversal sites to step 9 — `cli.rs:295` (**`--help`**, user-visible),
   `docs/.../modes/passthrough.md:71`, `parallel.rs:90`, `parallel.rs:263`. Leave the archived
   changelog entry alone.
6. **§7/T2:** add a `--cores 2` + gzip arm (specialty or `--demux`) so `gzp::ParCompress`'s teardown
   is covered, and assert decoder-level integrity rather than record count alone — a trailerless gzip
   still yields the right record count under most decoders, so T2 as written can pass for the wrong
   reason on the exact failure it exists to catch.
7. **I2:** drop "complete" from §4 step 3 (teardown errors are swallowed at every layer, so the claim
   is not delivered), and add the `.partial`-after-abort row to the §4 table — `ParCompress::drop`'s
   `.unwrap()` makes a double-panic abort reachable, which skips `Drop for PendingOutput`.
8. **I4:** add "no fallible operation between the first commit and the last" to step 6.
9. **I5 (O2):** decide the suffix on the glob-visibility argument, not legibility. My recommendation
   is a dot-prefixed sibling, `.{file_name}.partial`, so `*` globs and Nextflow `publishDir` patterns
   never see a leftover.
10. **I6:** fix §11 to give the injectivity reason (not test churn) for keeping `.partial` out of the
    pre-flight, correct the test name to `integration_preflight_tripwire.rs`, and add one line to that
    file's "Known-unplanned writes" header so a red `.partial` is not "fixed" by adding a candidate.
11. **§7:** say which uncovered arms are conscious omissions. §11's "`--demux` is the gap" is
    incomplete — `--clumpify`, `--hardtrim3`, `--clock`, `--implicon`, `--passthrough` and
    `--retain_unpaired` are all uncovered too. Give `--demux` a T3-style case regardless: it is the
    only stage that reads a committed output and then writes more (measured: reaches the bail, leaves
    103 B).
12. **T4:** specify that the re-run must target the **same output path** (same input filename), or it
    compares unrelated files and passes under Option A too.

### Optional

13. **A1:** add the recoverability argument to §3 — B leaves a salvageable `.partial` where A deletes
    the data outright. At Buckberry scale that is a stronger argument than the one §3 leads with, and
    it is the answer to "why not just unlink?".
14. **A3:** add `O_TMPFILE` + `linkat` to the "Ruled out" list in one clause (Linux-only; the matrix
    includes macOS) so the "why isn't the residue actually zero?" question is foreclosed.
15. **A2:** note the enum-instead-of-`Box<dyn Write>` refactor as a follow-up that would make teardown
    errors propagable, explicitly **out of scope** for #428.
16. **§4:** claim the second benefit of B that the plan misses — a concurrent reader of the previous
    output sees it survive intact until the atomic rename, where today `File::create` truncates it in
    place.
17. **§8.2:** narrow to "never partial *from a mid-stream bail*"; a disk-full during
    `report::write_report_header` does leave a truncated report.
18. **§6:** one sentence that a refusal now costs a full trailer write (and a compressor-thread join
    on `--cores N`) into a file about to be unlinked, so nobody files the slower error path as a
    regression.
19. Add the "`FastqWriter` is only ever constructed from a path" constraint to §8 — it is what Option
    B depends on, and it is a constraint on future work rather than an observation.

### Not recommended

- **Filing anything.** N1–N5 are real and inconsequential; none has a user-visible symptom. C2 is the
  only finding that describes a live user-facing defect, and it is a scope correction inside #428, not
  a separate issue. Nothing here meets the bar.
