# CODE REVIEW A — #428, commit `0a6f3a3`

**Reviewer:** A (independent; no coordination with Reviewer B)
**Under review:** `0a6f3a3` on `fix/428-no-partial-output-on-refusal`, parent `285afba`
**Plan:** `PLAN.md` v2, including §12 implementation notes
**Verified in-tree:** `cargo test --release` → **701 passed, 0 failed** (macOS, full log, not tailed);
`cargo clippy --all-targets --release -- -D warnings` → exit 0, zero warnings.

---

## Summary

The design is right and the implementation matches it. `PendingOutput` is a small, correct RAII
guard; the `take(pending)` → `drop(self)` → `commit()` order in both `FastqWriter::finish` and
`BamWriter::finish` does what its comments claim; the commit-point invariant holds at every one of
the four in-run consumers, on the success path and on every error arm I could reach; and no output
byte changes, so the `validation` job's byte-identity invariant is intact.

Nothing here is a Critical. One user-facing docs sentence makes a data-integrity claim that the
implementation deliberately does not deliver, and should be corrected before merge. Beyond that, my
substantive findings are about **what the new tests do not pin**, and two of them I proved by
mutation rather than inferred:

- Setting `released = true` *after* the rename instead of before — the plan's own safety-critical
  rule — leaves **701/701 green**.
- Replacing `w.finish()?` with `drop(w)` for the two `--retain_unpaired` writers in `run_paired`
  leaves the suite green; the same mutation on clock, implicon and passthrough turns **8 tests red**.

Both were run in a worktree pinned to `0a6f3a3` with a private `CARGO_TARGET_DIR` (cold 1m36s build,
204 crates — the mutants were compiled, not cached). A positive control (making `Drop` leak the
temporary) turned **12 tests red**, so the green results above are meaningful rather than a harness
that cannot fail.

**Verdict: safe to merge once H1 is fixed** (a one-clause docs edit). H2, M1, M3 and M4 are each a
few lines of test and I would take them in this PR, but none of them blocks: every code path I could
reach behaves as specified.

---

## Errors / logic

### H1 — the docs promise that a full disk leaves nothing behind. It does not. (High)

`docs/src/content/docs/guide/outputs.md:124`:

> …so a refusal — a malformed read name a million records in, a truncated input, **a full disk** —
> leaves nothing behind…

A full disk is the one case in the list this change does not cover, and it fails in the worst
direction. ENOSPC *mid-stream* is fine — `write_all` errors, `?` propagates, nothing is published.
ENOSPC *at trailer time* is swallowed by `GzEncoder::drop` / `bgzf::Writer::drop`
(`let _ = try_finish()`), so `finish()` returns `Ok(())`, `commit()` runs, and a truncated file is
published under the final name **with exit 0**. That is exactly the hazard PLAN §5 step 2(b) labels
"structurally invisible to any success-path test" and §10 O5 defers.

The code is honest about this. `FastqWriter::finish`'s doc comment says "not that the bytes are
complete", and §4 item 4 weakened the guarantee from "complete" to "closed" precisely so this claim
would not be made. The docs then re-strengthen it, and pick the deferred case as the headline
example. A user who reads that sentence will trust a run that filled the disk.

**Fix** — `docs/src/content/docs/guide/outputs.md:124`, delete the clause:

```
so a refusal — a malformed read name a million records in, a truncated input — leaves nothing behind,
```

(Everything else in the paragraph is accurate.)

### H2 — a missed commit on the `--retain_unpaired` outputs is invisible to the suite (High, coverage)

Mutation M8, in a pinned worktree: replace `w.finish()?` with `drop(w)` at
`src/main.rs:1634-1639` (`unpaired_w1`, `unpaired_w2`) and the whole suite stays green. Nothing
asserts that `*_unpaired_1.fq` / `*_unpaired_2.fq` are produced.

The tripwire looks like it should catch it and structurally cannot: it asserts `created ⊆ planned`
only, and `tests/integration_preflight_tripwire.rs:1081` (`comparison_accepts_overplanning`) pins
that a *missing* output is legal — "planned may legitimately exceed created". `paired_trim_retain_unpaired`
(`:660`) therefore passes with both files absent.

The same mutation applied to the three other arms PLAN §7 declared "consciously uncovered" turns
8 tests red, so the reliance was justified for those:

| Site mutated | Result |
|---|---|
| `specialty.rs:332-333` clock | red — `clock_multipair`, `test_clock_happy_path_extracts_umi_and_clips` |
| `specialty.rs:427-428` implicon | red — `implicon_pair`, `test_implicon_happy_path_umi_8`, `test_implicon_custom_umi_length_6` |
| `main.rs:1631-1633` passthrough | red — `paired_trim_passthrough_fastqc`, `passthrough_cli_cores_1_smoke`, `passthrough_output_follows_r1s_directory` |
| `main.rs:1634-1639` unpaired | **green — no test detects it** |

This is PLAN §11 remaining risk (b) made concrete, and it also corrects the argument behind
deviation 2 in §12: what actually guards the ~20 mechanical sites is *tests*, not the compiler.
Deleting `flush()` forces each site to call *something*, but `drop(w)` compiles, and where no test
asserts the file exists, nothing at all protects it. Under this change the symptom of a missed commit
is worse than before, too — the file is absent rather than present-and-empty.

**Fix** — one assertion. `tests/integration_partial_output.rs` is the natural home; a run with
`--paired --retain_unpaired` over `test_files/` R1/short-R2 asserting both `_unpaired_N.fq` exist
closes it, and it belongs here because this PR is what made their absence the failure mode.

### M1 — the `released`-before-rename rule has zero coverage (Medium)

Mutation, in a pinned worktree — `src/io.rs:727-737`, `commit()` rewritten to set the flag after the
rename succeeds instead of before it:

```rust
std::fs::rename(&self.tmp, &self.final_path).with_context(|| { … })?;
self.released = true;
Ok(())
```

**701/701 pass.** PLAN §5 step 1 gives this rule a whole block, and rightly: with the flag set late,
a failed rename returns `Err`, `self` drops with `!released`, and `Drop` **deletes the temporary** —
the exact inverse of §4's promise that a failed rename leaves the data in place, and a data-loss
outcome rather than a cosmetic one. The rule is correct in the committed code; it is simply
unverified, so a future refactor can invert it silently.

The positive control rules out a blind harness: making `Drop` a no-op turned 12 tests red (all the
`refused_*_publishes_nothing` cases — `assert_nothing_published`'s second half is doing real work).

**Fix** — one test, and it is cheap. Put a *directory* at the final output path, run an otherwise
successful trim, and assert (a) non-zero exit, (b) the `.partial` still exists and is non-empty,
(c) stderr names both paths. `rename(file, existing_dir)` fails with `EISDIR`/`ENOTDIR` on both CI
platforms, so it works on `ubuntu-latest` and `macos-latest` alike, and it is the only test that
would exercise `commit()`'s error path at all.

### M2 — peak disk usage for an in-place re-run is now old + new, and nothing says so (Medium)

§6 prices the extra `rename(2)` and the slower error path, and states "no extra copy, no extra
memory". It does not price disk. Before, `File::create` truncated the previous output before the
first record was read, so a re-run's peak footprint was one copy of the output. Now the previous
output survives until the rename, so **a re-run over an existing result needs room for both**.

That is inherent to Option B and is the same mechanism as the recovery property §3 sells it on — but
for the Buckberry-scale files the performance story is built on, a re-run that used to fit can now
hit ENOSPC. And per H1 that is the failure mode which publishes a trailerless file with exit 0. One
sentence in the PR description, and arguably in the docs paragraph, so nobody files it as a
regression or discovers it at 500 GB.

### M3 — `refused_clump_only_publishes_nothing` never checks why the run failed (Medium, test)

`tests/integration_partial_output.rs:157-174` is the only case that does not go through `refuse()`,
and it asserts only `!output.status.success()`. `--clump_only --output-format ubam` is a combination
this project rejects at the CLI layer for a long and growing list of neighbours (`cli.rs::validate`
§3.4a), so if it ever joined that list the test would pass while testing nothing — an absence
assertion over a run that never opened a writer.

**Fix** — assert the reason, same as `refuse()` does:

```rust
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("read name contains whitespace"),
        "expected the whitespace refusal, got: {stderr}"
    );
```

### M4 — `refuse()` does not pin that the bail was mid-stream (Medium, test)

`refuse()` (`tests/integration_partial_output.rs:69-90`) asserts the refusal *reason* but not its
*position*, though its own doc comment claims "refused at record 2 (i.e. inside the trimming loop,
past writer creation)". The error carries the ordinal — `bam.rs:395`/`:452` wrap it as
`"{path}: BAM record {n}"` — and the test this helper was factored out of still asserts it
(`integration_ubam.rs:320`, `"record 2"`).

Without it, moving the whitespace check into `sanity_check_any` would make **every** absence
assertion in the file pass without a writer ever existing. That is the vacuity risk in its sharpest
form, and it is one `contains` away from being closed. Note `refuse()` already returns the stderr and
**no caller uses the return value** — so either the ordinal assertion goes inside `refuse` (taking it
as a parameter, since `OFFENDER` bails at record 2 and `OFFENDER_BULK` at 5001) or the `-> String`
should go.

---

## Logic — what I checked and found correct

**Commit-point invariant (lead's risk 1).** Verified by control flow, not line order, for all
consumers and all arms:

| Consumer | Commit | Ordering |
|---|---|---|
| `clump_only.rs:375` `fastqc::run` | `:349` `pending_output.commit()?` | ✅ nothing fallible between `out.flush()?`/`drop(out)` and the commit |
| `clump_only.rs:555-556` `fastqc::run` | `:521-522` both commits | ✅ both flushes precede both commits; the two commits are adjacent |
| `clump_only.rs:828` `metadata(&output_path)?` | `:823` `writer.finish()?` | ✅ the sharp one — the `?`-propagating stat on the success path |
| `clump_only.rs:1072` `metadata`, `:1096` fastqc | `:1069` `writer.finish()?` | ✅ |
| `main.rs:1459` fastqc, `:1474` `demux::demultiplex` | `:1334` `writer.finish()?` | ✅ report at `:1431` also after |
| `main.rs:2239/2396/2519` fastqc + report headers | `:2154/2318/2447` `writer.finish()?` | ✅ all three uBAM drivers |

`clump_only_single` additionally appends an empty gzip member when `output_bytes == 0` (`:341-343`)
*before* `out.flush()?` and the commit, so the legitimately-empty output is still published — T7's
overshoot guard holds structurally, not just by test.

**Drop ordering and `released` (lead's risk 2).** Both claims hold. `finish()` takes `pending`
before `drop(self)` in both writers, so the guard cannot remove a temporary that is about to be
published, and the trailer is on disk before the rename. `commit()` sets `released` before the
rename, so a failed rename keeps the data (untested — M1). Field declaration order is as the
comments say: `fastq.rs:578` `writer` before `pending`, `bam.rs:526-531` `inner` before `pending`.
On the error path `finish()`'s `?` fires before `take()`, so `self` drops with `pending: Some` and
the temporary goes — correct in both writers.

**The removed `flush()` (lead's risk 3).** Nothing needed a mid-stream flush. All 18 call sites now
end in `finish()`; `grep` cross-check of `FastqWriter::create` / `BamWriter::create` against
`.finish()` accounts for every non-test site, and there is no `return Ok` between any create and its
finish (checked in `specialty.rs`, `main.rs`, `demux.rs`, `clump_only.rs`). The only remaining
`flush()` calls are on `BufWriter`/`File` locals that mark a commit point, which is what step 4 asked
for.

**`parallel.rs`'s closure (lead's risk 4).** Borrow discipline is fine — `open_pending`'s `&mut`
borrow of `pending_outputs` ends at `:224` under NLL, before the `for guard in pending_outputs` at
`:312`. Error arms are all correct: a `?` on any of the five `open_pending` calls drops the guards
already pushed; the `Err(e)` arm at `:269-274` breaks to `:305` and returns before the commit; and
`reader_join_err.or(first_error)` returns before `:311`. Locals drop in reverse declaration order, so
the five `File` handles close before `pending_outputs` is dropped — correct, and it matters only on
Windows (see L5). `run_single_end_parallel` is the same shape: the reader-join `Err` arms at
`:986-987` both precede `drop(out); pending_output.commit()?` at `:990-991`.

**`--demux`'s fan-out (lead's risk 5).** Nothing fallible sits between the first and last commit —
`writers.drain()` + `sort_by` + N `finish()?` + the NoCode `finish()?`, then the summary. The stage
runs *after* the trimmed output was published, so a demux failure cannot un-publish it, which is
correct per §4.7. One real behaviour change hides here — see L6.

**Byte-identity (lead's risk 6).** No change. Every writer emits identical bytes through an unchanged
encoder; only the path they first land on differs, and `rename` is a directory operation. The three
new success-path tests pin the two trailer mechanisms and the BGZF EOF marker, `BamWriter::finish`
still calls `try_finish()` before dropping `inner` exactly as before, and `File::create` on the
temporary uses the same mode as it did on the final path. Two consecutive SE FASTQ runs md5 to
`bb256c8b53b36c6fc2c6f753a3c422be`. The `validation` job's five arms (SE, PE, hardtrim5, clock,
demux) all go through `FastqWriter`, which is unchanged below the file handle. The one thing that
does change is the published file's **mode** when a pre-existing output had a non-default one — md5
comparisons do not look at mode, so CI is unaffected (L4).

---

## Low

**L1 — `partial_output_name` is lossy, so the injectivity claim holds only for UTF-8 names.**
`io.rs:679-682` goes through `to_string_lossy().into_owned()`, which maps distinct invalid UTF-8 byte
sequences onto the same replacement character. Two final names differing only in invalid bytes would
therefore share one temporary. The tripwire header now rests that argument on injectivity, so it is
worth making unconditional:

```rust
    let mut tmp = std::ffi::OsString::from(".");
    tmp.push(final_path.file_name().unwrap_or_default());
    tmp.push(".partial");
```

Exact bytes, same allocation count, and `format!` disappears.

**L2 — `PendingOutput::tmp_path()` is dead code.** `io.rs:740` is never called; `commit()` formats
`self.tmp` inline. It is `pub`, so neither rustc nor clippy flags it. The plan's step-1 sketch listed
it; the implementation did not need it. Delete.

**L3 — the create error names the wrong path.** `io.rs:715-716` reports "Failed to create output:
{final}", but the failing syscall targeted `.{final}.partial`, nine bytes longer. On a name within
nine characters of `NAME_MAX` the user gets ENAMETOOLONG against a path that is demonstrably short
enough, with nothing to explain it. Name both paths, as `commit()` already does.

**L4 — two undocumented behaviour changes at the output path, both verified empirically.**
A **symlink** at the output path is now replaced by a regular file instead of being written through:
after a re-run, the output was `-rw-r--r--` and the symlink's former target still held the old
bytes. And a **write-protected (mode 0444) existing output is now silently replaced** — the re-run
exited 0 and the file came back as 0644. The old code could not do that: `O_TRUNC` on a 0444 file
returns EACCES (confirmed at the shell), so `File::create` refused the run. The same mechanism resets
any non-default mode on a re-run. All three are inherent to publish-by-rename and on balance
defensible; none is written down anywhere.

**L5 — guard-vs-handle drop order is load-bearing only on Windows.** It is guaranteed by declaration
order at four raw-`File` sites (`clump_only.rs:297`, `:440-441`, `parallel.rs:214`+`:220-224`,
`:970`) and commented at only the two writer sites (`fastq.rs:577`, `bam.rs:524-525`). Assumption 7
documents it for the writers. No action while the CI matrix is POSIX-only — unlinking an open file is
fine there — but the asymmetry is worth knowing if Windows ever appears.

**L6 — duplicate barcodes in a `--demux` samplesheet now lose a planned file outright.**
`read_barcode_file` (`demux.rs:31-102`) validates alphabet and length but not uniqueness, and
`writers` is keyed by barcode (`:195`), so a second entry with the same barcode drops the first
writer. That used to leave an empty file at the first sample's final name; now it leaves nothing.
Distinct sample names give distinct paths, so the pre-flight does not catch it. Pre-existing latent
bug whose symptom this change alters from "empty file" to "missing file" — below the filing bar on
its own, and arguably the more honest symptom, but the underlying missing validation is the real
defect if anyone wants it.

**L7 — comment length, three places.** House rule is one line, two maximum, fact not reasoning,
current state only:

- `demux.rs:252-255` — four lines, and the fourth is rationale ("sorting keeps which files exist
  deterministic if one rename fails"). That belongs in the commit message.
- `fastq.rs:653-657` — the third paragraph of `finish()`'s doc ("Teardown errors below this point
  are unreachable: … not that the bytes are complete") pre-empts an objection, which the rule sends
  to the commit message or the PR thread. The first two paragraphs earn their place: they state an
  ordering constraint a future editor would otherwise break.
- `io.rs:673-675` — the `publishDir` / glob rationale for the dot prefix is three lines of "why we
  chose this", which O2 already records in the plan.

The opposite case, which should **stay**: `io.rs:694-700`'s "**Commit where the writer is closed,
never later**" block. That is an invariant sitting on the thing it constrains, which is what the rule
asks for.

---

## Notes — real, no action, and explicitly not new issues

Per the filing bar: one issue is already planned for the swallowed teardown errors across
`flate2` / `gzp` / `noodles-bgzf` (PLAN §9). N1 and N2 are the same family — an error discarded at
teardown or in a worker, producing a truncated file under the final name — so they belong as extra
bullets on that one issue, not as new ones.

**N1 — SE parallel worker errors are still not observed.** `parallel.rs:947-950`: a worker prints
`"Worker error: …"` and breaks. The main loop then ends because the channel closed, `expected` may be
short of the batch count, there is no `first_error` on this path, and `pending_output.commit()?`
publishes a truncated output with **exit 0**. Not a regression — the old code wrote the same bytes to
the same final name — and practically unreachable, since everything fallible in
`process_single_batch` writes into a `Vec`. Worth noting for the asymmetry: the paired path fixed
exactly this by putting `Result<PairedBatchResult>` in the channel (the "B-Crit-1 fix" at `:88`), and
the single-end path never got the same treatment.

**N2 — a worker panic in the paired parallel path commits before the panic surfaces.** The main loop
exits via a closed channel with `first_error = None`, the reader joins `Ok(Ok(()))`, the commit runs,
and only then does `thread::scope` re-raise — a truncated output under the final name with exit 101.
Same pre-existing status as N1; the guards themselves are never leaked, because the closure's locals
drop normally before `scope` panics.

**N3 — SIGINT/SIGKILL leaves a hidden `.partial`.** Documented (O2) and now explained in the docs
paragraph. The final name is still never created, so the guarantee holds. No action.

**N4 — the four §12 deviations were the right calls.** Deviation 1 (dropping the
`Option<Box<dyn Write + Send>>`) is strictly better: it keeps the ordering guarantee via `take` +
`drop(self)` while removing an `expect`/panic path from `write_record`, which is on the hot loop.
Deviation 2 is the finding of the implementation and it changed the plan's safety argument, though
not as far as §12 claims — see H2. Deviation 3 is justified by §2.4's own reasoning (paired batching
only makes a paired residue *smaller*, so the SE 5000-record fixture is the stronger test of a large
plausible residue), and the paired arms still get absence cases. Deviation 4 is right: `FastqReader`
is the decoder the tool ships and needs no external test dependency.

**N5 — the prior plan reviews were folded in.** Spot-checked the four items I could verify against
the code: the commit-point invariant is implemented as an invariant (§5 step 4 → the doc block on
`PendingOutput`), `released`-before-rename is present, `--clumpify` and `--hardtrim3` both have
cases, and the tripwire header carries the injectivity argument plus the "a red `.partial` means a
missed commit" clause. The `--help` text at `cli.rs:291-296` and both `parallel.rs` comments were
updated, so step 9's four locations are all done.

---

## Recommendations, in the order I would take them

| # | Priority | Action |
|---|---|---|
| H1 | **High — before merge** | `docs/.../guide/outputs.md:124`: delete `, a full disk`. It names the one case §10 O5 defers, and that case exits 0. |
| H2 | High | One assertion that `--retain_unpaired` publishes `_unpaired_1/2.fq`. Proven uncovered by mutation; the tripwire cannot cover it. |
| M1 | Medium | One test for `commit()`'s error path — a directory at the final path. The `released`-before-rename rule is currently unverified. |
| M3 | Medium | Add the stderr assertion to `refused_clump_only_publishes_nothing`. |
| M4 | Medium | Assert the record ordinal in `refuse()`, or drop its unused `-> String`. |
| M2 | Medium | One sentence in the PR: a re-run now needs room for the old output and the new one at once. |
| L1–L3 | Low | `OsString` in `partial_output_name`; delete `tmp_path()`; name both paths in `create`'s error. |
| L4 | Low | One line somewhere that publish-by-rename replaces a symlink and resets a non-default mode. |
| L7 | Low | Trim the `demux.rs` four-liner, `finish()`'s third doc paragraph, and `io.rs`'s glob rationale. Keep the `PendingOutput` invariant block. |
| N1–N2 | — | Fold into the already-planned teardown-errors issue as extra bullets. Do not file separately. |
| L5, L6, N3 | — | Report only. |
