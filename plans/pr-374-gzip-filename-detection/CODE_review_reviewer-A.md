# Code review — PR #374 `fix(io): read gzipped input whose filename does not end in .gz`

**Reviewer A** · base `dev` (`7201f03`) · head `7a74458` · +225/−58 over `CHANGELOG.md`, `src/fastq.rs`, `src/format.rs`, `src/main.rs`

Worktree used for all builds and runs:
`/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/db842fb4-f582-4e29-8459-cfbe1cf77ceb/scratchpad/wtA`
Pre-fix reference binary: `/Users/fkrueger/Github/TrimGalore/target/release/trim_galore` (`c5c1834`, two docs/CI commits behind `dev` — no source-path difference).

**Fixes applied: none.** This is an external contributor's branch; everything below is a recommendation. No file in either checkout was modified — `git status --short` in the worktree is empty. My experiments wrote only to a scratch directory (`scratchpad/bgztest`) and to `wtA/target/`.

---

## Summary

The diagnosis is correct, the mechanism of the fix is right, and it verifiably works on the paths it touches. I reproduced the bug and confirmed the fix end-to-end, and confirmed the trimmed output is byte-identical to a correctly-named `.fq.gz` run.

But **the fix is incomplete, and the incompleteness is not disclosed.** `--paired` at the default `--cores 1` still reads the input through the filename heuristic (`src/main.rs:1348-1349`), as does `--clump_only` with FASTQ output (`src/clump_only.rs:290,430,432`) and `--passthrough` (`src/main.rs:743,1319,1358`). I confirmed all three still fail on a `.bgz` input on the PR head.

For `--paired` this is worse than a missed case. Pre-fix, the run died in `sanity_check` before anything was created. Post-fix, `sanity_check_any` correctly passes, adapter auto-detection correctly scans the compressed file, the output directory is created, **two zero-byte `_val_{1,2}.fq` files are written**, and only then does the same `stream did not contain valid UTF-8` appear. The PR's own CHANGELOG entry cites "before any output directory was created" as a property of the old failure; on the most common invocation in the tool it has traded that property away without fixing the bug.

The six new tests all assert real behaviour and all pass, but none of them covers a paired run or the binary at all — which is exactly why the gap survived nine green CI checks.

`cargo fmt --all -- --check` clean · `cargo clippy --all-targets --release -- -D warnings` clean · `cargo test` 375 unit + 87 integration tests, 0 failures.

**Merge recommendation: request changes.** Not because what is here is wrong, but because the headline claim overstates what ships, and the paired path regresses in failure shape. The remaining call sites can be fixed in roughly five lines (§Critical-1), which is less work than documenting the exception.

---

## Verification performed

| Case | Pre-fix (`c5c1834`) | PR head (`7a74458`) |
|---|---|---|
| SE, plain-gzip content named `.fq.bgz` | `Error: stream did not contain valid UTF-8` | **works**, 10000 reads |
| SE, genuine BGZF-framed FASTQ named `.fq.bgz` | (n/a) | **works**, 10000 reads |
| **PE `--paired`, default `--cores 1`, `.bgz`** | error, empty output dir | **still errors**, + two 0-byte `_val_*.fq` |
| PE `--paired --cores 2`, `.bgz` | error | **works** |
| `--clump_only`, `.bgz` | error | **still errors**, + 0-byte `_clumped.fq` |
| `--hardtrim5`, `.bgz` | error | **works** |
| `--passthrough x.bgz` (R1/R2 proper `.gz`, `--cores 2`) | error | **still errors** |
| SE, plain FASTQ misnamed `.fq.gz` | `Error: invalid gzip header` | **now accepted**, output re-gzipped |
| PE `--cores 1`, plain FASTQ misnamed `.fq.gz` | `invalid gzip header` | **still** `invalid gzip header` |
| Plain FASTQ on a pipe (`gzip -dc … \| tg /dev/stdin`) | `doesn't seem to be in FastQ format` | **identical** |
| Gzip stream on a pipe | `invalid gzip header` | **identical** |

Byte-identity of trimmed reads, md5 of decompressed `_trimmed.fq`:

```
12fde5c81a6780ecfa57e54b7953894a   pre-fix,  g_R1.fq.gz  (correctly named)
12fde5c81a6780ecfa57e54b7953894a   PR head,  g_R1.fq.gz  (correctly named)
12fde5c81a6780ecfa57e54b7953894a   PR head,  R1.fq.bgz   (same bytes, .bgz name)
039657a84d6093a051627b734670358c   sentinel (same file + one appended line)
```

The sentinel line proves the comparison can fail; `diff -q` reported the sentinel as differing. Process substitution is blocked in this sandbox, so every comparison went through real temp files.

**"Byte-identical duplicates" claim — verified true.** `git show dev:src/format.rs` (lines 252-285) against the `main.rs` bodies removed by this diff: the two pairs differ only in `crate::`-qualification of the `fastq`/`bam` paths and in the doc comment. Both took `&Path`, both dispatched on the same three-variant match, both threaded `preserve_tags` to the BAM arm only. Nothing changed silently in the deletion.

**uBAM discrimination — not weakened.** `detect_input_format` is untouched by this diff; the `BAM\1`-in-decompressed-payload discriminator is unchanged. `detect_bgzipped_fastq_is_fastq_not_bam`, both `detect_*_ubam_via_decompressed_magic` tests, and all 8 + 23 uBAM integration tests pass. My hand-built genuine-BGZF FASTQ classified as FASTQ and trimmed correctly.

**Thread-boundary crossing — sound.** `open_threaded_with` moves a `bool` (`Copy + Send`) into the `spawn` closure and applies it inside via `open_direct` → `open_reader` (`src/fastq.rs:314,382,367`). No shared state, no ordering concern. Confirmed by `format::tests::threaded_reader_decompresses_gzip_under_non_gz_extension` and by a real `--paired --cores 2` run on `.bgz`.

---

## Issues by area

### Logic

**L1 (Critical). The fix covers three of the eight production call sites that take a user-supplied input path.** The PR fixes `format::open_sync_reader`, `format::open_threaded_reader` and `main::sanity_check_any`. Still filename-derived:

| Site | Path | Reached by |
|---|---|---|
| `src/main.rs:1348-1349` | `FastqReader::open(input_r1/r2)` | **`--paired` at `--cores 1` — the default** |
| `src/main.rs:1358` | `FastqReader::open(p)` | `--passthrough`, serial |
| `src/main.rs:1319` | `FastqReader::open_threaded(p)` | `--passthrough`, parallel |
| `src/main.rs:743` | `FastqReader::sanity_check(pt_path)` | `--passthrough`, pre-flight |
| `src/clump_only.rs:290` | `FastqReader::open(input)` | `--clump_only` SE, FASTQ out |
| `src/clump_only.rs:430,432` | `FastqReader::open(input_r1/r2)` | `--clump_only --paired`, FASTQ out |

(`src/demux.rs:170` reads back Trim Galore's own trimmed output, where the filename genuinely is authoritative — correctly left alone. `src/specialty.rs:531` and the `parallel.rs` hits are test helpers.)

Observed, PR head, default flags:

```
$ trim_galore --paired -o out R1.fq.bgz R2.fq.bgz
Adapter: Illumina (AGATCGGAAGAGC)                     ← detection read the .bgz fine
  Output R1: …/out/R1.fq_val_1.fq
Error: processing pair 1 of 1 (…)
Caused by: stream did not contain valid UTF-8
$ ls -la out
-rw-r--r--  0  R1.fq_val_1.fq
-rw-r--r--  0  R2.fq_val_2.fq
```

Pre-fix the same command left `out/` empty. The mirror case is the same defect: a plain FASTQ misnamed `.fq.gz` now works in SE and still fails in PE with `invalid gzip header`. So the observable behaviour of `.bgz` (and of a misnamed `.gz`) now depends on `--cores`, which no user would predict.

**L2 (High). The CHANGELOG entry claims more than the diff delivers.** It says "The detected format is now passed down to all three" — true of the three *functions* named, but a reader will take away "gzipped input under a non-`.gz` name now works". It does not, for `--paired` (default), `--clump_only`, or `--passthrough`. Either fix those or scope the entry to the modes that work.

**L3 (Medium). `io::is_gzipped`'s doc comment is now false.** `src/io.rs:12-14`: "The same heuristic that `FastqReader` uses to decide whether to wrap the input in a gzip decoder, so output naming + reader behaviour stay consistent." After this PR the reader does *not* use that heuristic on the primary paths, and the two can now disagree — which is precisely the "known limitation" the CHANGELOG records. The comment documents an invariant the PR breaks; it should say so.

**L4 (Medium). Output naming on a `.bgz` input is newly reachable and ugly.** `io::strip_fastq_extensions` (`src/io.rs:425-451`) matches only `.fastq.gz`/`.fq.gz`/`.fastq`/`.fq`, then falls back to `file_stem()`. Observed on the PR head:

```
input  sample.fastq.bgz
output sample.fastq_trimmed.fq                   ← stem keeps ".fastq"
report sample.fastq.bgz_trimming_report.txt      ← report keeps the full name
```

MultiQC derives the sample name from the report filename, so this run reports as `sample.fastq.bgz` against a file called `sample.fastq_trimmed.fq`. Before this PR the run could not get far enough to produce either. Not a defect the PR introduced, but a user-visible consequence it unlocks; `.bgz`/`.bgzf` belong in the strip list.

**L5 (Medium). Deferring the output-compression decision is defensible; leaving it silent is not.** `src/main.rs:285` — `let gzip = !cli.dont_gzip && naming::is_gzipped(&cli.input[0]);`. I agree with the contributor that moving this to the detected format is a separate change: it would flip behaviour for a plain file misnamed `.gz` (which now reads as plain but still writes gzip — verified: `plain_R1.fq.gz` → `plain_R1_trimmed.fq.gz`), and that touches #245 parity. But a `.bgz` user silently gets uncompressed output and no `--gzip` escape hatch (`--gzip` is deprecated and ignored, `src/cli.rs:401-402,966`). One `NOTE:` line on stderr when the detected format is `FastqGz` and `is_gzipped(input[0])` is false would close the gap without touching the parity decision.

**L6 (Low). Accepting a plain file misnamed `.gz` is a real behaviour change, undocumented as one.** Verified: pre-fix `Error: invalid gzip header`, post-fix a successful SE run. This is a loosening (Perl 0.6.11 would also have failed), almost certainly welcome, but it means a truncated or accidentally-decompressed `.gz` no longer announces itself. The CHANGELOG mentions the case only as a reason to defer the *output* change, not as something this PR already changes on the input side.

### Efficiency

**E1. Neutral — no I/O added or removed.** `open_with` calls `File::open` exactly once, as the old `open` did; `open_direct` opens once via the new `open_reader`, as before. `detect_input_format` is untouched, so its two opens on the gzip path (`src/format.rs:208` and `:231`) are unchanged in count and position. There is no measurable cost to this diff.

**E2 (Low). `open_with` duplicates the `String::with_capacity(512)` that `open_direct` already produces.** `src/fastq.rs:275-281` vs `:382-385`. `open_with` could be `let (reader, line_buf) = Self::open_direct(path, is_gzip)?;`, which also removes the second literal `512`.

### Errors

**Er1 (Critical, same root cause as L1). The failure mode on `--paired` and `--clump_only` moved from pre-flight to mid-run.** `format::reject_bam_format_mismatch_in_pair`'s own doc (`src/format.rs:95-99`) states the project's position: the guard exists so failures land "before any output directory is created — which is the point", because the sites it replaced "fired only after adapter auto-detection had scanned the inputs, the output directory existed, and (on multi-pair input) earlier pairs had already been written to disk." That is exactly the shape `.bgz` on `--paired` now has. Zero-byte `_val_1.fq`/`_val_2.fq` left in an output directory is the kind of artefact a pipeline glob will happily pick up.

**Er2. Pre-existing non-seekable-input bug: the PR is neutral.** I reproduced it on both binaries and the failure is character-for-character identical:

```
gzip -dc R1.fq.bgz | trim_galore -o out /dev/stdin
  pre-fix : Error: Input file '/dev/stdin' doesn't seem to be in FastQ format (first line doesn't start with '@')
  PR head : Error: Input file '/dev/stdin' doesn't seem to be in FastQ format (first line doesn't start with '@')

cat R1.fq.bgz | trim_galore -o out /dev/stdin
  pre-fix : Caused by: invalid gzip header
  PR head : Caused by: invalid gzip header
```

Statically: `detect_input_format` is not in the diff, and `open_with`/`open_direct` each perform the same single `File::open` the functions they replace did. So the open count per run is unchanged and no new byte is consumed from a pipe.

One design observation. The PR threads the *verdict* across the detection boundary but not the *consumed prefix*, so it fixes the classification half of the double-open problem and leaves the stream-position half untouched. A `bool` cannot carry four already-read bytes; the eventual fix (the `PeekReader` wrapper `src/format.rs:229-230` already points at) will want to hand a `Read` downstream rather than a flag. Worth knowing that `open_with(path, bool)` is not the shape that fix will need — it is not on the way to it, nor in its way.

**Er3 (Low, pre-existing, untouched). `detect_input_format` requires `nread == 4`** from a single `MultiGzDecoder::read` (`src/format.rs:234,240`). `Read::read` is permitted to return short; `read_exact` would be correct. Unreachable in practice for a real BAM (first BGZF block is ~64 KB), and outside this diff — noted only because the recommended follow-up work lands next to it.

### Structure

**S1 (High). `is_gz_filename` is a second copy of `io::is_gzipped`.** `src/fastq.rs:37-39` and `src/io.rs:19-21` are the same expression, and `io.rs`'s own doc says they are meant to be one concept. The PR gave the heuristic a name — the right instinct — but created a duplicate instead of reusing the existing name. `fastq.rs` should call `crate::io::is_gzipped`, or the predicate should move to `fastq.rs` and `io.rs` re-export it. As it stands a future change to one silently diverges from the other.

**S2 (High). Six public entry points where three would do, and the three weak ones are still the default choice.** `open`/`open_with`, `open_threaded`/`open_threaded_with`, `sanity_check`/`sanity_check_with`. The un-suffixed member of each pair is the footgun, has the shorter name, and is what the remaining unfixed call sites use. Per `CLAUDE.md`, the library crate exists "so unit tests can reach internals" — there is no external API-compatibility reason to keep the weak trio `pub`. See Critical-1 for the version of this that also fixes L1.

**S3 (Medium). `open_with(path, bool)` is boolean-blind, and `fmt == InputFormat::FastqGz` is written out three times.** `src/format.rs:265,283` and `src/main.rs:36`. A transposed or inverted argument compiles silently. Two cheap improvements: take `InputFormat` (or a two-variant `Compression`) instead of `bool`; and put the predicate on the enum as `InputFormat::is_gzip(self) -> bool` implemented with a `match`, so adding a fourth variant (the pluggable-I/O epic has zstd/BINSEQ in view) forces every site to be revisited instead of silently defaulting to "plain".

**S4 (Medium). Enshrining `sanity_check(&p).is_err()` in a test is the wrong call.** `src/fastq.rs:1001-1011`. The assertion is real — it fails today for the documented reason — but it pins the *broken* behaviour of a `pub` function. Someone implementing Critical-1 will see this test go red and reasonably read it as a regression they caused. A test should not make the correct fix look like a break. Delete the assertion (the doc comment on `open` already carries the warning), or, better, delete the function it protects.

**S5 (Low). Doc comments narrate bug history in the past tense.** `src/fastq.rs:517-520` ("on a `.fq.bgz` input it read the compressed bytes as text and failed with…"), `src/format.rs:299-302` ("Before this change the run died with…"), `src/fastq.rs:961-967`. `CLAUDE.md`'s comment convention asks for one line stating the fact, with the archaeology in the commit message. I'll note the surrounding code is itself heavily commented in this register (`src/format.rs:113-122`, `src/main.rs:1293-1301`), so this is a maintainer's call rather than a violation — but the past-tense framing specifically will read as stale once the bug is old, in a way "the verdict comes from the caller, not the filename" would not.

**S6 (Low). Stale struct doc.** `src/fastq.rs:246-252` still describes only "**Direct** (`open`)" and "**Threaded** (`open_threaded`)"; the `_with` variants are now the recommended entry points and are unmentioned.

### Test quality

All six new tests assert real behaviour; none passes trivially. Specifically:

- `format::tests::sync_reader_decompresses_gzip_under_non_gz_extension` and its threaded twin are the genuine regression tests — they exercise the factories that base `dev` got wrong, and they would fail on base `dev`.
- `fastq::tests::open_with_*` / `sanity_check_with_*` test an API that does not exist on base `dev`, so the "REGRESSION" label on `open_with_reads_gzip_under_non_gz_extension` (`src/fastq.rs:961`) overclaims slightly — it is a unit test of a new parameter, which is fine, just not a regression test.
- `open_with_honours_a_plain_verdict_on_a_gz_name` is the right mirror-image test and would catch an inverted flag.

**T1 (High). Nothing covers a paired run, and nothing drives the binary.** `git grep -n bgz -- tests/` returns nothing. All coverage sits at the `FastqReader`/factory level, and both unfixed hot spots (`main.rs:1348-1349`, `clump_only.rs:290`) bypass the factories entirely — which is why nine green checks did not notice. A `tests/integration_gz_extension.rs` that runs the built binary on a `.bgz` fixture across SE, `--paired --cores 1`, `--paired --cores 2` and `--clump_only` would have caught L1 immediately, and is the test I would gate the merge on.

**T2 (Low). `gz_tmpdir` is a third copy of the same helper.** `src/fastq.rs:945-951` duplicates `format.rs`'s `fresh_tmpdir` (and `clump_only.rs` has its own). Existing convention, not introduced here; a shared `#[cfg(test)]` helper would be the tidy-up.

---

## Recommendations

### Critical

**C1 — Cover the remaining call sites, ideally by making `open` self-sufficient.** The five-line version, in `src/fastq.rs`, replacing `is_gz_filename`:

```rust
/// True iff the first bytes of `path` are the gzip magic. Content-based, so
/// `.fq.bgz` and a plain file misnamed `.gz` are both classified correctly.
/// Unreadable or short files fall back to "not gzip" — the caller's own open
/// will produce the real error.
fn probe_is_gzip(path: &Path) -> bool {
    let mut buf = [0u8; 2];
    File::open(path)
        .and_then(|mut f| f.read(&mut buf))
        .map(|n| n == 2 && buf == [0x1F, 0x8B])
        .unwrap_or(false)
}
```

`open`, `open_threaded` and `sanity_check` then call this instead of the filename, and `main.rs:1348-1349`, `main.rs:1319/1358/743` and `clump_only.rs:290/430/432` are all fixed with no edits. Three properties make this safe where delegating to `format::detect_input_format` would not be:

- it cannot misfire on a plain FASTQ (first byte `@`, never `0x1F`);
- it never introduces a new error — `detect_input_format` *bails* on empty and unrecognised files, which would newly break `demux.rs:170` reading back a legitimately empty all-filtered output, and the `parallel.rs` read-back tests;
- BAM-vs-FASTQ discrimination stays in `format.rs` where it belongs; this answers only "do I need a gzip decoder".

Cost is one open + 2-byte read per reader construction, against reading the whole file. It also removes the need for the `_with` overloads on the sync path, though keeping them for callers that already hold a verdict is reasonable.

If the contributor prefers to keep the change minimal, the acceptable minimum is: pass the detected verdict into `run_paired` and `clump_only` (both already have the format in hand — `main.rs:1302-1311` calls `detect_input_format` on both pair members for the BAM backstop, so the verdict is *literally already computed and discarded* two dozen lines above the broken `FastqReader::open`), plus §H1 below.

### High

**H1 — Do not leave zero-byte outputs behind.** Whichever route C1 takes, the paired and clump-only paths should reach their decision before `FastqWriter::create`. `main.rs:1302-1311` is already the natural place: it loops over both pair members calling `detect_input_format` for the uBAM backstop, immediately before the readers are opened.

**H2 — Scope the CHANGELOG to what ships.** If C1 lands, the current wording becomes accurate. If it does not, name the modes that still require a `.gz` name.

**H3 — Add binary-level coverage.** `tests/integration_gz_extension.rs` per T1. This is the single highest-value addition to the PR.

**H4 — Reuse `io::is_gzipped` rather than adding `is_gz_filename`** (S1), and fix the now-false invariant comment at `src/io.rs:12-14` (L3).

**H5 — Make the weak trio non-public** (S2), or delete it if C1 lands.

### Medium

- **M1** Take `InputFormat` rather than `bool`, and move the gzip predicate onto the enum as a `match` (S3).
- **M2** Drop the `sanity_check(&p).is_err()` assertion (S4) — it makes the correct fix look like a regression.
- **M3** Emit a one-line stderr `NOTE:` when the detected format is gzip but the output will be plain because the filename says so (L5). Keeps the #245 parity decision untouched while removing the silence.
- **M4** Add `.bgz` (and `.bgzf`) to `io::strip_fastq_extensions`, so a `.bgz` run does not produce `sample.fastq_trimmed.fq` and a MultiQC sample called `sample.fastq.bgz` (L4).

### Low

- **L-a** `open_with` should reuse `open_direct` and drop the duplicated `512` (E2).
- **L-b** Trim the past-tense bug narration in the new doc comments (S5); update the `FastqReader` struct doc to mention the `_with` variants (S6).
- **L-c** Note in the CHANGELOG that a plain file misnamed `.gz` is now accepted rather than rejected (L6).
- **L-d** Out of scope but adjacent: `detect_input_format` should use `read_exact` for the 4-byte payload probe (Er3), and the non-seekable/pipe double-open (Er2) deserves its own issue — the PR neither helps nor hurts it, and the `bool` handoff is not the shape that fix will need.
