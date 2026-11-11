# CODE REVIEW B — #428 `0a6f3a3` (`fix/428-no-partial-output-on-refusal`)

Reviewer B, fresh context. Parent `285afba`. Plan: `PLAN.md` v2 (§12 implementation notes read).
No file in the repository was modified; every fix below is a proposed replacement, not an edit.

---

## Verdict

**Safe to merge**, with one documentation line I would fix first (F1) and one measured regression
to decide on (F2). The design is right, the commit-point invariant holds under control-flow
analysis (not just line order), and nothing here changes output bytes on any success path, so the
`validation` job is not at risk.

What the change does *not* have is test coverage of the two ordering rules it rests on. I proved
that by mutation: both survive a full green suite. That is a coverage finding, not a defect — the
code is correct today — but it means a future edit can silently undo either rule.

## Verification performed

| Check | Result |
|---|---|
| `cargo test --release` (main tree, read-only) | **701 passed, 0 failed** — matches §12 |
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --all-targets --release -- -D warnings` | clean |
| Mutation A — `released = true` moved *after* the rename (fails → `Drop` deletes the data) | **701 green: survives** |
| Mutation C — `finish()` commits *before* dropping the trailer-writing inner writer | **701 green: survives** |
| Mutation D — `run_paired`'s two `--retain_unpaired` `finish()` calls → `drop(w)` | **701 green: survives** |
| `--clump_only --output-format ubam` on the offender, by hand | refused at `BAM record 2`, nothing left in the output dir |
| ENAMETOOLONG band, by hand | measured, see F2 |

Mutations ran in a worktree pinned to `0a6f3a3` with a private `CARGO_TARGET_DIR`; it was reverted
and `git worktree remove`d. `git worktree list` and `git status` confirm the main tree is untouched.

## What I checked and found sound

- **The commit-point invariant, by control flow.** Every commit precedes every in-run consumer on
  every path, not merely by line number: `main.rs:1334` → `fastqc` `:1459` / `demux` `:1474`;
  parallel SE commits inside `run_single_end_parallel` (`parallel.rs:991`) before it returns to
  `:1329`; `clump_only.rs:349` → `write_clump_only_report` `:368` → `fastqc` `:375`; `:521-522` →
  `:550-551` → `:555-556`; `writer.finish()` `:823` → `metadata(&output_path)?` **`:828`** → report
  `:841` → `fastqc` `:848`; the same shape at `:1069` → `:1072` → `:1091` → `:1096`. The sharp one
  (the `?`-propagating stat) is genuinely covered.
- **No missed record-output site.** `grep File::create src/` leaves only reports
  (`main.rs:1431/1450/2213/2230/2596/2621`, `demux.rs:286`, `clump_only.rs:591`) and `mod tests`
  bodies. Every `FastqWriter::create` / `BamWriter::create` in the lib and binary has a matching
  `finish()` on its Ok path; I checked each of the 30 sites individually, and there is no early
  `return Ok` between any create and its finish.
- **Drop ordering and the `released` flag.** `FastqWriter` has a pre-existing `Drop` (`fastq.rs:670`)
  that flushes, so `finish()` cannot destructure; `pending.take()` + `drop(self)` is the only shape
  that works, and field order (`writer` then `pending`) closes the compressor before the temporary is
  removed. `BamWriter` has no `Drop`, so `drop(self)` is exactly field-order drops — `inner` then
  `pending`. Both comments describe what the code does.
- **`parallel.rs`'s `&mut` closure.** `pending_outputs` lives on the main thread only; no guard can
  reach a worker. On any `?` inside the scope closure the guards drop, the temporaries go, and the
  implicit `thread::scope` join still completes (the `result_rx` drop that unblocks it is inside the
  closure). Commit is after the reader join **and** after `reader_join_err.or(first_error)`, so no
  error arm can publish.
- **The removed `flush()`.** Every deleted call sat at teardown (immediately before a `drop` or the
  end of scope); no site wanted a mid-stream flush, and nothing in `trimmer.rs` could have called one
  or the build would fail. Deviation #2's reasoning is correct and was the right call.
- **Byte-identity.** No write path changed; the only new syscall is `rename(2)`. The flush sequence
  is unchanged too (explicit flush, then `Drop`'s flush, then the inner `try_finish` — same as
  before), so gzip and BGZF bytes are identical. `@PG CL:`, the reports and FastQC all still see the
  final path. Nothing asserts on the two error strings that were unified into
  `"Failed to create output"` (grep over `tests/`, `docs/`, `.github/` is empty). **SE, PE,
  hardtrim5, clock and demux are unaffected.**
- **Tripwire interaction.** `snapshot()` uses `read_dir`, so it sees dotfiles; a leaked `.partial`
  would land in `created` and fail `created ⊆ planned` on all 11 arms. That is a real extra net, and
  the new module-doc note describes it accurately. The injectivity argument holds, including under
  `collision_key`'s folding, because the marker is a pure prefix/suffix of the file name.
- **The four §12 deviations were each the right call.** #1 (no `Option<Box<dyn Write>>`) avoids an
  `expect` on the hot `write_record` path and gets the same ordering from field-drop order. #2 as
  above. #3 — 5000 SE records is a better large-plausible-residue fixture than a paired one, since
  `BAM_READER_BATCH_SIZE` batching only makes a paired residue *smaller*. #4 — reading back through
  `FastqReader` keeps the no-external-deps rule and uses the same decoder the tool uses.

---

## Findings

### F1 — High: `outputs.md` tells users to delete a `.partial` that may be the only copy of their data

`docs/src/content/docs/guide/outputs.md` (new section):

> A leftover `.<name>.partial` therefore means the process was killed rather than that it finished.
> It is safe to delete; the next run at the same output path overwrites it.

There is a second way to get a leftover, and it is the one the design deliberately created:
`PendingOutput::commit` sets `released = true` *before* the rename precisely so a **failed rename
leaves the temporary on disk**, and the error says `"the trimmed data is intact at <tmp>"`
(`io.rs:729-739`). The docs and the error message now contradict each other, and the docs are the
side that loses data. The trigger is rare (a directory sitting at the output path → `EISDIR`, an
immutable file, the output directory vanishing mid-run), but the whole point of the `released`
ordering is that this case preserves the data.

**Proposed replacement for the second paragraph** (one sentence, current state):

```markdown
A leftover `.<name>.partial` means the process was killed mid-run, or that the final rename failed — in that case the error names it as holding the trimmed data. Otherwise it is safe to delete; the next run at the same output path overwrites it.
```

### F2 — Medium: a legal output name in a 9-character band now fails with `File name too long`

`io.rs:713-717` creates `.{file_name}.partial`, nine bytes longer than the final name. Where the
final name is 247–255 bytes, the temporary exceeds the 255-byte `NAME_MAX` and the run dies at
output-creation time. Measured on this machine (APFS):

```
$ trim_galore --hardtrim5 10 -o D  <235-n stem>.fastq.gz     # output name 253 bytes
Error: Failed to create output: …/nnnn….10bp_5prime.fq.gz
Caused by:
    File name too long (os error 63)
```

and the same for a plain SE trim with `--no_report_file`. I confirmed by direct syscall that a
254-byte final name **is** creatable in that directory and a 263-byte temporary is not, so at
`285afba` these runs succeeded.

The band is only reachable where no report is written — with reports on, `{input}_trimming_report.txt`
is already too long at those lengths, so the old code failed anyway, just later. So: `--hardtrim5`,
`--hardtrim3`, `--clock`, `--implicon`, and any run with `--no_report_file`.

It fails loudly, before any work, with no data loss, and ~235-character sample names are unusual —
so my recommendation is **accept it and put one clause in the PR description**, not to complicate the
naming. No derived-name scheme can avoid it (any prefix or suffix adds at least one byte), and the
only fully-safe alternative — shortening the base name — is what the pre-flight's injectivity
argument depends on not happening. If you would rather not regress at all, the minimal fix is a
fallback triggered by the error rather than by a length guess:

```rust
// in PendingOutput::create, replacing the single File::create
let file = match std::fs::File::create(&tmp) {
    Ok(f) => f,
    // A nine-byte-longer sibling can exceed NAME_MAX where the final name does not.
    Err(e) if e.raw_os_error() == Some(63) || e.kind() == std::io::ErrorKind::InvalidFilename => {
        // …shorten and retry; note this breaks the injectivity claim in
        // tests/integration_preflight_tripwire.rs, which would need rewording.
        return Err(e).with_context(|| format!("Failed to create output: {}", final_path.display()));
    }
    Err(e) => return Err(e).with_context(|| format!("Failed to create output: {}", final_path.display())),
};
```

I am **not** recommending that shape — I include it to show the cost. My advice is the PR clause.

### F3 — Medium: both ordering rules the design rests on are invisible to the suite

Proven, not inferred:

- **Mutation A.** Rewrite `commit` so `released = true` runs only after a successful rename (i.e. the
  exact defect PLAN §5 step 1 sets in bold, where a failed rename deletes the data `Drop` was
  protecting). **701 tests still pass.** No test ever provokes a rename failure.
- **Mutation C.** Rewrite `finish` to commit *before* `drop(self)`. **701 tests still pass** — and
  this one is worth understanding rather than just patching: on POSIX the rename does not invalidate
  the open descriptor, so the trailer still lands in the same inode, now under the final name. T2
  cannot catch this reordering, and no test can, unless the process dies inside the window.

Consequence for the source comment at `fastq.rs:648-651`: *"renaming before that runs would publish a
truncated file"* is only true if the process is killed between the rename and the drop. The current
order is the right one — it is the one that survives a SIGKILL in that window — but the comment
states a stronger mechanism than exists.

**Recommended test** — closes the A gap directly, portable to both CI platforms (`rename(2)` returns
`EISDIR`/`ENOTEMPTY` when the destination is a non-empty directory). Add inside `src/io.rs`'s
`mod tests`:

```rust
#[test]
fn a_failed_commit_leaves_the_temporary_on_disk() {
    let dir = std::env::temp_dir().join("tg_pending_commit_fail");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let final_path = dir.join("out.fq");
    let (pending, mut file) = PendingOutput::create(&final_path).unwrap();
    let tmp = pending.tmp_path().to_path_buf();
    use std::io::Write;
    file.write_all(b"@r\nA\n+\nI\n").unwrap();
    drop(file);

    // A non-empty directory at the final path makes the rename fail.
    std::fs::create_dir(&final_path).unwrap();
    std::fs::write(final_path.join("blocker"), b"x").unwrap();

    assert!(pending.commit().is_err(), "renaming onto a directory must fail");
    assert!(tmp.exists(), "a failed commit must leave the trimmed data at the temporary");
}
```

While there: `partial_output_name` has no unit test at all, and its `file_name() == None` branch
silently yields `..partial`. Two asserts (`.../x.fq` → `.../.x.fq.partial`, and a bare `x.fq` →
`.x.fq.partial`) would pin the naming that the tripwire's injectivity note depends on.

### F4 — Medium: a missed commit on `--retain_unpaired` would ship silently

Mutation D replaced `main.rs:1635`/`:1638`'s `w.finish()?` with `drop(w)` — so `--retain_unpaired`
publishes no unpaired files at all — and **701 tests stayed green**. `--retain_unpaired` appears in
`tests/` only in `integration_preflight_tripwire.rs`, and that tripwire asserts `created ⊆ planned`
only; over-planning is legitimate by design (`:436-438`), so it structurally *cannot* see a missing
output. `main.rs:1898`/`:1901` (the interleaved-uBAM twin) is in the same position.

This is the plan's own §11 risk (b) — "a missed one is a silent hole" — and `finish()` consuming
`self` does not mitigate it, exactly as deviation #2 discovered for `flush()`. The passthrough writer
is fine: `tests/integration_passthrough.rs:133-157` asserts existence *and* record count for all
three outputs.

**Recommended**: one assertion on the serial `--retain_unpaired` arm — the fixture already exists
(`write_fastq_short` in the tripwire strands reads deliberately). Cheapest form, in
`tests/integration_passthrough.rs`'s style: run `--paired --retain_unpaired` on a pair whose R2 is
short, then

```rust
assert!(dir.join("s_R1_unpaired_1.fq").exists(), "unpaired R1 must be published");
assert_eq!(count_records(&dir.join("s_R1_unpaired_1.fq")), 40);
```

This is also the T7 direction (a legitimately-produced file must still appear), which is open as
task #32.

### F5 — Medium: `refused_clump_only_publishes_nothing` does not pin why the run failed

`tests/integration_partial_output.rs:157-174` is the one case that bypasses `refuse()` and asserts
only `!status.success()`. An `assert_nothing_published` after an unrelated early failure — a CLI
validation change, a future refusal of `--clump_only` with `--preserve-tags` — passes vacuously, and
that arm's writer is `clump_only.rs:1006`, one of only two `BamWriter` sites outside `main.rs`.

I verified by hand that it currently fails for the right reason (`Error: …: BAM record 2` →
`BAM read name contains whitespace`), so this is hardening, not a bug.

**Proposed addition after line 172:**

```rust
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("BAM record 2") && stderr.contains("read name contains whitespace"),
        "expected the mid-loop whitespace refusal, got: {stderr}"
    );
```

Related, and cheaper: `refuse()` at `:69-90` asserts the message but not the *position*. Its own doc
comment claims "refused at record 2 (i.e. inside the trimming loop, past writer creation)", which is
the property that keeps all eleven absence assertions non-vacuous — if that check ever moved to
`sanity_check_any`, every one of them would pass with no writer ever created.
`integration_ubam.rs:319-323` already pins `BAM record 2`; the new file should too:

```rust
    assert!(
        stderr.contains("BAM record 2"),
        "must be refused inside the trimming loop, past writer creation, got: {stderr}"
    );
```

### F6 — Low: a mid-sequence commit failure deletes the outputs that had not committed yet

`parallel.rs:312-314`, `main.rs:1629-1639`, `demux.rs:258-261`: on `guard.commit()?` returning `Err`,
the remaining guards drop and their temporaries are **removed**. So the failing output's data is
preserved (F1) but every later output in the fixed order loses its `.partial` too. PLAN §4's
edge-case row also asked for a "hard error naming what committed and what did not"; the message names
only the failing pair.

Rename failure is rare enough that I would not restructure for it. If you want the promise met
literally, the minimal change is to release the remaining guards before propagating — e.g. in
`parallel.rs`:

```rust
        for (i, guard) in pending_outputs.into_iter().enumerate() {
            if let Err(e) = guard.commit() {
                // Later outputs keep their temporaries rather than being deleted.
                return Err(e.context(format!("{i} of the run's outputs were already published")));
            }
        }
```

except that this does not work as written — the un-iterated guards are dropped by the `?`/`return`
anyway. Doing it properly needs a `release_all()` on the collection, which is more machinery than the
fault deserves. **My recommendation: leave the code and drop the "naming what committed and what did
not" clause from the plan**, since it is not what ships.

### F7 — Low: peak disk usage on an in-place re-run is now old + new

`File::create` used to truncate the previous output before the first record was read; now the previous
output survives until the rename, so a re-run in the same output directory needs room for both copies
at once. On a nearly-full volume a re-run that used to succeed can now hit `ENOSPC`. This is the exact
flip side of the property the PR is bought for — and for the Buckberry-scale files §3 uses to argue
for B, doubling peak footprint is the same order of consideration as the salvageable `.partial`. Worth
one clause in the PR description beside the "refusals are slower" note, not a docs change.

### F8 — Low: rename-to-publish replaces a symlink at the output path instead of writing through it

If a user pre-creates the output name as a symlink to another volume (a real, if obscure, trick for
keeping large outputs off a small filesystem), the old code followed the link and wrote to the target.
Now the bytes land in `.name.partial` beside the *link* — on the small volume — and the rename replaces
the link itself. No action; noting it because it is a semantic change that no test would catch. The
hard-link case changes in the opposite, better direction: a previously hard-linked copy keeps the old
content instead of silently seeing the new bytes.

### F9 — Low: `FastqWriter::flush` was a `pub fn` on a published crate

`trim-galore` publishes a library (`src/lib.rs`, no `publish = false`), so removing `pub fn flush` and
adding `pub fn finish` is a breaking change to the library surface at a patch/minor version. Practical
impact is near zero for a `command-line-utilities` crate, and the CHANGELOG should **not** carry it
(one sentence, user-facing). One clause in the PR body is right.

### F10 — Low: comment length, and one comment that claims too much

House rule is one line, current state, no rationale. Three of the new comments exceed it:

- `src/fastq.rs:648-658` — an eleven-line doc comment: two paragraphs of rationale plus the
  "teardown errors are unreachable" argument, which is the PR's justification rather than the API's
  contract, and contains the overstatement in F3. Proposed:

  ```rust
  /// Close the writer and publish the output under its final name. Consumes
  /// `self` so a caller cannot forget.
  ///
  /// `pending` comes out before `drop(self)`: the trailer is written by the
  /// inner writer's own `Drop`, and the final name must not appear before it.
  /// The guarantee is that the name appears only after the writer is closed,
  /// not that the bytes are complete — teardown errors are unreachable through
  /// `Box<dyn Write>`.
  ```

- `src/demux.rs:252-255` — four lines. Proposed: `// Publish in barcode order, then NoCode, with nothing fallible in between.`
- `src/bam.rs:648-651` — three lines; it reads as current-state invariant, so I would leave it, but
  it could collapse to `/// `try_finish` before the rename, and `pending` out before the drop.`

### F11 — Note: the SE parallel path still publishes a short output with exit 0 on a worker error

`parallel.rs:947-950` only `eprintln!`s a worker error and breaks the loop; the commit at `:991` then
publishes whatever was written, exit 0. Pre-existing (the same bytes reached the same final name
before), and out of scope. I flag it only because the new doc line on the *paired* function
(`:89-91`, "No output file is published unless every writer closes cleanly") is accurate there and
invites a reader to assume the SE twin has the same property. It does not. This belongs in the
already-planned §9 issue about swallowed teardown errors rather than a second one.

### F12 — Note: two concurrent runs at the same output path now share one temporary path

The derivation is deterministic, so both runs write `.x_trimmed.fq.partial`. The first to finish
renames it away; the second's rename then fails `ENOENT` while the message says "the trimmed data is
intact at <tmp>", which by then is false. Both runs corrupt each other today too (they share the final
path), so this is not a regression — but it is the one case where the new error message can lie.

---

## Plan-to-code audit

Everything in PLAN §5 is implemented, with these differences from the text:

| Plan says | Code does | Judgement |
|---|---|---|
| §5 step 1: `partial_output_name` gets a "no suffix may be appended to a committed name" note above it | The doc comment states the glob/extension reasoning instead | Equivalent; fine |
| §5 step 4: "**Do not remove the existing `flush()` calls** … They are what marks the commit point" | Kept everywhere it mattered (`clump_only`, `parallel`) — only `FastqWriter::flush` the *method* went | Correct reading of the instruction |
| §5 step 6 / §4 edge row: error "naming what committed and what did not" | Only the failing pair is named | F6 — recommend amending the plan |
| §7 T3: "no output file **and** no leftover temporary — assert both halves" | Both asserted, via `assert_nothing_published` | Done |
| §7 T7: confirm an existing all-filtered test asserts *presence* | Covered for `clump_only` (empty input reads back through `read_all`), not for a trim whose reads are all filtered | F4/task #32 |
| §7 T10 | `fastqc_report_is_named_from_the_published_output`, plus a no-dotfile assertion | Done, and stronger than asked |

Prior-review spot-check: PLAN_REVIEW_A/B's substantive findings did land in v2 — the commit-point
invariant (Assumption 5 → §5 step 4), `--clumpify`/`--hardtrim3` in §2.4, the `released`-before-rename
rule, the de-interleaver diagnosis, and O2's dot-prefix answer are all present in the code as well as
the prose.

## Issue-filing

Nothing here clears the bar for a new issue. F2 is the only finding with a user-visible symptom I can
name, and it is narrow enough to belong in the PR description. F11 folds into the already-planned §9
issue on swallowed teardown errors across `flate2`/`gzp`/`noodles-bgzf`. F1/F3/F4/F5 are edits to this
PR.
