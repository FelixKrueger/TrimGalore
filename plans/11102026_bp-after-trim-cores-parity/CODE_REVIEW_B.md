# CODE REVIEW B — #541 `total_bp_after_trim` parity between serial and parallel paths

**Reviewer:** B (independent; no contact with Reviewer A or the plan-manager)
**Target:** uncommitted working-tree change at `/Users/fkrueger/Github/TrimGalore`, branch `dev`, base `c0bbee4`
**Diff:** 3 files, +222/−5 — `src/parallel.rs`, `CHANGELOG.md`, `docs/src/content/docs/reference/changelog.md`
**Date:** 2026-11-10

---

## Summary

**The change is sound. No Critical or High findings. I recommend it as-is**, with six Low notes
below, none of which blocks, plus one Medium finding (**M1**) that is *not* part of this change: a
second `--cores`-dependent figure in the paired report, which I hit while verifying Validation 2.
It should not be fixed inside this diff.

The fix is the minimal correct one: three expressions change value, nothing moves, and
`TrimResult` was already in scope at all three sites. I re-derived the plan's claims rather than
accepting them, and every one I checked held. Two of the brief's four "weakest evidence" questions
had a cheap empirical answer, so I answered them with the binary rather than with the tests that
are themselves under review:

- **The absolute anchors are right.** On a synthetic fixture identical to the tests' own, the
  release binary reports `Total written (filtered): 960 bp` (SE) and `960` / `990` (PE R1/R2) at
  both `--cores 1` and `--cores 4`, with `Total basepairs processed: 1,155`. `15 × 30 + 15 × 34 =
  960` and `30 × 33 = 990` are therefore correct, independently of either code path.
- **The R2 fixture picks up no spurious adapter match.** The R2 trimming report reads
  `Trimmed: 0 times` at `error_rate 0.1` / `min_overlap 1`, and R2's `Total written (filtered)`
  equals its `Total basepairs processed` (990 = 990). The SE/R1 fixture reads `Trimmed: 15 times`,
  i.e. the 15 clean 34 bp reads also match nothing. Adapter trimming really does yield exactly
  30 bp.
- **The whole report, not just the one line, is now core-count invariant.** On
  `BS-seq_10K_R1.fastq.gz --clip_R1 5` and `SRR24766921_RRBS_R1.fastq.gz --rrbs`, the full
  `*_trimming_report.txt` is byte-identical between `--cores 1` and `--cores 4` (`diff` exit 0),
  and the JSON differs only in its `command_line` field. That is a stronger statement than the
  implementer's single-figure check.

I also reproduced all four mutation checks in an isolated copy of the source tree (not the shared
working tree), and each failed as claimed, in the right test and at the right assertion.

### What I verified by running

| Check | Result |
|---|---|
| `cargo test` (own `CARGO_TARGET_DIR`) | exit 0 — **844 passed, 0 failed, 0 ignored**; both new tests appear by name in the log and were not filtered out |
| `cargo fmt --all -- --check` | exit 0, empty output |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0, **zero** warnings — and the invocation provably lints this crate's own code, see "Notes on method" |
| Validation 1 — SE figure | `644,525 bp (99.2%)` and `490,527 bp (96.2%)` at **both** `--cores 1` and `4` |
| Validation 2 — PE figure | `--paired --clip_R1 5 --clip_R2 5` on `BS-seq_10K_R{1,2}`: R1 `644,525`, R2 `645,073` at both core counts (my own run, matching the implementer's); synthetic PE R1 `960` / R2 `990` at both |
| Validation 3 — records | `--rrbs` output decompressed and `cmp`-identical between core counts (9,729 records); synthetic SE output `cmp`-identical |
| Validation 4 — `total_bp_written` | JSON `594,545` at both core counts |
| Mutation `:1157` → `record.seq.len()` | SE test **fails** at `parallel.rs:1587` (whole-struct parity) |
| Mutation `:542` → `r1.seq.len()` | PE test **fails at the R1 assertion**, `parallel.rs:1693` |
| Mutation `:543` → `r2.seq.len()` | PE test **fails at the R2 assertion**, `parallel.rs:1699` |
| Mutation serial **and** parallel together | parity passes (both 720); the **absolute anchor fires** at `parallel.rs:1598`, `left: 720 / right: 960` |
| Adjacent check — full paired report across core counts | R1 report byte-identical; **R2 report differs by three lines** → finding **M1** |
| Working tree after my checks | `git diff --stat` still 3 files, +222/−5 (mutations were applied to a copy under my scratchpad, never to `/Users/fkrueger/Github/TrimGalore`) |

### What I verified by reading

- **The three changed sites are exhaustive for the parallel path.** `src/parallel.rs` has exactly
  three `trimmer::trim_read` call sites (`:495`, `:496`, `:1133`); `classify_paired` and
  `classify_single` have exactly one caller each (`process_pairs` at `:628`, `process_reads` at
  `:1194`); `--clumpify` changes only the *reader* dispatch (`read_single_clumpy` / `read_pairs_clumpy`)
  and `--passthrough` only adds a write, so both still funnel through the same two functions.
  Repo-wide, `total_bp_after_trim` is accumulated at nine sites — six serial in `src/trimmer.rs`
  (`:350`, `:536`, `:537`, `:672`, `:780`, `:781`) and these three — and all nine now read
  `bp_after_cutadapt`. `src/demux.rs` and `src/specialty.rs` touch neither field.
- **Placement is unchanged and matches the serial paths.** Both families accumulate after the
  `--discard_untrimmed` early return and before `filters::filter_*`; `TrimStats::merge`
  (`report.rs:70`) is a field-wise sum, so a per-read value change cannot double-count across
  batches.
- **`total_bp_written` cannot be perturbed.** It is written only in the filter `Pass` arms
  (`parallel.rs:554-555`, `:1166`), which the diff does not touch, and the new tests assert
  whole-struct `TrimStats` equality, so a change to it would fail them.
- **No consumer subtracts or re-derives the field.** Its only reads are `report.rs:635-636`
  (`format_number` + `percentage` against `total_bp_processed`). Since `bp_after_cutadapt ≤`
  post-quality length `≤` input length, the percentage stays ≤ 100% — same bound as before — and
  there is no subtraction anywhere that could underflow.
- **Empty `adapters_r2` is handled in exactly one place.** The fallback lives inside `trim_read`
  (`trimmer.rs:114-118`), and both `run_paired_end` (`trimmer.rs:484`) and the parallel
  `classify_paired` (`parallel.rs:496`) call `trim_read(r2, config, true)`, so there is no
  serial/parallel divergence to have. R2 is in fact scanned with the **R1** adapter list; the
  fixture simply contains no match for it (confirmed by `Trimmed: 0 times`).
- **No golden fixture needs updating.** `grep -rn "Total written" tests/` returns nothing; no test
  asserts report text or a report hash, so the plan's "no fixture changes" claim holds.
- **The CHANGELOG mirror is exact.** `CHANGELOG.md` from line 3 and
  `docs/src/content/docs/reference/changelog.md` from line 9 are byte-identical (`diff` exit 0),
  including the new entry.
- **Every flag the CHANGELOG entry names exists.** `--polyA` and `--polyG` are real clap aliases
  (`cli.rs:370`, `:382`), `--trim-n` is the canonical long form (`cli.rs:118`), and `--clumpify`
  does require `--cores >= 2` (`cli.rs:744-747`).

---

## Issues by area

### Logic

**Nothing wrong found in the fix itself.** The three sites now read the same quantity the serial
paths read, from a `TrimResult` that was already bound; `bp_after_cutadapt` is a `Copy` `usize` and
the records are not re-read, so there is no borrow or aliasing consequence. **M1** below is a
pre-existing defect in the neighbouring lines, not a consequence of this change.

**L1 — the new comment's shorthand is narrower than the quantity it describes**
`src/parallel.rs:541`, `:1156`

```rust
// Post-adapter, pre-RRBS length, so the cutadapt section matches Perl v0.6.x (#232).
```

The captured value (`trimmer.rs:165`) precedes RRBS **and** poly-G, poly-A, N-trimming **and**
fixed clipping. "pre-RRBS" names only the first of those five. That matters here specifically
because the discriminator in the two new tests is *clipping*, not RRBS: a reader who checks the
comment against the test it is guarded by sees `clip_r1`/`clip_r2` in the fixture and a comment
that says "pre-RRBS", and the two do not obviously describe the same invariant. The test doc
comments use a third shorthand, "post-adapter, pre-clip" (`:1542`, `:1597`).

*Failure scenario:* someone adds a sixth post-adapter stage, reads "pre-RRBS" as exhaustive, and
puts the new stage above the capture point at `trimmer.rs:165`.

*Suggestion (style only, still one line):*
`// Post-adapter length, before RRBS, poly-X, N-trim and clipping (#232).`

Confidence: **high** that it is imprecise; **low** that it will cause harm, because the full
contract is one `TrimResult` field away at `trimmer.rs:64-71`.

**L3 — the discrimination guard admits a second way to be satisfied**
`src/parallel.rs:1609-1612`, `:1715-1720`

```rust
assert!(stats_serial.total_bp_after_trim > stats_serial.total_bp_written, ...);
```

`total_bp_written` excludes filtered reads, so this inequality is satisfied by post-adapter
trimming **or** by any read being filtered out. Today `parity_test_config()` sets
`length_cutoff: 0` (`parallel.rs:1433`) and nothing is filtered (I confirmed 30/30 records written
on both sides), so it does mean clipping. If a future edit gave `parity_test_config()` a non-zero `length_cutoff`, the
guard would keep passing with the four clip values removed — the exact re-vacuuming it exists to
prevent.

This is **not a hole**, because the absolute anchors (960 / 990) would still fail in that scenario;
the guard is redundancy. I note it only because the plan's Step 7 framed the guard as the thing
that stops re-vacuuming, and it is the anchors that actually do. Worth recording that the
implementer's choice of `> total_bp_written` is *stronger* than the plan's suggested
`< total_bp_processed`, which adapter trimming alone would satisfy.

Confidence: **high** (verified `length_cutoff: 0` at `parallel.rs:1433` and 30/30 records written).

**M1 — a second `--cores`-dependent report figure survives, five lines below the edited ones
(pre-existing, out of scope, but reportable)**
`src/parallel.rs:571-576` vs `src/trimmer.rs:584-595`

Verifying Validation 2 on the real fixture turned this up. With plain
`--paired --clip_R1 5 --clip_R2 5` on `BS-seq_10K_R{1,2}.fastq.gz` and **no `--retain_unpaired`**,
the R2 trimming report gains three lines at `--cores 4` that are absent at `--cores 1`:

```
Unpaired read 1 kept:                 0
Unpaired read 2 kept:                 4
```

and the JSON agrees: `"r2_unpaired": 4` at `--cores 4`, `0` at `--cores 1` (`pairs_removed: 4` on
both). No `*_unpaired_2.fq.gz` exists in either output directory — nothing was kept, because no
rescue writer was ever opened.

Mechanism: the parallel `TooShort` arm increments `r1_unpaired` / `r2_unpaired` from the filter's
`r1_ok` / `r2_ok` unconditionally (`parallel.rs:571-576`), while the serial arm increments only
inside `if let Some(ref mut w) = unpaired_r1 … && r1_ok` (`trimmer.rs:584-595`). `report.rs:469`
gates the whole block on `> 0`, so the line appears or vanishes with the core count.

This is **not introduced by this change** — the diff does not touch `pair_stats` — and it is
already documented as known at `parallel.rs:2269-2276`, where the passthrough fixture deliberately
makes both mates short so `r1_ok == r2_ok == false` and the parity assertion survives. I raise it
because it is the same defect class as #541 (a trimming-report figure that moves with `--cores`),
it lives in the same function, and unlike #541 it asserts something that is simply false: reads
"kept" when no unpaired file was written. It also explains why **L6** matters.

*Recommendation:* out of scope for this change — do not widen the diff. Worth its own issue, which
is Felix's call under the issue-filing bar; the reproducer above is complete, and the fix is to
gate the two increments on the rescue writers as the serial path does. I filed nothing.

Confidence: **high** — reproduced at both core counts, in both the text and JSON reports, with the
mechanism read in both paths.

### Efficiency

**No change, confirmed.** `bp_after_cutadapt` is a `usize` already computed and stored in
`TrimResult` for every read on every path; the read replaces a `String::len()` call. No
allocation, no extra pass, no change in complexity or memory. The two new tests add ~0.05 s to the
lib suite (the filtered run of both completed in 0.01–0.05 s).

### Errors

**None found.** No new panics, no unwraps, no overflow or underflow path (see "no consumer
subtracts" above), no change to error propagation. `fs::write` / `FastqWriter::create` in the new
tests use `?` consistently with the module.

### Structure

**L2 — the SE test duplicates the existing parity test's fixture generator verbatim**
`src/parallel.rs:1551-1564` vs `:1482-1495`

The 30-read, two-sequence loop is character-for-character the same in both tests (I diffed them).
A shared `fn write_parity_fixture(path: &Path)` would keep the clean-config run and the
clipped-config run on *the same* input, which is the relationship that makes the pair meaningful.
As written, a future edit to one fixture silently decouples them, and the clean-config test keeps
passing on a different input than the one the clipped test anchors.

*Failure scenario:* someone lengthens the clean read to 40 bp in one test to probe a length filter;
the two tests now disagree about what "the parity fixture" is, and `EXPECTED_BP_AFTER_TRIM` in the
other test no longer corresponds to anything the clean-config test exercises.

Confidence: **high** (duplication verified by reading both).

**L4 — both new tests exercise a single batch**

30 reads against `BATCH_SIZE = 4096` (`parallel.rs:33`) means one batch per run, so for this field
the tests pin per-read accumulation and a single merge, not multi-batch accumulation. The existing
`test_parallel_serial_trim_stats_parity` says as much about itself (`:1478-1481`). `TrimStats::merge`
is a plain sum and my 10,000-read end-to-end runs (3 batches) came out byte-identical, so this is
empirically fine — but no *unit* test pins it.

*Named input that would close it:* the same clipped config over a fixture of >4,096 reads. Cheap:
raise the loop bound in the SE test from 30 to 5,000 and scale the two constants, or add a third
test.

Confidence: **high** (`BATCH_SIZE` read; batch count inferred from record counts).

**L6 — the paired test does not compare `PairValidationStats`**
`src/parallel.rs:1692-1703`

It asserts `serial.0 == parallel.0` and `serial.1 == parallel.1` but leaves `.2` uncompared. With
`length_cutoff: 0` the `TooShort` arm is unreachable, so asserting `.2` today would pass — but see
**M1**: it would *not* pass under a non-zero length cutoff, for a reason that has nothing to do
with this change. Adding the third `assert_eq!` as the fixture stands is one line and costs
nothing; adding it along with a length cutoff would turn the test red.

### Comment and changelog style (text this change introduces)

The two source comments are one line each, state the current fact, and do not narrate the change —
compliant with the project rule, modulo **L1** above.

The test doc comments run 3–5 lines each. They state invariants and justify the fixture rather than
recounting the bug's history, and they match the register of the surrounding test module (the
existing parity test carries a 7-line doc comment). I read this as within house style.

Two inline comments are worth a second look, both in the "justification in source" category the
project rule pushes toward the commit message:

- `:1594-1596` — "Hand-derived: an expectation read back from either path would reproduce the
  defect it exists to catch." This is an instruction to the next editor (*do not regenerate this
  constant from output*) attached to the thing it constrains, which the rule explicitly allows.
  Defensible; I would keep it.
- `:1453-1455` (the helper doc) — "Both R2 values are required: `stats_r2.total_bp_after_trim` is
  a separate site." Same shape: a forward-looking invariant, not history. Keep.

**L5 — the CHANGELOG entry carries root-cause rationale**
`CHANGELOG.md:8-14` and its mirror

The entry is a bold headline plus one ~72-word sentence that explains *why* the figure moved ("The
parallel path accumulated the final record length into `total_bp_after_trim` while the serial paths
used the post-adapter, pre-clip length introduced for #232…"). The global rule is one sentence per
entry with rationale going to the commit message or PR thread.

Against that: the file's own house style is multi-sentence paragraphs that carry exactly this kind
of rationale, and this entry is **shorter than every 2.4.0 entry above it**. "Follow the file's
existing format exactly" and "one sentence" pull in opposite directions here, and the entry
satisfies the former. I am flagging the tension, not asking for a rewrite. If you do trim it, the
clause that earns its place for a reader deciding whether the upgrade affects them is the flag list
and the `--clumpify` consequence; the serial-vs-parallel mechanism is the part that belongs in the
commit message.

One small internal inconsistency: the entry spells the flags `--polyA` / `--polyG`, while the 2.4.0
entries in the same file use `--poly_a` / `--poly_g`. Both are valid (`cli.rs:370`, `:382`); only
one is the canonical long form.

Confidence: **high** on the facts; the priority call is editorial and yours.

---

## What these tests still cannot catch

Named states, not generic gaps. None of these is introduced by this change.

1. **`--discard_untrimmed` at `--cores >= 2`.** The contract that a discarded read contributes to
   *neither* field depends on the accumulation sitting below the early return (`parallel.rs:533-538`,
   `:1150-1154`). Both new tests set `discard_untrimmed: false`, so moving the accumulation above
   that return would leave them green. Input that would catch it: the clipped config with
   `discard_untrimmed: true` over the existing mixed adapter-bearing/clean fixture, asserting
   `total_bp_after_trim == 15 × 30` rather than 960.
2. **`--rrbs --paired --cores 4`.** The new configs set `rrbs: false`, so the one combination that
   motivated `bp_after_cutadapt` in the first place — the RRBS 2 bp truncation, plus the auto-set
   `clip_r2 = 2` and the `rrbs_r2_clipped_5prime` counter that commit 3996fc5 fixed — has no
   serial-vs-parallel parity test. Structurally low-risk (the capture point is in the shared
   `trim_read`, pinned by `trimmer.rs:1076`), and my end-to-end `--rrbs` run was byte-identical
   across core counts, but single-end only.
3. **uBAM output at `--cores >= 2`.** The new tests drive the FASTQ entry points only. uBAM output
   is single-threaded today (`--cores` silently ignored), so there is nothing to diverge — but if
   the pluggable-I/O epic ever parallelises it with its own accumulation, these tests will not see
   it. `#529` as it stands routes through `classify_*` and inherits the fix, which matches the
   plan's §7.
4. **A read consumed entirely by adapter trimming** (`bp_after_cutadapt == 0`) is absent from both
   fixtures. The arithmetic is trivial (adds 0 on both paths) and `filter_single_end` with
   `length_cutoff: 0` would pass a zero-length record, so this is a low-value gap — but it is the
   one input shape where "post-adapter length" and "final length" coincide at zero while every
   other read diverges.

---

## Recommendations

**Critical:** none.

**High:** none.

**Medium** — not a defect in this change, and **not to be fixed inside it**:

- **M1** — `r1_unpaired` / `r2_unpaired` also move with `--cores`, and at `--cores >= 2` the paired
  report claims unpaired reads were "kept" when `--retain_unpaired` was never passed and no
  unpaired file exists. Pre-existing, already noted at `parallel.rs:2269-2276`. Deserves its own
  issue; reproducer in §Logic.

**Low** — all optional; the change can ship without any of them:

1. **L1** — widen the two one-line comments (`parallel.rs:541`, `:1156`) from "pre-RRBS" to the full
   set of post-capture stages, or point at `TrimResult::bp_after_cutadapt`. Cheapest fix, best
   ratio of the six.
2. **L2** — factor the duplicated 30-read fixture generator into one helper shared by
   `test_parallel_serial_trim_stats_parity` and `test_parallel_serial_bp_after_trim_parity_with_clipping`.
3. **L4** — if you want multi-batch coverage for this field, raise one fixture above
   `BATCH_SIZE = 4096` and scale the constants.
4. **L5** — editorial call on the CHANGELOG entry's second sentence; at minimum consider
   `--poly_a` / `--poly_g` for consistency with the rest of the file.
5. **L3** — optionally tighten the discrimination guard to name the clip values directly (e.g.
   `assert!(config.clip_r1.is_some() && config.clip_r2.is_some())`), or leave it: the absolute
   anchors already carry that weight.
6. **L6** — one extra `assert_eq!(serial.2, parallel.2)` in the paired test. It passes on today's
   fixture and would pin `PairValidationStats` against future drift — but do **not** pair it with a
   length cutoff until **M1** is fixed, or the test goes red for an unrelated reason.

If you take only one, take **L1**.

---

## Notes on method

- I ran everything under my own `CARGO_TARGET_DIR` in the session scratchpad, never the
  repository's `target/`. The repo binary was older than `src/parallel.rs` (22:06 vs 22:08), so I
  built my own rather than trusting it.
- All mutation checks were applied to a **copy** of `src/`, `Cargo.toml`, `Cargo.lock`, `build.rs`,
  `tests/` and `examples/` under the scratchpad, so the shared working tree was never mutated. It
  is unchanged: `git diff --stat` reports 3 files, +222/−5, and `git status --porcelain` lists only
  the three expected files.
- **The clippy result is not a skipped unit.** `cargo clippy` prints only one status line per
  package, and for this crate it reads `Compiling trim-galore` (the build script needs codegen)
  rather than `Checking trim-galore`, so the log alone does not prove the lib was linted. I proved
  it in the copy: injecting `let _probe = result.bp_after_cutadapt.clone();` next to `:1157` makes
  the identical invocation fail with `error: using clone on type usize which implements the Copy
  trait … could not compile trim-galore (lib)`. The repository's exit 0 with zero warnings is
  therefore a real pass. (Probe removed; the copy was restored.)
- My fourth mutation (`stats.total_bp_after_trim += record.seq.len()`) matched by text and so also
  hit the uBAM serial site `trimmer.rs:672`, one line more than the implementer's. The conclusion
  is unaffected: the parity assertion passed at 720 on both sides and the absolute anchor fired.
- Per-run figures I report as "verified" came from reading the `*_trimming_report.txt` /
  `.json` files produced by my own build, not from the plan.
- I wrote no source file. The only file I created in the repository is this report.
