# PLAN — `r1_unpaired` / `r2_unpaired` parity between the serial and parallel paths

**Issue:** [#544](https://github.com/FelixKrueger/TrimGalore/issues/544)
**Branch base:** `dev` @ `86aa301`
**Date:** 2026-11-11
**Revision:** 2 — see §12. Revision 1 specified a fixture that made two of its own validations
unsatisfiable, named a host for the new assertions that the CLI forbids, and asserted a gating
behaviour contradicted by its own measurements. All three were caught by dual plan review.
**Sibling:** #541 / PR #543 — same defect class, adjacent code.

---

## 1. Goal

Make `PairValidationStats::r1_unpaired` / `r2_unpaired` count reads the run **actually wrote**, on the
parallel path as on the serial one and as Perl 0.6.11 does, so `Unpaired read N kept` does not appear
at `--cores >= 2` for a run that wrote no unpaired file.

Secondary, and the reason this needs planning rather than patching: **both** tests that assert these
fields are blind to them, for two different reasons. Both have to change or the fix ships unverified.

Out of scope: #542, #529, and the `--retain_unpaired` rescue behaviour itself, which is correct.

## 2. Context

### The counters are a Perl-parity surface, not just an internal convention

Perl 0.6.11 (`git show 0.6.11:trim_galore`, lines 2325-2341) counts inside `if ($retain)`, and
increments immediately after the four `print UNPAIRED1` statements:

```perl
if ($retain){
    if ( length($seq_1) >= $length_read_1){
        print UNPAIRED1 $id_1;
        ...
        ++$read_1_printed;
    }
```

Its report line (`:2396`) reads `Number of unpaired read 1 reads printed`. So the quantity is
"records written", gated on rescue being requested, by upstream definition. The serial Rust path
mirrors it; the parallel path does not. **This is a departure from Perl, not merely an internal
inconsistency** — which matters, because Perl parity is a hard invariant here.

### Where the counters surface

| Consumer | Location | Gated on `> 0`? |
|---|---|---|
| Text report block | `src/report.rs:469-479` | **yes** |
| Summary lines | `src/main.rs:1902-1909`, `:2123-2130` | **yes** |
| JSON `r1_unpaired` / `r2_unpaired` | `src/report.rs:1073-1074` | **no** — `json_int` is unconditional |

Revision 1 claimed all consumers gate on `> 0`. The JSON does not: the keys are always present and
the fix changes a machine-read **value** from 4 to 0, while the text block stops being emitted.
Measured both ways. This distinction belongs in the CHANGELOG entry.

### The divergence

Serial (`src/trimmer.rs:584-595`) couples the counter to the write:

```rust
if let Some(ref mut w) = unpaired_r1.as_deref_mut()
    && r1_ok
{
    w.write_record(&r1)?;
    pair_stats.r1_unpaired += 1;
}
```

Parallel (`src/parallel.rs:571-576`) increments in `classify_paired` on the filter's verdict alone,
then returns `PairOutcome::Unpaired { r1_ok, r2_ok }` for its caller to act on. `classify_paired` is
never given the writer state, so it cannot condition on it.

### Why the fix has an unambiguous home

The call graph is narrow and fully traced; both reviewers re-verified it independently:

- `classify_paired` has **exactly one** caller, `process_pairs` (`src/parallel.rs:602`), and
  `PairOutcome::Unpaired` exactly one consumer.
- `process_pairs` takes `writer_up_r1: Option<&mut W>` / `writer_up_r2` (`:613-614`), and its
  `Unpaired` arm (`:642-662`) already writes under `if let Some(ref mut w) = writer_up_r1 && r1_ok`.
- Both `process_pairs` call sites (`:386` gzip, `:415` plain) pass those as `Some` exactly when
  `retain_unpaired` (`:370-379`, `:430-438`), which is `unpaired_r1_path.is_some()` (`:108`), traced
  back to `src/main.rs:1706-1713`.
- `process_paired_batch` (`:338`) has one caller, so **`--clumpify` cannot bypass `classify_paired`**.
- `run_paired_end_to_bam`'s own `TooShort` arm (`src/trimmer.rs:805-809`) neither rescues nor counts,
  so the uBAM path is unaffected.
- There is **no fourth increment site** anywhere in the crate.

So `Some(writer) && r*_ok` is exactly the serial predicate. Moving the two increments into that arm
fixes every route, with no signature change. `PairValidationStats::merge` (`src/report.rs:148-149`)
is a field-wise sum, so per-worker counting aggregates correctly.

### The filter can rescue both mates at once

`src/filters.rs:88-107` — and its own Perl-parity note citing `0.6.11:trim_galore:2325-2343` — is
explicit that when at least one mate fails `--length`, **each** read is routed independently on its
own length against `--length_1` / `--length_2`. There is no requirement that the other mate failed:

```rust
if r1_short || r2_short {
    return PairFilterResult::TooShort {
        r1_ok: r1.len() >= unpaired.r1,
        r2_ok: r2.len() >= unpaired.r2,
    };
}
```

So `r1_ok == r2_ok == true` is reachable and legitimate. Gating rescue on "the other mate failed" was
the beta.5 regression #245 fixed, which is why the fixture in §5 must include that shape — a
degenerate "fix" that reintroduced it would otherwise pass.

Note the two thresholds are distinct: `length_cutoff` (`--length`) discards the pair;
`unpaired.r1` / `unpaired.r2` (`--length_1` / `--length_2`) decide per-side rescue.

### Measured on `86aa301`

| Run on `BS-seq_10K_R{1,2}.fastq.gz` | `--cores 1` | `--cores 4` | unpaired files |
|---|---|---|---|
| `--paired` | lines absent | `read 2 kept: 4` | none |
| `--paired --clip_R1 5 --clip_R2 5` | lines absent | `read 2 kept: 4` | none |
| `--paired --clip_R1 5 --clip_R2 5 --retain_unpaired` | `read 2 kept: 4` | `read 2 kept: 4` | 2 written, R2 file holds exactly 4 records |

Rows 1 and 3 were independently re-derived by both reviewers. Row 3 settles it: 4 is the *correct*
rescue count and both core counts reach it, so the parallel path reports a true number for a question
nobody asked. Clipping is irrelevant (row 1 has no post-adapter trimming), so **any** paired run at
`--cores >= 2` where one mate survives the length filter and the other does not is affected.
`--retain_unpaired` runs are already correct at every core count, so the fix must not disturb them.

### Why no test catches it — there are two, not one

| Test | Line | Asserts | Blind because |
|---|---|---|---|
| `test_passthrough_serial_parallel_parity` | `:2360` | `PairValidationStats` at `:2338-2342` | `gen_multiome_rows` (`:2283`) makes **both** mates short, so `r1_ok == r2_ok == false` and neither path increments. Its doc comment (`:2275-2282`) says so outright. |
| `test_parallel_serial_bp_after_trim_parity_paired` | `:1622` | `PairValidationStats` at `:1705-1710` | `parity_test_config()` sets `length_cutoff: 0` (`:1433`), so the `TooShort` arm is **unreachable** and the counters are always 0. |

Revision 1 knew only about the first. The second was added by PR #543 two days ago, on a code
reviewer's suggestion — an assertion that cannot fail for the very fields this issue is about.

## 3. Behavior

1. On the parallel path, `r1_unpaired` increments only when an unpaired R1 writer is open **and**
   `r1_ok`; likewise R2. Identical to `src/trimmer.rs:584-595` and to Perl 0.6.11.
2. Without `--retain_unpaired`: both counters 0 at every `--cores`, the text block absent, both JSON
   values 0.
3. With `--retain_unpaired`: identical counts at every `--cores`, each equal to the records actually
   written to the corresponding unpaired file.
4. No other `PairValidationStats` field, no `TrimStats` field, and no output record changes.

### Edge cases

| Case | `r1_ok` / `r2_ok` | Expected |
|---|---|---|
| No `--retain_unpaired`, one mate rescuable | mixed | Writer `None`; no increment, no write, both paths |
| `--retain_unpaired`, only R1 rescuable | `true` / `false` | `r1_unpaired` increments; `r2_unpaired` does not |
| `--retain_unpaired`, only R2 rescuable | `false` / `true` | `r2_unpaired` increments; `r1_unpaired` does not |
| **`--retain_unpaired`, both rescuable** (#245 shape) | `true` / `true` | **Both** increment, both written — `src/filters.rs:88-107` |
| `--retain_unpaired`, neither rescuable | `false` / `false` | Neither increments — the case the current fixture pins |
| Both mates pass `--length` | — | Not `TooShort`; arm unreached |
| `--passthrough` + `--retain_unpaired` | — | **Rejected by `Cli::validate`** (`src/cli.rs:948-952`, test `:2021`) — verified, not assumed |
| `TooLong` / `TooManyN` | — | Do not construct `Unpaired`; untouched |
| `--cores 1` | — | Worker pool bypassed; serial already correct |

## 4. Signature

No signature changes. Two statements move between functions within `src/parallel.rs`.

## 5. Implementation outline

### Step 1 — remove the unconditional increments

`src/parallel.rs:571-576`: delete the two `if r1_ok { … += 1 }` / `if r2_ok { … += 1 }` blocks from
`classify_paired`'s `TooShort` arm. The arm still sets `pairs_removed`, `too_short`, and returns
`PairOutcome::Unpaired { r1_ok, r2_ok }` unchanged.

### Step 2 — count where the write happens

`src/parallel.rs:649-658`, inside each existing guard, after the write — matching both the serial
path and Perl's `++$read_1_printed` placement:

```rust
if let Some(ref mut w) = writer_up_r1
    && r1_ok
{
    r1.write_to(*w)?;
    pair_stats.r1_unpaired += 1;
}
```

and the same for R2.

### Step 3 — rewrite the arm's comment

The comment at `:643-648` explains that the rescue writes are no-ops without `--retain_unpaired`,
without noticing the counter had already fired upstream. Replace with one line:

```rust
// Counted only when the rescue writer is open, so the figure matches the records written.
```

### Step 4 — add a dedicated paired fixture with all four shapes

**Do not mutate `test_parallel_serial_bp_after_trim_parity_paired`.** Its `bp_after_trim` anchors
assume nothing is filtered (`length_cutoff: 0`); giving it a non-zero cutoff would invalidate them
and entangle #541's coverage with this fix. Add a **sibling** in the same module, reusing its
serial-versus-parallel driving structure (`run_paired_end` against `run_paired_end_parallel`, which
is non-passthrough and therefore may legally set `--retain_unpaired`).

The fixture needs **four** record shapes, because each exercises a different branch, and with
**unequal counts per shape** so that R1 and R2 expectations differ and a copy-paste error between
them fails:

| Shape | Construction | Yields | Exercises |
|---|---|---|---|
| A | R1 below `--length_1`, R2 long | `r1_ok=false, r2_ok=true` | `r2_unpaired` only |
| B | R2 below `--length_2`, R1 long | `r1_ok=true, r2_ok=false` | `r1_unpaired` only |
| C | both below `--length`, both above `--length_{1,2}` | `r1_ok=true, r2_ok=true` | both — the #245 shape |
| D | both below `--length_{1,2}` | `r1_ok=false, r2_ok=false` | neither |

Use a config with non-zero `length_cutoff` and explicit `UnpairedLengths`, chosen so each shape lands
in its intended branch. Revision 1 specified shape A alone, which leaves `r1_unpaired` at 0 on both
paths — making Validation 7 unsatisfiable and Step 6's R1 anchor `0`, a value a missing increment
also produces.

### Step 5 — assert under both rescue modes

In the new sibling, drive serial and parallel **twice** over the same fixture:

1. **Without `--retain_unpaired`** — assert `serial.2 == parallel.2` and that both counters are 0 on
   both paths. This is the regression under test.
2. **With `--retain_unpaired`** — assert `serial.2 == parallel.2` and that each counter equals the
   records actually written to its unpaired file.

The second run is what prevents a "fix" that simply zeroes the counters from passing as parity.
It cannot live in `test_passthrough_serial_parallel_parity`, because that combination is
CLI-rejected; the sibling host exists precisely to make it legal.

### Step 6 — anchor the expected counts

`serial == parallel` pins the paths to each other and to nothing true. Assert hand-derived absolute
counts from the shape table: `r1_unpaired` = (B + C counts), `r2_unpaired` = (A + C counts), with the
two deliberately unequal. Derive them from the fixture definition, never by reading back a run.

### Step 7 — guard the fixture's discriminating power

```rust
assert!(
    serial.2.r1_unpaired > 0 && serial.2.r2_unpaired > 0,
    "fixture must rescue on both sides, or this test cannot fail"
);
```

In the `--retain_unpaired` run, where non-zero is correct. Requiring **both** to be non-zero is what
catches a future fixture edit that silently drops shape A or B.

### Step 8 — un-blind the passthrough fixture

Both reviewers independently checked all five `gen_multiome_rows` callers (`:2396`, `:2477`, `:2589`,
`:2633`, `:2680`) and found **none** depends on the both-mates-short shape; lengthening R2 on a
short-R1 row changes only `r2_ok`, not the pair verdict or the drop rate. So edit the shared
generator rather than adding a variant — a variant would leave
`test_passthrough_serial_parallel_parity` blind, which is the state this plan exists to end.

Since passthrough forbids `--retain_unpaired`, that test then pins the **no-rescue** case: counters
must be 0 on both paths even though mates are rescuable. Delete the now-false claim from the doc
comment at `:2275-2282` and state what the fixture does.

### Step 9 — CHANGELOG

One sentence under `### Unreleased`, in the file's format, mirrored into
`docs/src/content/docs/reference/changelog.md`. It must name **both** effects: the text block stops
appearing, and the JSON values change from a spurious count to 0.

### Step 10 — mutation checks

Revert the R1 increment alone, then the R2 increment alone, then **both**. The first two must each
fail; the third must fail the `--retain_unpaired` assertion, or "always zero" would pass as parity.

## 6. Efficiency

No change. Two `usize` increments move inside guards already evaluated for the writes.

## 7. Integration

**Reads:** `writer_up_r1` / `writer_up_r2` and `r1_ok` / `r2_ok`, all already in scope.
**Writes:** `pair_stats.r1_unpaired` / `r2_unpaired` only.
**Downstream:** for no-rescue runs at `--cores >= 2` the fix **removes** the text report block
(`src/report.rs:469-479`, `src/main.rs:1902-1909`, `:2123-2130`) and **changes the JSON values**
4 → 0 (`src/report.rs:1073-1074`, ungated). Both bring `--cores >= 2` into line with `--cores 1` and
with Perl.

**Perl parity:** the CI matrix runs `--cores 1`, so it exercises the already-correct path and stays
green either way. It proves no collateral damage, not the fix.

**Step 8 touches a shared fixture** — blast radius verified nil by two reviewers, but re-run the full
suite rather than assuming.

## 8. Assumptions

1. The counters mean "records written to the unpaired files", gated on rescue being requested.
   **Verified against Perl 0.6.11** (`:2325-2341`, `:2396`), not inferred from our own code.
   **Fixed rule.**
2. `writer_up_r*` is `Some` exactly when `--retain_unpaired`. Verified `:108`, `:370-379`,
   `:430-438`, `src/main.rs:1706-1713`.
3. `classify_paired` has one caller, `PairOutcome::Unpaired` one consumer, no fourth increment site.
   Verified by both reviewers.
4. `PairValidationStats::merge` is a field-wise sum. Verified `src/report.rs:148-149`.
5. `--passthrough` with `--retain_unpaired` is rejected at the CLI layer. **Verified**
   (`src/cli.rs:948-952`, test `:2021`) — revision 1 took this from a comment.
6. Editing `gen_multiome_rows` does not disturb its other four callers. Verified by both reviewers;
   re-run the suite regardless.

## 9. Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| 1 | No-rescue run core-invariant | `--paired`, `--cores 1` and `4` | Text block **absent** from both; JSON both 0 |
| 2 | Same with clipping | add `--clip_R1 5 --clip_R2 5` | As above |
| 3 | Rescue run unchanged | add `--retain_unpaired`, `--cores 1` and `4` | `read 2 kept: 4` at both; 2 files written |
| 4 | Counts match the files | `gzip -dc` each unpaired file, count records | Equals the reported figure |
| 5 | JSON values move | compare `r1_unpaired`/`r2_unpaired` before and after, no-rescue at `--cores 4` | 0/4 → 0/0, keys present throughout |
| 6 | Records unchanged | SHA-256 of trimmed R1/R2 before and after | Identical |
| 7 | **R1 site covered** | Revert the R1 increment alone | A test fails — requires shape B or C in the fixture |
| 8 | **R2 site covered** | Revert the R2 increment alone | A test fails — requires shape A or C |
| 9 | "Always zero" is not parity | Revert **both** | The `--retain_unpaired` assertion fails |
| 10 | #245 not reintroduced | Shape C present and asserted | Both counters increment for one record |
| 11 | Suite and gates | `cargo test`; `fmt --check`; `clippy --all-targets --release -D warnings` | 845+ pass, both clean (local rustc 1.95 vs CI 1.99 — necessary, not sufficient) |

Validations 7–9 are the ones revision 1 could not have satisfied.

## 10. Questions or ambiguities

**Open (assumption taken, does not block):**

- Exact `length_cutoff` and `UnpairedLengths` values for the Step 4 fixture. Any set that puts the
  four shapes in their intended branches works; pick round numbers and record the derivation beside
  the fixture.
- Whether Step 8's edit should also add an assertion to the passthrough test that the counters stay
  0, or rely on its existing whole-struct comparison. Taken: the existing comparison suffices once
  the fixture can produce rescuable mates.
- Whether the docs site describes these counters anywhere needing a matching edit. Taken: no — but
  grep `docs/` for "unpaired" before finishing.

**No critical ambiguities.** Semantics are fixed by Perl, the call graph is traced, every assumption
is now verified rather than quoted.

## 11. Self-Review

**Efficiency.** Nil impact.

**Logic.** The fix is a move, not a reimplementation: the destination guard is already the serial
predicate. Confirmed the uBAM and clumpify routes cannot bypass it.

**Adjusted in revision 2.** The test plan changed more than the fix did — see §12. The fix itself
(Steps 1–3) is unchanged from revision 1 and was confirmed correct by both reviewers.

**Edge cases.** §3 now includes the both-rescued #245 shape, which revision 1 omitted and which a
degenerate fix would otherwise slip past.

**Integration.** The fix removes report lines *and* changes JSON values; revision 1 asserted only the
former, on a false claim about gating.

**Remaining risks.**

1. The Step 4 fixture must put four shapes in four branches. If a shape lands in the wrong branch the
   test may pass while covering less than intended — which is why Step 7 guards both sides and
   Validations 7–8 revert each site separately.
2. CI cannot catch a regression here; Steps 4–8 are the whole guard.
3. Step 8 edits a fixture shared by five tests. Two reviewers checked it; the suite is the backstop.

## 12. Implementation notes

Implemented 2026-11-11 against `dev` @ `86aa301`. All ten steps done. Diff confined to three files,
per `git diff --numstat`: `src/parallel.rs` (+203/−32), `CHANGELOG.md` (+7),
`docs/.../reference/changelog.md` (+7). **`src/trimmer.rs` is untouched** — it was mutated during the
checks below and reverted exactly, which `git diff` confirms.

| Step | Outcome |
|---|---|
| 1 — remove unconditional increments | `classify_paired`'s `TooShort` arm now only sets `pairs_removed` / `too_short` |
| 2 — count at the write | Increments added inside both `if let Some(w) && r*_ok` guards in `process_pairs` |
| 3 — comment | One sentence, citing #544 and Perl 0.6.11 (wraps to two physical lines) |
| 4 — fixture | New sibling `test_parallel_serial_unpaired_counter_parity`; four shapes (3×A, 5×B, 2×C, 4×D) plus 6 passing pairs, `length_cutoff: 30`, `UnpairedLengths { r1: 20, r2: 20 }` |
| 5 — both rescue modes | Loop over `[false, true]`, asserting `serial.2 == parallel.2` in each |
| 6 — anchors | Hand-derived `r1_unpaired = B+C = 7`, `r2_unpaired = A+C = 5`, deliberately unequal |
| 7 — guard | Both counters non-zero **and** unequal, in the rescue run |
| 8 — un-blind the passthrough fixture | `gen_multiome_rows` R2 now always 48 bp; doc comment rewritten |
| 9 — CHANGELOG | One sentence naming both the vanished text block and the JSON values, mirrored to the docs copy |
| 10 — mutation checks | Four performed; see the table below, plus a fifth added after review |

### Lengths chosen

`seq(n)` cycles `ACGT`, and every length avoids `n % 4 == 1` so no sequence ends in `A` — which at
`min_overlap: 1` would match the Illumina adapter's first base and perturb the lengths. Hence 26 bp
for shape C rather than 25.

### Mutation checks

| Mutation | Result |
|---|---|
| Parallel R1 increment removed | **Failed** — parity assertion, retain mode |
| Parallel R2 increment removed | **Failed** — parity assertion; `r2_unpaired` serial 5 vs parallel 0 |
| **All four increments removed (serial *and* parallel)** | Parity assertion **passed** at 0/0; the **absolute anchor** fired. This is the degenerate "zero everything" fix, and only Step 6 catches it. |
| Original bug reintroduced in `classify_paired` | **Two** tests failed — the new one *and* `test_passthrough_serial_parallel_parity`, which was blind before Step 8. Proves the un-blinding worked. |

The third and fourth are the ones that justify Steps 6 and 8 respectively.

**Deviation from Step 10 as written.** The plan said to revert "both" parallel increments. Done, but
it proves little: serial still counts, so the parity assertion fails exactly as the single-side
reverts do. The informative mutation removes all four, which makes `serial == parallel` hold at 0/0
so that only the anchor can fire. Step 10's wording understated what was needed; the stronger check
was run.

### Verification

| Check | Result |
|---|---|
| V1 / V2 — no-rescue, text block | **Absent** at `--cores 1` and `4`, with and without clipping |
| V3 — rescue run | `Unpaired read 2 kept: 4` at both core counts |
| V4 — counts match files | report 4, file 4 records, both core counts |
| V5 — JSON | keys present, both values 0 at `--cores 1` and `4` (was 0/4 at `--cores 4`) |
| V6 — records unchanged | `53d619cd…`, 9,996 records, identical to pre-fix |
| V11 — suite and gates | `cargo test` exit 0, **845 passed / 0 failed** (844 + 1 new); `clippy --all-targets --release -D warnings` exit 0, 0 warnings; `fmt --check` exit 0 |

Clippy initially failed with two `needless_borrows_for_generic_args` in the new test (`&dir.join(…)`
passed to a generic parameter); both fixed, then all three gates re-run clean.

**Not done, deliberately:** no commit, branch or PR — not requested. Local rustc is 1.95 against CI's
1.99, so the clippy pass is necessary but not sufficient.

**Clarification on §1 vs Step 4.** §1 says "both tests have to change"; Step 4 says not to mutate
`test_parallel_serial_bp_after_trim_parity_paired`. Both hold: the two *blind* tests are the
passthrough one (changed in Step 8) and the #541 paired one — and the latter is un-blinded not by
editing it but by a sibling that covers the branch it cannot reach. #541's anchors stay undisturbed.

**§10's `docs/` grep, recorded:** run across `docs/`, `wiki/` and `README.md`. The only matches for
`Unpaired read`, `r1_unpaired` or `r2_unpaired` are the changelog lines this change adds, so no
documentation edit is needed. The assumption taken in §10 was correct.

### Post-review pass

Dual code review (`CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`) and the coverage audit (`COVERAGE.md`,
**Verdict: COMPLETE**) all recommended the change as-is: **no Critical, no High**. Four changes were
applied from their findings and the full gated set re-run.

| Change | Source |
|---|---|
| `gen_multiome_rows`'s R2 now ends in `T` (`GCAT`-cycled, still 48 bp) | **B.** The previous `TGCA`-cycled string ended in `A`, so the adapter's leading `A` matched a one-base overlap and the filter saw **47** bp where the comment claimed 48. Verified directly: that record trims 48 → 47 with `Reads with adapters: 1 (100.0%)`. The trailing-`A` discipline was written into the new test's helper and then not applied to the fixture it edited. |
| `length_1` → `length_2` in that comment | **Both reviewers.** R2's rescue threshold is `--length_2`; the wrong name contradicted the correct sibling comment four lines below. |
| The fixture-discrimination guard replaced by two `const` assertions on the shape totals | **B (L-1).** The old `assert!(r1 > 0 && r2 > 0 && r1 != r2)` ran *after* the 7/5 anchors had passed, so it was unconditionally true — the "assertion that cannot fail" pattern this plan indicts in §2. The replacement is evaluated at compile time against the constants themselves; verified it fires by setting the two totals equal, which fails the build with `evaluation panicked: shapes must rescue on both sides with unequal totals`. |
| Test scaled to 4,200 pairs and looped over `gzip in [false, true]`, plus an assertion that no rescue file exists without `--retain_unpaired` | **B (M-1).** At 20 pairs with `gzip = false` the test covered neither the `GzEncoder` rescue sinks nor cross-batch `PairValidationStats::merge` of these fields — i.e. what a real `--paired --retain_unpaired --cores 4` run does. `PAIRS > BATCH_SIZE` is now itself a `const` assertion. |

Re-verified afterwards: `cargo test` 845 passed / 0 failed; `clippy --all-targets --release
-D warnings` exit 0 with 0 warnings; `fmt --check` exit 0; `src/trimmer.rs` still unmodified; diff
confined to the same three files (`src/parallel.rs` +211/−32, both changelogs +7).

### Follow-up found during review — recorded, not filed

**`--clumpify --retain_unpaired` is legal and untested** (Reviewer A). It is rejected in combination
with uBAM output, `--clump_only` and `--passthrough`, but not with `--clumpify` — which forces
`--cores >= 2` and therefore always takes the path this issue concerns. The clumpify tests pass
`clump_layout: Some` with unpaired paths `None`; the new test does the inverse, so the combination is
uncovered. A verified it by hand at 7/5 with 7 and 5 records on disk, the same multiset as the
non-clumpify run in a different order, as the clusterer intends. Not a defect; a gap.

### Harness traps worth knowing (Reviewer A)

- The `Unpaired read N kept` block is written to the **R2** report. Grepping the R1 report returns
  nothing even on a correct rescue run, which reads as absence.
- The package is `trim-galore` (hyphen) while the binary is `trim_galore` (underscore). Grepping a
  clippy log for the underscore form makes a clean run look like a false green.

## 13. Revision History

**Revision 2 — 2026-11-11.** Dual plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`). Both
reviewers verified the call graph independently and both found revision 1's test plan broken in the
same three places. Every correction below was re-verified against the source before being applied.

| Change | Why |
|---|---|
| **Step 4: four fixture shapes with unequal counts**, replacing "some records have only R1 short" | **Critical, both reviewers.** Shape A alone leaves `r1_ok == false` for every pair reaching the arm, so `r1_unpaired` stays 0 on both paths: Validation 7 could not fail, and Step 6's R1 anchor would be `0` — a value a missing increment also produces. |
| **Steps 4–7 move to a new sibling of `test_parallel_serial_bp_after_trim_parity_paired`**, not the passthrough test | **Both reviewers.** `--retain_unpaired` with `--passthrough` is CLI-rejected, so revision 1's Step 5 named a host that cannot run it, and §10 offered an "implementer's choice" that did not exist. The sibling is non-passthrough. Not mutating the existing test, whose `bp_after_trim` anchors assume `length_cutoff: 0`. |
| **§2: the JSON is not gated on `> 0`** | **Both reviewers.** `src/report.rs:1073-1074` emits unconditionally. Revision 1 claimed otherwise while its own measurements showed the keys present with value 0. Changes what Step 9 must say. |
| §2: a **second** blind parity test identified (`:1622`, assertion `:1705-1710`, blind via `length_cutoff: 0`) | **Both reviewers.** Added by PR #543 two days ago on a reviewer's suggestion — an assertion that cannot fail for these fields. Revision 1 knew only of the passthrough one. |
| §2: Perl 0.6.11 grounding added (`:2325-2341`, `:2396`) | Reviewer A. Perl counts inside `if ($retain)`, increments after the `print`, and calls the figure "printed" — so Assumption 1 is upstream-sanctioned and this is a **Perl-parity** defect, not an internal inconsistency. |
| §3: **both-rescued (#245) shape added**; Step 4 shape C; Validation 10 | Reviewer B. `src/filters.rs:88-107` makes `r1_ok == r2_ok == true` reachable. Gating rescue on "the other mate failed" was the beta.5 regression #245 fixed; without this shape a fix reintroducing it would pass. |
| §2: `process_paired_batch` and `run_paired_end_to_bam` routes closed; "no fourth increment site" | Reviewer B, and A respectively. Confirms `--clumpify` and uBAM cannot bypass the fix. |
| Step 8: **edit `gen_multiome_rows` directly**, dropping revision 1's hedge toward a variant | **Both reviewers** checked all five callers: blast radius nil. The hedge would have left the passthrough test blind — the state this plan exists to end. |
| Assumption 5 promoted from unverified to verified | Both reviewers (`src/cli.rs:948-952`, test `:2021`). Revision 1 quoted a source comment. |

Steps 1–3 — the fix — are unchanged and were confirmed correct by both reviewers.
