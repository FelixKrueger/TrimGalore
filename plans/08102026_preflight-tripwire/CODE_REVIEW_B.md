# Code Review B — #414 preflight tripwire

**Reviewer:** Code Reviewer B (independent)
**Branch:** `feature/414-preflight-tripwire`
**Date:** 2026-08-10

## Scope

- `src/io.rs` — the hook (+68 lines, working tree)
- `tests/integration_preflight_tripwire.rs` — **51** `#[test]`s
- Read for context: `PLAN.md`, `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`, `CLAUDE.md`,
  `src/main.rs`, `src/cli.rs`, `src/parallel.rs`

> **Housekeeping, not a defect:** the branch has **no commits**. `git diff dev...HEAD`
> and `git log dev..HEAD` are both empty; the whole change is in the working tree
> (`M src/io.rs`, `?? tests/integration_preflight_tripwire.rs`). Anyone reviewing from
> the diff alone would see nothing.

> **The test file moved twice while I reviewed it.** Another editor in this shared tree
> landed three changes mid-review: the T3 hash set broadened from `guarded` to all of
> `before`; a new `rest.is_empty()` guard in `compare` plus its self-test
> (`comparison_rejects_an_allowance_that_takes_everything`, 50 → 51 tests); and the
> root-containment prefix gained a `MAIN_SEPARATOR`. All three are improvements, one of
> them closes a hole I had walked past (see "What I missed"), and one is a fix I was about
> to recommend. **Everything below is verified against the tree as of my final read;**
> line numbers are as-of that read and may drift again.

---

## Summary

The deliverable does what it claims. I did not take that on trust: I mutated
`src/main.rs` to add a writer with no candidate and confirmed the suite goes red with an
actionable message, and I walked every early return in `compare()` hunting for an input
that reaches `Ok(())` while a real miss exists.

The one substantive gap left is a level up from the assertion: **nothing detects a 12th
pre-flight site or a new dispatch arm.** `every_site_is_reached_by_some_case` proves 11
hand-listed shapes reach 11 *distinct* sites; it does not prove 11 is *all* of them. Its
own docstring — "turns 'all 11 sites are covered' from prose into an assertion" —
overstates what it asserts, and the gap is the same coverage-rot mechanism as the five
shipped defects, one level up. The fix is cheap and I verified the oracle it needs is
reliable (H1).

Below that: two Mediums about unasserted prose claims, and a short Low list.

**Verification I ran myself**

| Check | Result |
|---|---|
| `cargo test --test integration_preflight_tripwire` | **51 passed, 0 failed**, 11.1 s (debug) |
| same, before the concurrent edits | 50 passed, 0 failed, 9.6-10.4 s across three runs |
| `cargo test --release --test integration_preflight_tripwire` | 50 passed (pre-edit state), 1.33 s |
| `cargo test --release --lib dump_planned` | 1 passed |
| `cargo clippy --all-targets --release -- -D warnings` | clean, **21.9 s compile** — not a cache hit |
| Clippy really lints the test target | Planted `v.len() == 0`: clippy exits `error: could not compile trim-galore (test "integration_preflight_tripwire")`. Probe removed |
| **Mutation test** | red, with the right message — below |
| Release binary, env unset | 0 dump lines |
| Release binary, env `=1` | 3 dump lines, site `src/main.rs:1006` (**relative**) |
| `-V` / `--version` | `-V` has no provenance line; `--version` emits `3008b37 — macos/aarch64 — built …` |

### The mutation test — can the tripwire fail when it should?

Yes. I added a writer at `src/main.rs:1013` with no candidate:

```rust
run_single_file(&cli, input, &config, gzip, output_dir, &adapters_r1)?;
std::fs::write(input.with_extension("MUTANT"), b"x")?; // REVIEW-B MUTATION
```

`se_trim_two_inputs` went red:

```
se_trim_two_inputs: created but never planned: ["a.MUTANT", "b.MUTANT"]
Add each to the candidate list at src/main.rs:1006. Do NOT add them to the FastQC
allowance and do NOT narrow the candidate list — those are how this tripwire dies.
The 6 candidates that arm did plan: ["a.fastq_trimming_report.json",
"a.fastq_trimming_report.txt", "a_trimmed.fq", "b.fastq_trimming_report.json",
"b.fastq_trimming_report.txt", "b_trimmed.fq"]
```

Reverted; `git diff src/main.rs` empty. Two things beyond "it fails": the site is
reported **relative**, and the full candidate list ships with the failure, so the next
maintainer does not have to reconstruct the arm.

---

## What I missed

Worth recording, because it bears on how much weight to put on my "sound" verdict below.
I asked whether the FastQC allowance could swallow a *non-FastQC* write and concluded it
could not (only `fastqc::run` produces `*_fastqc.{html,zip}`). I did **not** ask the
other version of the question: what if the allowance swallows the *entire* created set?
`rest` goes empty, and T1 — `rest \ planned_keys` — is then trivially satisfied. The
concurrent editor closed exactly that, at `:355-361`:

```rust
// An allowance that swallowed the whole write set would leave T1 comparing empty sets.
if rest.is_empty() {
    return Err(format!(
        "the FastQC allowance took every created file ({allowed:?}), so T1 has nothing \
         left to check — narrow `is_fastqc_artifact`, do not widen it"
    ));
}
```

I checked the new guard for false reds and it is safe: `rest` can only be empty when
`active` is true (when `active` is false the partition puts everything in `rest`, and
`created` is already known non-empty), and no real arm writes *only* FastQC artefacts —
every FastQC-active case also writes trimmed output. Green across all 51.

---

## Where the code is sound

Checked specifically, hunting false negatives.

**`compare()`'s early returns — all seven walked, in the current form.**
`created.is_empty()` → `Err`; the two FastQC self-policing arms → `Err`; the new
`rest.is_empty()` → `Err`; `Arm::Unplanned` with a non-empty dump → `Err`;
`planned.is_empty()` on `Arm::Planned` → `Err`; two sites → `Err`; a candidate outside
the root → `Err`. The only unconditional `Ok(())` reached before T1 is the
`Arm::Unplanned` arm at `:371-381`, gated on `planned.is_empty()`. There is no argument
tuple that reaches `Ok(())` with a populated `rest` without T1 having run over it.

**No vacuous assertion left that I can find.** The three ways a subset relation goes
vacuous are each blocked: `created` empty (which matters because `main.rs:215-230` routes
clap's `DisplayHelpOnMissingArgumentOrSubcommand` to **stdout with exit 0**, so a case
whose positionals got swallowed by a flag would pass the exit-0 gate too), `planned`
empty, and now `rest` empty.

**`parse_planned()` under-yielding fails safe.** Dropping dump lines can only shrink
`planned`, which makes T1 *stricter*; dropping all of them trips the "hook did not fire"
branch. One theoretical exception, L3.

**`fastqc_active()` mirrors production; the plan's claim holds.** Verified independently
rather than from the plan: `cli.rs:337` is
`#[clap(long = "fastqc_args", allow_hyphen_values = true)]` with **no** `requires`, and
all five production gates (`main.rs:1458`, `:1777`, `:2241`, `:2398`, `:2521`) read
`cli.fastqc || cli.fastqc_args.is_some()`. So `--fastqc_args --quiet` genuinely runs
FastQC with `cli.fastqc == false`, and `se_trim_fastqc_args_only` tests a real property —
its non-empty-allowance assertion would have failed otherwise, and it passes.

**`snapshot()` sees what a run can write.** `DirEntry::file_type()` does not follow links,
so a created symlink lands in the file set rather than being skipped or recursed into;
files under directories created during the run are picked up by the stack walk;
`unwrap_or_else(panic)` on `read_dir` means an unreadable subtree fails loudly rather
than silently shrinking the set. The one thing it cannot see is a write outside the root
— M2.

**Temp-directory tags are unique.** Scripted over both `Tripwire::new` and bare
`case_root` call sites: 51 tag uses, 51 distinct, zero duplicates. Nine pairs are in a
*prefix* relation (`clump_only_se` / `clump_only_se_fastqc`, `hook_gate` /
`hook_gate_off`, `se_trim_fastqc` / `se_trim_fastqc_args_only`, …). Harmless:
`remove_dir_all` targets one exact path, and the `_{pid}` suffix sits between the tag and
any extension of it, so no case's directory name is a prefix of another's — the
shared-PID hazard the doc comment warns about does not bite.

**Both `parallel.rs` writer blocks are exercised — including the SE one.** I went looking
for the gap "plain SE `--cores > 1` is untested" and it is not one. `File::create(output)`
at `parallel.rs:956` is reached through the single `run_single_end_parallel` call site
(`main.rs:1321`), which `se_trim_clumpify` (`--clumpify --cores 2`) drives; the clumpy
path swaps the *reader* (`read_single_clumpy` vs `read_single_round_robin`), not the
writer. The paired block (`parallel.rs:212-216`, five handles) is covered by
`paired_trim_retain_unpaired_cores` and `paired_trim_passthrough_cores`. Everything at
`:1308`, `:2024-2026` and `:2109-2111` is inside `#[cfg(test)]` (module starts `:1217`),
so it is not production surface.

**`#[track_caller]` is sound and reproducibility-safe.** `preflight_output_collisions` is
only ever called directly — the 11 `main.rs` sites plus the `io.rs` unit tests — so there
is no fn-pointer or trait-object shim to reason about. `Location::file()` yields the
relative `src/main.rs`, confirmed empirically from the release binary, so nothing absolute
is embedded. I did **not** run the full double-build reproducibility job: the only new
embedded data are the `#[track_caller]` location constants (relative, verified) and two
`&'static str`s — no wall-clock, hostname or absolute path — so the invariant cannot be
threatened by this diff. Flagging that I reasoned rather than measured.

**The dump's placement is right.** It sits at the top of `preflight_output_collisions`,
before the collision checks, so a run that goes on to `bail!` still dumps its candidate
list — which is what makes the site reachable in a failure report.
`a_refused_run_writes_nothing` exercises that path and confirms the refusal still fires
with the hook on, and that nothing is written.

**`Arm::Unplanned` asserts more than I first credited.** Weak by construction — there is
no candidate list to compare against — but the live assertions are: exit 0, *something*
was created, *nothing* was dumped (so the day either orphan gains a pre-flight the case
goes red and says "promote the case"), and no pre-existing file changed. Note `hash_file`
panics on a missing file, so *deleting* a pre-existing file is caught too, not only
rewriting it. For those two arms that is a fair account of what can be asserted, and the
file says so plainly rather than implying more.

---

## Findings

### Critical

None.

### High

**H1 — Nothing detects a 12th pre-flight site, so the suite's central coverage claim can
go stale silently.** `tests/integration_preflight_tripwire.rs:1182-1267`:

```rust
/// Turns "all 11 sites are covered" from prose into an assertion. Collects the site
/// each representative arm dumps and requires 11 distinct values …
    let shapes: [(&str, &[&str]); 11] = [ … ];
    …
    assert_eq!(sites.len(), 11, "expected 11 distinct pre-flight sites, got {sites:?}");
```

What is asserted is *11 hand-listed shapes reach 11 distinct sites*. What the docstring
claims is *11 is all the sites there are*. Add a 12th `preflight_output_collisions` call
site to `main.rs` tomorrow — a new dispatch arm, or a split of an existing one — and this
test stays green while that arm has no case, no T1, and no tripwire. A new arm with *no*
pre-flight (the shape of the two existing orphans) is equally invisible; those two are
covered only because someone hand-listed them.

This is the mechanism behind the five shipped defects — coverage silently stops extending
to a writer — moved one level up, from "a candidate is missing from a list" to "an arm is
missing from the suite". Since resisting exactly that rot is the file's stated purpose,
the docstring's claim is worth making true rather than softening.

I verified the oracle a fix would need is reliable: `src/main.rs` mentions
`preflight_output_collisions` on exactly 11 lines and **all 11 are calls** — zero
doc-comment or prose mentions (`grep -c` → 11;
`grep -vc 'naming::preflight_output_collisions('` → 0). So:

```rust
/// The site table above is complete only if `main.rs` has no site the table misses.
#[test]
fn the_site_inventory_is_complete() {
    let src = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("src/main.rs"),
    )
    .unwrap();
    assert_eq!(
        src.matches("preflight_output_collisions(").count(),
        11,
        "main.rs gained or lost a pre-flight site; add a case for it to \
         every_site_is_reached_by_some_case and update this count"
    );
}
```

That closes the "12th site" half. The "new arm with no pre-flight" half cannot be closed
by a scan and belongs with the deferred follow-up for the orphan arms — worth one sentence
in "## Limits" so the gap is on record beside the case-folding one.

### Medium

**M1 — Each case's arm identity is unasserted, so a case can silently become redundant.**
The section headers are prose: `// ── site 2670 — run_specialty_paired, serving three
modes ──` covers `clock_multipair`, `implicon_pair`, `clump_only_paired_split_dirs` and
`clump_only_paired_no_report_file`; `// ── sites 661 / 706 / 752 / 809 — clump_only ──`
covers four more. Nothing checks that any individual case reached the site its header
claims. If a flag combination's routing changes — `--clump_only --paired` migrating
between `run_specialty_paired` and a clump-only site, say — the case keeps passing against
whichever arm it now lands on, the header still names the old one, and the arm it used to
be the only case for loses coverage.

`compare` already has the site in hand and reports it on failure; it just is not asserted.
Either (a) give `Tripwire` an `expect_site_of(&str)` group key and assert that all cases
sharing a key dump the same site — no line numbers, so unrelated `main.rs` edits cannot
break it; or (b) if that is more machinery than it is worth, reword the headers to say
"expected site" so nobody reads them as verified. What to avoid is leaving them as
unverified assertions in comments, since the file's whole argument is that prose claims
about coverage rot.

**M2 — A write outside the case root is invisible, and unlike the planned side there is no
guard.** `compare` explicitly protects the planned side (`:409-421`):

```rust
// A candidate resolving outside the case root would be invisible to `snapshot`,
// making T1 vacuous for it.
```

The created side has no equivalent. `snapshot(root)` only walks under `root`, so a writer
that escapes — an absolute `--output_dir`, `--fastqc_args "-o /abs/dir"`, a `..`-relative
default — writes a file no case can see, and T1 passes for it. No arm can escape *today*
(every fixture path in the suite is relative and `current_dir` is the root), which is why
this is Medium and not High. But "## Limits" documents the case-folding hole in five
careful lines and does not mention this one, and this is the likelier of the two to be
introduced by a future change. Either add it to "Limits", or assert in `run()` that no arg
is an absolute path — a two-line guard that would also catch a case accidentally written
against an absolute fixture path.

### Low

**L1 — Wording now lags the T3 broadening.** Two spots:

- `:293` `"{}: input {} was modified by the run"` — the hash set is every pre-existing
  file now, not only inputs. "pre-existing file" would be accurate.
- `:219-222` the `guarded` field doc still reads as the definition of the hashed set; its
  only remaining job is the staged-existence assertion at `:259-265`.

I did not apply these — see "Fixes applied".

**L2 — The hook activates on an *empty* env value.** `src/io.rs:178`:

```rust
if std::env::var_os(DUMP_PLANNED_ENV).is_some() {
```

Verified on the release binary: `TRIM_GALORE_DUMP_PLANNED_OUTPUTS= trim_galore a.fastq`
prints 3 dump lines. `FOO=` is a common idiom for "cleared" in shell wrappers and Docker
`ENV`, so a user could get stderr noise with no way to explain it. The doc comment does say
"Set this to any value", so this is as-specified rather than a bug; if you want the
stricter gate it is `.is_some_and(|v| !v.is_empty())`. Not a correctness risk — stderr
only, and it cannot change an output file.

**L3 — `parse_planned` drops a dump line glued to a preceding partial stderr line.**
`strip_prefix` requires the prefix at position 0, so if anything ever `eprint!`s without a
trailing newline immediately before the pre-flight, the first dump line is lost. This
fails *safe* everywhere except `Arm::Unplanned`, where losing all lines would make a newly
added pre-flight look absent. Nothing writes a partial line there today, and `eprintln!`
holds the stderr lock for the whole line so threads cannot interleave. Theoretical;
recorded because it is the one place under-parsing is not safe.

**L4 — RESOLVED during review.** Root containment was a bare `starts_with` on collision
keys with no separator check; the concurrent editor added
`let root_prefix = format!("{root_key}{}", std::path::MAIN_SEPARATOR);` at `:411`. Correct
fix, and it makes the check say what it means.

**L5 — `planned_keys` is built and then recomputed.** `:404-407` collects the key set,
then `:412-421` recomputes `collision_key(&root.join(p))` for each candidate. Reuse the
set, or build it after the loop. Immaterial at these sizes; it just reads as two passes
that should be one.

**L6 — `hook_is_inert_when_env_unset_and_that_goes_red` passes a snapshot where `compare`
expects a created set.** `:1177-1179`:

```rust
let created = snapshot(&root2);
compare(&root2, &created, &[], &["a.fastq"], Arm::Planned, true)
    .expect_err("no dump must fail the comparison, not pass it");
```

`snapshot` is the whole tree including the staged input, not `after \ before`. Harmless —
the `expect_err` is satisfied by the empty-dump branch regardless — but it reads as if the
third argument were a created set, and the next person to touch it may believe that.

**L7 — `case_root` silently degrades.** `:78` `std::fs::canonicalize(&d).unwrap_or(d)`,
where the doc comment two lines above explains that *not* canonicalising makes the
comparison pass vacuously on macOS. A silent fallback to the thing the comment calls
wrong should be `.expect("canonicalise case root")`.

**L8 — the hook is public API on a published crate.** `DUMP_PLANNED_ENV`,
`DUMP_PLANNED_PREFIX` and `dump_planned_lines` are `pub` in `src/io.rs`, and
`trim-galore` publishes to crates.io, so a test hook lands in the rendered docs and the
semver surface. `#[doc(hidden)]` on the three would keep them usable from `tests/` (same
crate) while keeping the hook out of public docs. Judgement call — the doc comments do say
"Test hook", which may be enough.

**L9 — no commits on the branch** (see the note at the top).

---

## Fixes applied

**None.** Three candidates and why I left each:

- The T3 hash-set broadening (`guarded` → all of `before`). I probed it — 51/51 pass at
  unchanged runtime — then reverted my probe, and a concurrent editor landed the same
  change with its own comment before I could recommend it. See R1.
- L4's separator fix — also landed concurrently.
- L1's two wording lines. One-word changes, but the file was under active concurrent edit
  and racing an `Edit` into it is how a fix gets silently lost.

Everything I did apply was a temporary probe, each reverted and verified clean: the
`main.rs` mutation, the T3 broadening probe, and the planted clippy violation.
`grep -rn "REVIEW-B\|MUTANT\|review_b_clippy" src tests` → no residue; `git diff src/main.rs`
empty; `git diff --stat src/io.rs` → 68 insertions, hook only.

---

## Recommendations, prioritised

1. **H1 — assert the site inventory is complete.** Add the source-scan test above (or
   equivalent) so a 12th pre-flight site cannot be added without touching this suite, plus
   one line in "## Limits" recording that a new arm with *no* pre-flight is still invisible.
   The only finding I would hold the branch for.
2. **M1 — assert or soften per-case arm identity.** Group cases by expected site and assert
   it, or reword the section headers so they do not read as verified.
3. **M2 — document (or guard) the outside-root blind spot** in "## Limits", beside the
   case-folding limit already there.
4. **L1 — two wording lines** to catch up with the T3 broadening.
5. L2, L5, L6, L7 at leisure; L3 and L8 are notes, not requests. L4 is done.

## On the plan-review fold-ins

Checked rather than assumed.

- **C1 (T2 unfalsifiable; promote T3).** Folded in — v1's T2 is gone and T3 is the
  byte-hash check. Note the plan is internally inconsistent about scope: §3.3's headline
  says "no **pre-existing file's** bytes changed (all 15 arms)" while §3.2 step 3 and §4's
  signature say the hashed set is **`guarded_inputs`**, and §14 records "`guarded_inputs`
  named explicitly as the hashed set". The implementation first followed §3.2 (narrower);
  with the concurrent edit it now follows §3.3 (broader), which is the stronger reading and
  matches PLAN_REVIEW_B's C1 wording — **R1**. I verified the broadening is free: all tests
  pass, runtime unchanged. It matters because the narrow version was complete only by the
  unasserted convention that every staged file appears in `guarded`; that convention held
  for all 41 arm cases when I checked, but nothing enforced it. Worth a line in §13 as a
  deviation, since §13 currently reads "Nothing in the plan was skipped".
- **C2 (put site identity in the dump).** Folded in via `#[track_caller]` +
  `Location::caller()`, and it works: the mutation test named `src/main.rs:1006`
  unprompted, with no hand-maintained arm→site table anywhere.
- **C3 (scope the FastQC non-empty invariant).** Folded in as `fastqc_capable`, stated in
  the module doc with the #421 link, and the stale-marker direction is asserted too
  (`(true, false) if !allowed.is_empty()`), which C3 did not ask for.
