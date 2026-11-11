# PLAN_REVIEW_B — #414 empirical exact-set tripwire

**Reviewer:** Plan Reviewer B (independent, fresh context)
**Target:** `plans/08102026_preflight-tripwire/PLAN.md` (v1, 2026-08-10)
**Tree:** `dev` @ `3008b37` — every line reference below re-read from that tree
**Verdict:** sound mechanism, unusually well evidenced. Two Critical defects — one
advertised assertion is unfalsifiable, and the dump lacks the site identity that §3.4
and the coverage claim both need. Seven Important items, all cheap. Nothing here
challenges a settled decision.
**Status:** COMPLETE.

---

## 0. Verification log — what I checked and what held

Stated plainly, because several of the plan's load-bearing claims turned out to be
exactly right and should not be re-litigated by the implementer.

**Verified correct:**

- **The 11 sites and their line numbers.** `grep -n preflight_output_collisions
  src/main.rs` at `3008b37` returns exactly 466, 497, 661, 706, 752, 809, 927, 1006,
  2047, 2109, 2670. Every number in §2.4's table is right.
- **Site 2670 serves three modes.** `run_specialty_paired` is `main.rs:2654`; the
  pre-flight is at `:2670`; callers are `--clock` (`:527`), `--implicon` (`:549`) and
  `--clump_only --paired` (`:605`). §2.4's claim holds.
- **A1 — there is no third write-capable dispatch path.** I traced every leaf of
  `main()` rather than trusting the grep: `hardtrim5`→466, `hardtrim3`→497,
  `clock`→2670, `implicon`→2670, `clump_only`{FASTQ:{paired→2670, SE→661},
  UBam:{paired N=1→706, paired N>1→752, SE→809}}, `--output-format ubam`→
  `run_ubam_output`{paired N=1→**orphan**, paired→2047, SE→2109}, `paired` N=1→
  **orphan**, `paired`→927, SE→1006. Every non-test `File::create` /
  `FastqWriter::create` / `BamWriter::create` in `src/` sits downstream of one of
  those 15. The two shapes that could have been a third orphan are closed:
  `--clump_only --paired` with one input is rejected at `main.rs:597`, and
  `--clock`/`--implicon` go through `validate_paired_input` (`cli.rs:1004`, `:1007`),
  which requires an even input count. A1 is sound.
- **The matrix does reach all 11 sites.** 466→cases 18,19 · 497→20,21 · 661→24,26 ·
  706→28 · 752→29 · 809→27 · 927→11–17 · 1006→1–10 · 2047→32 · 2109→30,31 ·
  2670→22,23,25 · orphans→33,34. The row count is 34 and matches V5.
- **`ci.yml:69` / `:76`** are exactly `cargo nextest run` and `cargo test --release`.
  A4's reasoning is right and is the single most important assumption in the plan.
- **The `root.join(p)` arrangement is not just right, it is necessary.**
  `collision_key` → `lexical_normalise` → `std::path::absolute` (`io.rs:47`)
  absolutises against the **test process's** CWD, not the child's. A raw
  `collision_key(dumped_relative_path)` would key every candidate against the wrong
  root and T1 would pass vacuously. §3.2 step 1 (canonicalised root) plus §3.3
  (`root.join(p)`) close this. Good catch by the plan author.
- **FastQC's write set.** On every probed arm it is exactly `<stem>_fastqc.html` +
  `<stem>_fastqc.zip`; `-o DIR` relocates both; and fastqc-rust creates `DIR`
  itself (probe `fqc_odir`, exit 0, `qc/` present), so case 4 needs no `mkdir`.
- **Hook inertness.** Nothing in `src/` reads the environment today except
  `env_logger::init()` (`main.rs:207`). With the var unset the dump is one `var_os`
  and no output, so the `validation` md5 job and the `reproducibility` job are
  genuinely untouched. The dump lands on stderr before any thread is spawned or any
  file opened, so there is no interleaving hazard from `--cores` workers or FastQC.
- **Exactly one pre-flight call per run, on every pre-flighted arm.** All 11 sites
  sit outside their per-pair loops (`2666-2670`, `869-931`, `731-752`, `989-1006`,
  `2026-2047`, `2090-2109`). `parse_planned` therefore never has to disambiguate
  interleaved dumps from two arms. **This is an invariant the parse silently relies
  on and the plan does not state — worth writing down, because the day someone moves
  a pre-flight inside a loop the parse starts unioning two arms' lists and T1 gets
  quietly weaker.**
- **§2.3's orphan-arm safety conclusion is correct** — see O-b for why the stated
  reason is not the load-bearing one.

---

## 1. Logic review

### C1 (Critical) — T2 can never fail. It is vacuous under every case in the matrix.

§3.3:

> **T2 — no created path aliases an input** (all 15 arms): for every created path,
> `collision_key` differs from every input's `collision_key`. Independent of the
> pre-flight's own check, because it looks at what was actually *written* rather
> than at what was planned. This is the assertion that gives the two orphan arms
> real coverage.

`created` is defined in §3.2 step 6 as `after \ before`. An output that overwrites
an input is written to a path that **was already in `before`** — inputs exist before
the run by definition. So it is not in `created`, and T2 cannot see it. The one
scenario T2 claims to guard is precisely the one its own set arithmetic excludes.

The only residual way T2 can fire is a *case-only* alias on a case-sensitive
filesystem (input `Sample.fq`, output written to `sample.fq` — a new path, so in
`created`, and fold-equal under `collision_key`). No case in the matrix constructs
one; and on the 13 pre-flighted arms a case that did would be refused by the
pre-flight, so it would fail §3.2 step 5's exit-0 requirement rather than reach T2.

Consequences:

- T3 is doing all the work, on all 15 arms. The sentence "This is the assertion that
  gives the two orphan arms real coverage" is false; T3 is.
- The orphan arms' actual coverage reduces to: the run exits 0, `planned` is empty,
  and no input's bytes changed. That is thin but it is what was decided — the defect
  is the plan asserting more than it delivers, not the decision.

**Fix.** Compute the written set properly and run T2 over it:

```
modified = { p ∈ before : hash(p) changed }
written  = created ∪ modified
```

T2 over `written` is non-vacuous, directly catches #409's harm, and subsumes T3
(inside a fresh case root the only pre-existing files *are* the inputs, so
`modified ≠ ∅` ⟺ an input was clobbered). Either restate T2 on `written` and keep
T3 as the named special case, or delete T2 and promote T3 to "no pre-existing file's
bytes changed" — which is strictly stronger than the current pair and one assertion
shorter.

### C2 (Critical) — the dump carries no site identity, so §3.4 cannot be implemented and a case can silently exercise the wrong arm.

Two consequences of one gap.

**(a) §3.4's failure message is unimplementable as specified.** It requires:

> `<path>` was created but never planned. Add it to this arm's candidate list in
> `src/main.rs` (site `<N>`).

§3.1's dump format is `TRIM_GALORE_PLANNED\t<path>` — no site, no arm. So `<N>` has
to come from a hand-written arm→site table in the test file. That table is exactly
what §1 forbids:

> The expectation is derived by **observation**, never by hand-enumeration —
> hand-enumerating relocates the original knowledge problem into a test file, where
> it rots the same way.

And it rots on the first edit to `main.rs`: all eleven numbers shift.

**(b) Nothing binds a case to the arm it is meant to exercise.** Sites 706 and 752
differ only by input count; 2670 serves three modes; 466 and 497 by one flag. If a
fixture or flag slip routes case 28 to site 752, T1 still passes — against 752's
candidate list — and the plan's central claim ("all 11 arms") becomes untrue at
runtime with no signal at all. This is a false-negative on the *coverage* claim, not
on any single assertion, which makes it the harder kind to notice.

**Fix — one change closes both, and it touches none of the 11 call sites.** Put
`#[track_caller]` on `preflight_output_collisions` and emit the caller in the dump:

```rust
#[track_caller]
pub fn preflight_output_collisions(...) -> Result<()> {
    if std::env::var_os(DUMP_PLANNED_ENV).is_some() {
        let at = std::panic::Location::caller();
        for (path, _) in planned {
            eprintln!("{}{}:{}\t{}", DUMP_PLANNED_PREFIX, at.file(), at.line(), path.display());
        }
    }
    ...
```

This gives exact `src/main.rs:466`-style attribution for free — the candidate-
assembly code stays untouched, which is the constraint §2.1 is built around.
`Location::file()` is the source-relative path cargo passes to rustc, so it is
deterministic and the reproducibility job is unaffected. Stable since 1.46; the
crate floor is 1.88 (`Cargo.toml:5`).

Then: have each case assert *which* site dumped (recorded once per case, updated by
the compiler-free but loud route of a failing assert), or — to avoid pinning line
numbers at all — assert only that **all dumped lines in one run carry the same
site**, and print the site in the T1 failure message. The second form gives §3.4 its
`<N>` with zero maintenance and still catches the "two arms' lists got unioned"
regression flagged in §0.

### C3 (Important, but it undercuts §12's central claim) — the FastQC self-policing rule is false on four of the fifteen arms.

§3.3 states the rule unconditionally:

> when `--fastqc` is set the allowed set must be **non-empty** (otherwise FastQC
> never ran and the allowance is masking a broken case)

`src/specialty.rs` contains **no reference to `fastqc` at all** (grep is empty), and
`run_paired_ubam_single_file` (orphan arm 1, `main.rs:1810-2003`) has no
`fastqc::run` call either — the complete set is `main.rs:1459`, `:1778`, `:1784`,
`:1796`, `:2242`, `:2399`, `:2522` plus `clump_only.rs:374`, `:554`, `:555`, `:847`,
`:1095`. Nor does `--fastqc` carry any `conflicts_with` (`cli.rs:330-331`). So
`--hardtrim5 20 --fastqc a.fastq` is accepted, exits 0, and writes zero FastQC files.

On arms 1, 2, 5, 6 and orphan 14, therefore, the invariant is false. No case in the
matrix trips it today, but the first person to add `hardtrim5_fastqc` or
`clock_fastqc` gets a red test for a reason that has nothing to do with output
collisions — and the obvious fix is to weaken the rule. That is the exact death §12
says §3.4 exists to prevent:

> The failure mode is not a false negative, it is a *maintainer weakening the
> allowance* when a legitimate new writer appears. Three things guard it: […] it
> must be non-empty when active and empty when not […]

**Fix.** Scope the invariant to the arms where `fastqc::run` is reachable, and say so
in the test: "`--fastqc` is a silent no-op on `--hardtrim5/3`, `--clock`,
`--implicon` and on single-file interleaved paired FASTQ output; a case must not pair
them." One sentence, and it turns a landmine into documentation. (Whether the silent
no-op is itself a bug is a separate question and out of scope for #414.)

### I1 (Important) — V2's negative-control edit, taken literally, does not compile and over-deletes.

V2:

> Comment out the report candidates at `main.rs:919-925`

The tree at `3008b37`:

```
918            if !cli.no_report_file {
919                for input in [&chunk[0], &chunk[1]] {
920                    let src = naming::OutputSource::Input(input.clone());
921                    candidates.push((naming::report_name(input, Some(&pair_dir)), src.clone()));
922                    candidates.push((naming::json_report_name(input, Some(&pair_dir)), src));
923                }
924            }
925            planned.extend(candidates);
```

Two problems with `919-925` as an *instruction*:

- It starts at 919, so once 924 goes the `if !cli.no_report_file {` at 918 is
  unclosed — brace mismatch, no build.
- It ends at 925, which is `planned.extend(candidates);`, not a report candidate.
  Deleting that drops **every** candidate for the arm, primaries included, so
  `planned` is empty and T1 fails naming all six created paths. V2's stated
  expectation — "Fails naming all four report paths" — is then not what you observe,
  and that difference is precisely the signal that tells you the control did what you
  meant rather than something broader.

The correct range is **918–924**, and with it V2's expectation holds exactly: four
report candidates disappear while `write_paired_reports` still writes all four files,
so T1 fails naming four paths.

Worth fixing because V2 is the single step that proves the tripwire can go red at
all — a compile error there is the worst place for the implementer to lose time.

**Explicitly not a finding:** §10's other five citations are sound. The plan is
consistently citing the *namer-call lines* rather than enclosing blocks —
`1745-1751` and `2369-2375` each span exactly the four
`report_name`/`json_report_name` calls, `2040-2045` spans the loop plus its braces,
and `1962` / `2492` are exactly the first hand-rolled `format!` in each block. Only
`919-925`'s end is off by one, and only V2 turns it into something you execute.

### I2 (Important) — no non-emptiness floor, so T1's liveness is guarded only by accident.

Nothing asserts `planned` is non-empty on an `Arm::Planned` case, and nothing asserts
`created` is non-empty on any case. Today the arithmetic saves you: if the hook dies
or `parse_planned` breaks, `planned` is empty, `created \ allowed` is not, and T1
fails. But that is an accident of the current shape, not a designed guard — and it
stops holding the moment a case's run writes only allowed FastQC files, or writes
nothing at all.

Make both explicit:

- `Arm::Planned` ⇒ `planned` non-empty (the hook fired and the parse worked).
- all arms ⇒ `created` non-empty (the run actually did something).

Both are one line, and together they make "this tripwire is live" a property the
suite asserts rather than one it happens to exhibit — the point of
`feedback_harness_false_negatives`, which §5 step 5 already invokes for `compare`
but not for the wiring around it.

### I3 (Important) — §5 step 2's unit test cannot be written as described.

> **`src/io.rs` unit test** — `dump_env_emits_one_line_per_candidate`: with the var
> set, a two-candidate list produces two prefixed lines; unset, none.

Three problems, compounding:

1. Edition is 2024 (`Cargo.toml:4`), so `std::env::set_var` is an `unsafe fn`.
2. `.github/workflows/ci.yml:55-56` records the opposite as a maintained invariant:
   *"no `static mut` / `OnceCell` / `env::set_var` globals exist — the codebase is
   nextest-clean."* Adding one makes that comment false.
3. Under `cargo test --release` all unit tests share one process and run as parallel
   threads, so the variable leaks into the ~15 other `preflight_output_collisions`
   unit tests at `io.rs:1129-1310`. They only inspect `Result`, so nothing fails
   today — a silent, invariant-breaking change is worse than a loud one.

The plan's own fallback ("checking the dumped-line count matches the known candidate
count for that arm") reintroduces a hand-written count.

**Fix.** Do V1 entirely in the integration file, with no in-process env mutation and
no `unsafe`: run one arm **twice**, once with `DUMP_PLANNED_ENV` set and once
without, and assert prefix-line count `> 0` and `== 0` respectively. That is a
stronger test than the unit version — it proves the hook fires *through the real
binary*, which is the only place it matters — and the "inert when unset" half is the
half that actually needs asserting, since the fires-when-set half is implied by 32
green cases.

### I4 (Important) — O1 excludes `--cores`, but `--cores > 1` selects a different writer implementation.

O1:

> `--rename`, `--length_1/2`, `--cores`, and `--dont_gzip` are excluded because they
> add no writer.

`main.rs:1316`, `:1548` and `:1847` all branch `if cli.cores > 1 || cli.clumpify`
into `parallel.rs`, whose paired writer opens its own five handles at
`parallel.rs:212-216` — primaries, passthrough, and both unpaired files — entirely
separate code from the serial block at `main.rs:1586-1604`. `--cores` does not add a
writer; it swaps in a **second writer implementation** behind the same candidate
list.

The default is `--cores 1` (`cli.rs:347`), and every paired case in the matrix runs
at the default. Cases 13 (`--retain_unpaired`) and 14 (`--passthrough`) — the two
that exercise the most conditional writers — therefore only ever observe the serial
block. Two writer implementations kept in step with one hand-assembled candidate
list by human diligence is precisely the drift shape the tripwire exists to police,
and it is the half the matrix does not look at.

Add `--cores 2` variants of cases 13 and 14 (two extra cases), or simply run the
paired cases at `--cores 2`. Cost is negligible on 40-read fixtures.

### I5 (Important) — a duplicated case tag can turn a real miss into a pass.

`case_root(tag)` inherits `integration_output_collision.rs::tempdir` (`:31-36`),
which names the directory `tg_coll_{tag}_{process::id()}` and opens with
`remove_dir_all`. Under `cargo test --release` all 34 cases share one PID and run as
parallel threads, so uniqueness rests entirely on `tag`.

A duplicated tag does not merely collide — one case deletes another's tree *mid-run*.
And because `created = after \ before`, a deletion that removes an unplanned file
converts a genuine T1 miss into a pass. With 34 hand-written tags that is a live
risk. State that tags must be unique (ideally derive them from the test function
name), or key the directory on tag + thread id.

### I6 (Important) — T2/T3's input set must be `guarded_inputs`, not `cli.input`.

§3.2 step 3 says "a content hash of every input" without naming the set.
`guarded_inputs` (`main.rs:74-83`) is `cli.input` **plus** the `--passthrough` carrier
**plus** the `--demux` barcode file. Cases 5, 6, 14 and 16 have one of each, and a
barcode file or carrier clobbered by a namer is exactly the harm class #389 and #409
were. Say `guarded_inputs` explicitly, or the implementer will reasonably read
"input" as the positional list.

---

## 2. Assumptions

A1–A9 all check out on the facts, with two refinements.

- **A1** — verified independently; see §0. Sound.
- **A2** — "`--fastqc` is the only source of unplanned writes." True across the
  probed set. The residual is `--fastqc_args "--svg"`, which sets
  `config.svg_output` (`fastqc.rs:81`) and was not probed; see O-e. Also note A2 is
  only "fixed by the self-policing allowance" on the arms where FastQC runs — C3.
- **A3, A5, A6, A7, A9** — sound. A9's two fixtures are both in `git ls-files
  test_files/` (`ubam_test.bam`, `ubam_paired_test.bam`). ✔
- **A4** — the most important assumption in the plan, and correct: `ci.yml:69`
  (`cargo nextest run`, debug) and `:76` (`cargo test --release`). A `cfg` gate would
  make the release run a false pass. ✔
- **A8** — correct as stated, with an unstated limit worth recording; see O-a.

**Implicit assumptions the plan does not surface:**

1. *Exactly one pre-flight invocation per run.* True on all 11 sites today (§0), and
   `parse_planned`'s "union all prefixed lines" silently depends on it.
2. *The dump is complete.* The hook sits at the top of the pre-flight, so it dumps
   even when the pre-flight goes on to `bail!`. Harmless here (§3.2 step 5 requires
   exit 0), but it means a dump is not evidence the run was accepted.
3. *`--basename` is not honoured on the orphan arms.* Verified — both use
   `Path::file_stem()` unconditionally (`main.rs:1820`, `:2424`) and never read
   `cli.basename`. This is load-bearing for §2.3 and is not written down; see O-b.

---

## 3. Efficiency

No concerns. 34 child processes on 40-read fixtures, five of them running bundled
FastQC; `snapshot` walks ≤ ~20 files; comparisons are `BTreeSet` over
`collision_key` strings. `.config/nextest.toml` sets `slow-timeout = 60s,
terminate-after = 4`, so a hung case is killed at 240 s rather than hanging the job.
Two notes:

- The default is `--cores 1` (`cli.rs:347`), so 34 concurrent single-threaded
  children will not oversubscribe a 2-core runner. If I4's `--cores 2` variants land,
  that is +2 children at 2 threads — still fine.
- Everything runs twice (debug nextest + release `cargo test`). The debug binary is
  markedly slower, but at 40 reads that is noise. O2's contingency list (drop cases
  6, 16, 26, 31) is reasonable and correctly protects 3 and 4.

---

## 4. Validation sufficiency

This is the area the plan cares most about, so I went looking specifically for ways
the tripwire could pass when it should fail.

**Found (detailed above):**

| | False-negative path | Severity |
|---|---|---|
| C1 | T2 cannot fire — an aliasing output is never in `created` | Critical |
| C2b | A case silently exercises a different arm; T1 passes against that arm's list | Critical |
| I2 | T1 vacuous if `created \ allowed` is empty; `planned` liveness unasserted | Important |
| I4 | The parallel paired writer implementation is never observed | Important |
| I5 | A duplicated case tag deletes a sibling's tree, hiding an unplanned file | Important |

**Checked and found *not* to be false-negative paths:**

- *Hook dies / parse breaks → empty `planned` → trivial pass?* No: T1 fails, because
  `created \ allowed` is non-empty. Fail-safe direction. (But by accident — I2.)
- *Files written outside the snapshot root?* No. Every namer is either CWD-anchored
  (`specialty::hardtrim_output_name` returns a bare filename when `output_dir` is
  `None`, `specialty.rs:452-455`) or input-parent-anchored, and all inputs and all
  `-o` targets in the matrix live inside the root.
- *`/tmp` vs `/private/tmp` passing vacuously?* Closed by §3.2 step 1, and the
  reason is stated correctly.
- *stderr interleaving corrupting a dumped line?* No — the dump precedes every
  thread spawn and every file open.
- *`-o qc` failing because `qc` does not exist?* No — fastqc-rust creates it (probe
  `fqc_odir`, exit 0).
- *Directories counted as writes?* Closed by A7; `ensure_output_dir` runs at
  `main.rs:420`, before all dispatch, so `-o out` exists before the pre-flight.
- *`--clumpify` temp files?* Probe-confirmed none, and it needs `--cores >= 2`.

**Does V2 prove what it claims?** Yes for what it tests, no for what the plan implies.
Removing `main.rs:918-924` (not 919–925 — I1) drops four candidates while the writer
still writes four files, so T1 goes red on the paired FASTQ arm. That proves the
whole chain — hook, dump, parse, compare — is live **once, on one arm, at plan time**.
It does not make that permanent: V3 exercises `compare` on synthetic sets and never
touches the hook or the parse. So after V2 is reverted, the only permanent guard that
the hook still fires is the accident in I2. Adding I2's two assertions is what makes
V2's finding durable, and it is cheaper than any alternative.

**V4** is a good test and works on the SE arm it targets. Note it is asserting the
inverse of the rule that C3 shows is false on the specialty arms — so V4 will keep
passing while the invariant it protects is wrong elsewhere.

**V6's "589 + new tests"** — I count 590 `#[test]` attributes across `src/` and
`tests/`. Close enough that the discrepancy is probably a `#[test]` inside a
doc-comment or a `cfg`-gated block; not worth chasing, but do not treat 589 as an
exact gate.

---

## 5. Alternatives (within the settled decisions)

1. **`#[track_caller]` on the pre-flight** — the main recommendation, C2. Gives exact
   site attribution with zero edits to the 11 call sites, which is the whole reason
   the env-var hook was chosen over a threaded flag.

2. **Dump to a file named by the env var, rather than to stderr.** The plan treats
   the variable as a boolean ("present with any value"). Treating the value as a
   path costs one line, removes the "parse a needle out of a noisy stderr stream"
   surface entirely, lets you use NUL separators and so retire A3's newline caveat,
   and means the dump cannot be perturbed by anything else writing to fd 2. Against:
   a gated file write in a released binary is a slightly larger surface than a gated
   `eprintln!`. Either is defensible; worth a deliberate choice rather than an
   inherited one.

3. **Assert the converse with a *knowable* allowance.** The plan rejects `planned ⊆
   created` for good reason (A5) — but the two legitimate over-planning sources are
   enumerable rather than open-ended: `_unpaired_N` under `--retain_unpaired`, and
   the `--demux` per-barcode set. So `planned \ created ⊆ {those two rules}` is
   assertable, and it would catch the mirror-image defect the tripwire currently
   cannot see: a **stale** candidate whose writer was removed, which causes spurious
   collision refusals (#398's shape, from the other side). Against: it is a second
   allowance, therefore a second decay surface, and #398 already showed that
   narrowing a candidate list is a tempting wrong fix. I would offer this to Felix
   as an explicit choice rather than fold it in.

4. **Record the observed write set per arm.** §7 claims "the observed write sets in
   §2.2 become the specification the deferred structural fix needs", but the tripwire
   asserts a subset relation and records nothing. Either print `created` per arm
   under `--nocapture`, or say plainly that the spike files are the record and the
   tripwire is not.

---

## 6. Action items

### Critical

1. **C1 — fix T2 or delete it.** As specified it can never fail: an output that
   overwrites an input is not in `created = after \ before`. Compute `written =
   created ∪ {p ∈ before : hash changed}` and assert T2 over that (it then subsumes
   T3), or drop T2 and promote T3 to "no pre-existing file's bytes changed". Correct
   §3.3's claim that T2 "gives the two orphan arms real coverage" — T3 does.
2. **C2 — put site identity in the dump.** `#[track_caller]` on
   `preflight_output_collisions` + `Location::caller()` in the dump line. Without it
   §3.4's `<N>` needs a hand-maintained arm→site table (which §1 forbids and which
   rots on every `main.rs` edit), and nothing detects a case that silently exercises
   the wrong arm.

### Important

3. **C3 — scope the FastQC non-empty invariant.** `src/specialty.rs` has no `fastqc`
   reference and no CLI conflict rejects the combination, so `--fastqc` is a silent
   no-op on `--hardtrim5/3`, `--clock`, `--implicon` and on orphan arm 1. State that
   in the test file so nobody adds a specialty `--fastqc` case and then "fixes" the
   invariant.
4. **I1 — V2's range is `main.rs:918-924`, not `919-925`.** As written the negative
   control does not compile (918's `if` is left unclosed) and it deletes
   `planned.extend(candidates)` at 925, which empties the whole candidate list
   instead of just the four reports — so V2's stated expectation would not match what
   you see. §10's other five citations are fine as they stand.
5. **I2 — add two non-emptiness assertions**: `planned` non-empty on `Arm::Planned`,
   `created` non-empty on all arms. This is what makes V2's one-off finding permanent.
6. **I3 — drop the in-process unit test.** `env::set_var` is `unsafe` in edition 2024
   and contradicts `ci.yml:55-56`'s nextest-clean invariant. Do V1 in the integration
   file: one arm, two runs, var set / var unset.
7. **I4 — add `--cores 2` variants of cases 13 and 14.** `--cores > 1` switches to
   `parallel.rs`'s independent writer block (`:212-216`); O1's "adds no writer" is
   wrong about it, and the paired parallel writer is currently unobserved.
8. **I5 — require unique `case_root` tags.** Shared PID under `cargo test --release`
   plus `remove_dir_all` means a duplicate tag deletes a sibling's tree mid-run,
   which can convert a T1 miss into a pass.
9. **I6 — say `guarded_inputs`** in §3.2 step 3, so the `--demux` barcode file and the
   `--passthrough` carrier are hashed.

### Optional

10. **O-a — record A8's limit.** `collision_key` folds case (`io.rs:80` →
    `to_ascii_lowercase`), so the tripwire cannot see a candidate that differs from
    the written path only in case — which on Linux is a real defect (writer writes
    `X`, pre-flight guards `x`, an input named `X` is unprotected). The choice is
    right; one line in the module doc stops it being rediscovered.
11. **O-b — replace §2.3's safety argument with the crisp one.** "A fixed point would
    need an infinitely long name" is a length argument. The load-bearing fact is that
    both orphan arms derive the stem with `Path::file_stem()` (`main.rs:1820`,
    `:2424`), which only ever removes a `.`-prefixed extension, while every appended
    suffix begins with `_` (`_val_1`, `_val_2`, `_val.bam`, `_unpaired_N`,
    `_R{1,2}_trimming_report.*`). So `stem + suffix` can never equal `file_name =
    stem + "." + ext` — for any input, with or without `-o`, and case-folding cannot
    change it. Also record that **neither arm reads `cli.basename`**: that is what
    ties the stem to the input, and if a later change honours `--basename` there the
    proof collapses and the follow-up pre-flight issue becomes a live data-loss
    defect rather than defence-in-depth. That sentence belongs in the follow-up issue.
12. **O-c — case 8 needs real gzip bytes.** The helper list has only `write_fastq`.
    A plain file named `.gz` would still exercise the naming branch (`is_gzipped` is
    extension-based, `io.rs:73`; format detection is content-based) but not the gz
    writer. Add `write_fastq_gz`, or take the shortcut deliberately and say so.
13. **O-d — probe `--fastqc_args "--svg"`** (`fastqc.rs:81`) once, to close A2 across
    the curated flag subset. O1's exclusion list does not mention it.
14. **O-g — keep the allowance keyed on the explicit flag, not on the args.** §4's
    `compare(..., fastqc: bool)` is right; §3.3's prose ("a created file is allowed
    when the arm's args contain `--fastqc`") is loose, and `--fastqc_args` *implies*
    `--fastqc` (`cli.rs:333`), so an args-sniffing implementation would misclassify a
    case that passes only `--fastqc_args`. V4 depends on the two being separable, so
    align §3.3's wording with §4's signature.
15. **O-e — add `--retain_unpaired` to case 33.** It is legal on orphan arm 1
    (`main.rs:1837`) and adds two writers, so the arm's full write set gets observed
    for the follow-up pre-flight issue. Case 34 cannot — `--retain_unpaired` is
    rejected with uBAM output (`cli.rs:629`).
16. **O-f — state the one-dump-per-run invariant** (§0) that `parse_planned` relies
    on, so moving a pre-flight inside a per-pair loop is a visible change rather than
    a silent weakening.

---

## 7. Summary

The mechanism is the right one and several of the hardest details are already
correct: the canonicalised-root/`root.join` handling of `lexical_normalise`'s
CWD-relative absolutisation, A4's refusal of a `cfg` gate, the subset direction, and
the choice to key on `collision_key`. A1 and the 11-site inventory hold under
independent tracing, and the case matrix does reach every site.

What needs fixing before implementation: T2 is unfalsifiable as written (C1); the
dump needs to carry its call site, both to implement §3.4 and to stop a case silently
testing the wrong arm (C2); the FastQC non-empty invariant is false on four arms and
so is a decay vector inside the anti-decay mechanism (C3); V2's edit range does not
compile (I1); and the plan's liveness rests on an accident that two one-line
assertions would turn into a guarantee (I2). I3 and I4 are the two places where the
plan's own framing — "no `env::set_var`", "`--cores` adds no writer" — is at odds
with what the tree actually does.
