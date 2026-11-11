# PLAN_REVIEW_A — #414 empirical exact-set tripwire

**Reviewer:** A (independent; no access to Reviewer B's work)
**Target:** `plans/08102026_preflight-tripwire/PLAN.md` (v1, 2026-08-10)
**Tree verified against:** `dev` @ `3008b37` (matches the plan's stated base)
**Verdict:** the mechanism is right and the plan is unusually well evidenced. Three
false-negative paths would let the tripwire pass when it should fail; all three are
cheap to close. Two settled decisions are executed cleanly; one (the V1 unit test)
collides with a documented CI invariant.

---

## 0. Verification log — what I checked against the tree

Everything below was read from `3008b37`, not taken from the plan.

| Plan claim | Result |
|---|---|
| 11 pre-flight sites at `466, 497, 661, 706, 752, 809, 927, 1006, 2047, 2109, 2670` | **Exact.** `grep -n preflight_output_collisions src/main.rs` returns precisely those 11 lines. |
| `io::preflight_output_collisions` at `src/io.rs:148`; `collision_key` `:80`; `pair_output_dir` `:91`; `PlannedOutput` `:137`; `OutputSource` `:104` | **All exact.** (`PlannedOutput` is the `pub type` two lines above `:148`; the plan's `:137` is right.) |
| Site `2670` = `run_specialty_paired`, serving `--clock`, `--implicon`, `--clump_only --paired` | **Correct.** `run_specialty_paired` at `main.rs:2654`; the three closures are at `:526` (clock), `:546` (implicon), `:601` (clump_only paired). |
| Two orphan arms: `main.rs:854` → `run_paired_ubam_single_file` (`:1810`), `run_ubam_output:2006` → `run_ubam_output_paired_single_file` (`:2413`) | **Correct.** Neither function contains a pre-flight call. |
| Hand-rolled `{stem}_R{N}_trimming_report.{txt,json}` at `main.rs:1962`, `:2492` | **Exact**, both. |
| Appendix consumers `919-925`, `1745-1751`, `2040-2045`, `2369-2375`, `1962-1971`, `2492-2501` | **All six exact.** |
| `--demux` summary *is* planned, `demux.rs:156` | **Correct** — `demux_output_paths` pushes `{base}_demultiplexing_summary.txt` at that line. |
| A4's `ci.yml:69` (`cargo nextest run`) and `:76` (`cargo test --release`) | **Exact.** |
| A9: `test_files/ubam_test.bam`, `ubam_paired_test.bam` committed | **Confirmed** via `git ls-files test_files`. |
| Risk 3: "the first `std::env::var` read in `src/`" | **Correct.** Only `build.rs` reads env today (`:18`, `:64`, `:65`). |
| `--extract` not in the `--fastqc_args` curated subset | **Correct** (`fastqc.rs:71-72`) — so no `_fastqc/` **directory** can appear and the two-suffix allowance rule is complete for the supported flag set. |

**A1 (the enumeration) holds — I tried to break it and could not.** Walking every
write-capable entry point reachable from `main()` gives exactly 15 dispatch paths:
`specialty::hardtrim5{,_to_bam}` (466), `hardtrim3{,_to_bam}` (497),
`specialty::clock` / `implicon` / `clump_only::clump_only_paired` (all 2670),
`clump_only_single` (661), `clump_only_paired_to_bam_one_pair` shape B (706) and
shape A (752), `clump_only_single_to_bam` (809), `run_paired` (927),
`run_single_file` (1006), `run_ubam_output_paired_two_files` (2047),
`run_ubam_output_single` (2109), plus the two orphans. `resolve_clump_layout` writes
nothing; `adapter.rs`'s `File::create` calls are test-only; `--clumpify` bins are
in-memory (probe confirms no temp files). **There is no third orphan.**

**The test matrix does reach all 11 sites.** Mapping: 466→18/19, 497→20/21, 661→24/26,
706→28, 752→29, 809→27, 927→11–17, 1006→1–10, 2047→32, 2109→30/31, 2670→22/23/25.
No site is unreached, and the 34-row numbering is internally consistent (rows 18 and
20 each name two tests and consume two numbers).

**T1 would have caught all five historical instances.** #383 (no pre-flight at all →
`planned` empty, `created` non-empty), #388 (paired reports unplanned → case 11), #391
(clump_only paired reports unplanned → case 25), #409 (SE uBAM reports unplanned →
case 30), #398 (reports planned in the mate's dir, written in R1's → case 11,
`paired_trim_split_dirs`, which is exactly the arm the plan chose for V2). That is a
real, checkable regression-coverage argument and it holds.

---

## 1. Logic review

### 1.1 The comparison relation is right

`created \ allowed ⊆ planned`, keyed on `io::collision_key(root.join(p))`, is the
correct construction and the plan's reasoning for it is sound. Two details I checked
that the plan gets right without stating them:

- `root.join(p)` works for **both** relative and absolute dumped paths, because
  `Path::join` on an absolute argument discards the base. So a candidate built from
  an absolute `--output_dir` still compares correctly.
- `collision_key` is `norm_path(lexical_normalise(p))` (`io.rs:80`), so `./x`, `x`,
  and `a/../x` all fold together. The dumped "as constructed" spelling therefore
  cannot disagree with the filesystem-walked spelling on punctuation alone.

### 1.2 The assertion direction and the FastQC allowance are correctly justified

Over-planning really is legitimate on the arms the plan names — the probe's
`pe_retain` case shows both `_unpaired_N` files created even with nothing stranded,
and `demux` plans one file per barcode. `planned ⊆ created` would be wrong. Agreed.

The FastQC rule is complete for what FastQC can actually emit: `--extract` is not in
the curated flag subset, so there is no directory case, and `fastqc::run` writes into
`--output_dir` or (unset) beside the file it is given — the probe's `fqc_odir` and
`ubam_fqc` rows confirm both branches, including a `.bam` primary and the
`--passthrough` carrier.

### 1.3 §2.3's non-aliasing proof is sound — and I found what it silently rests on

I tried to construct a counterexample. Every namer on both orphan arms builds
`dir.join(format!("{stem}_…", stem = input.file_stem()))` (`main.rs:1819-1830`,
`:1962-1971`, `:2422-2432`, `:2492-2501`). For an output to equal its input you would
need `X.ext == X_val_1.fq` for a shared prefix `X`, i.e. `.ext == _val_1.fq` — a `.`
must equal a `_`. Impossible. `--output_dir` does not help (it moves both sides or
neither), and with one input there is no output-vs-output pair. **The claim is
correct.**

But it is correct *only because `--basename` is not threaded into either arm* — I
verified there is no `basename` reference anywhere in `main.rs:1810-2004` or
`:2413-2560`. `--basename` decouples the output stem from the input filename, which
is precisely the premise of the proof: with basename honoured, `--basename ip` plus an
input named `ip_val_1.fq` produces `ip_val_1.fq` and the input is destroyed. The
plan's §2.3 argument therefore has a documented-invariant dependency it does not
name, and `--basename` on these arms is today a **silent no-op** (see O4).

### 1.4 T2 is structurally near-vacuous — the plan overstates it

> **T2 — no created path aliases an input** (all 15 arms) … "This is the assertion
> that gives the two orphan arms real coverage."

`created` is defined in §3.2 step 6 as "files under the root after the run, minus
`before`". An output that overwrites an input writes to a path that **already existed
in `before`**, so it is excluded from `created` by construction. T2 can never see the
#409 harm. On macOS CI (case-insensitive APFS) even a case-variant alias resolves to
the same inode and is likewise absent from `created`. T2 only has teeth on a
case-sensitive filesystem, for a fold-equal-but-not-identical path — i.e. on the Linux
runner, for one narrow shape.

T3 is doing all the real work on the orphan arms. That is fine — but the plan should
say so, because "T2 gives the orphan arms real coverage" is the sentence a future
maintainer will rely on when deciding those arms are covered. Either redefine T2 over
`after` (`for every input, no *other* path in after folds to its key`) or downgrade the
prose. See I5.

### 1.5 §2.3's arm-14 write set is incomplete

The table lists arm 14 as writing `ip_val_{1,2}.fq` plus the four reports.
`run_paired_ubam_single_file` also writes `{stem}_unpaired_{1,2}{ext}` when
`--retain_unpaired` is set (`main.rs:1837-1839`), and `--retain_unpaired` is rejected
only for uBAM *output*, so `--paired --retain_unpaired ip.bam` is a legal FASTQ-out
invocation. Two more writers on an arm with no pre-flight at all. Because arm 14 is
`Arm::Unplanned` this does not cause a false T1 pass — but the follow-up issue "add a
pre-flight to the orphan arms" would be planned from a list that is missing two paths,
and T3 never hashes against them. See I3.

---

## 2. Assumptions

A3, A5, A6, A7, A8, A9 are sound and I have nothing to add. A1 and A4 I verified
above and both hold, A4 with exactly the CI line numbers it cites. Three comments:

**A2 — "`--fastqc` is the only source of unplanned writes. Established by observation
across 40 combinations."** True for what was probed, but the 40 combinations never
crossed `--fastqc` with arms 14, 15, or any specialty arm. I closed that gap by
reading instead: `fastqc::run` is called from exactly six places — `main.rs:1459`
(`run_single_file`), `:1778`/`:1784`/`:1796` (`run_paired`, incl. the carrier),
`:2242` (`run_ubam_output_single`), `:2399` (`run_ubam_output_paired_two_files`),
`:2522` (`run_ubam_output_paired_single_file`) — plus the `cli.fastqc` pass-through
into the `clump_only` drivers (`:636`, `:671`, `:715`, `:778`, `:819`). It is **not**
called from `specialty::hardtrim*`, `clock`, `implicon`, or
`run_paired_ubam_single_file`. So A2 survives, but its corollary bites the plan's own
self-policing rule (I2), and `--fastqc` being a silent no-op on five arms is a real
product observation worth a separate issue.

**A2's gate is a disjunction, and the plan's rule is not.** Production runs FastQC on
`if cli.fastqc || cli.fastqc_args.is_some()` (`main.rs:1458`, `:1777`, `:2241`,
`:2398`, `:2521`). §3.3 says the allowance applies "when the arm's args contain
`--fastqc`". A case passing only `--fastqc_args` would produce allowed files while the
rule declares the allowed set must be empty → spurious red → the exact weakening
pressure §3.4 exists to resist. See I1.

**Unstated assumption: every case's run does real work.** Nowhere does the plan
require `created` to be non-empty, and `main.rs:219-226` routes
`DisplayHelpOnMissingArgumentOrSubcommand` to **stdout with exit 0**. This is the
plan's biggest hole; see C1.

**Unstated assumption: `inputs` for T2/T3 equals `guarded_inputs`.** `guarded_inputs`
(`main.rs:74-83`) is `cli.input` **plus** the `--passthrough` file **plus** the
`--demux` barcode file. If cases 5, 6, 14, 16 hash only the FASTQ inputs, the two
files the pre-flight was extended to guard are the ones T3 does not watch. See I7.

---

## 3. Efficiency

No concerns. ~34 subprocess launches × 2 profiles (nextest debug + `cargo test
--release`) ≈ 68 invocations on 40-read fixtures; the existing
`integration_output_collision.rs` already runs 57 such tests, so the precedent is
measured, not hoped for. `snapshot` walks ≤ ~20 files; `BTreeSet<String>` comparison
is trivially sized. The hook is one `var_os` per run.

The hook's cost claim in §6 ("one `var_os` lookup per run (11 sites, but one call per
run)") is accurate — every arm returns after its single pre-flight.

One robustness-not-speed note: under `cargo test --release` all 34 tests share one
process, so `tempdir`'s `std::process::id()` discriminator
(`integration_output_collision.rs:31-36`) is identical for all of them and uniqueness
rests entirely on the tag string. A duplicated tag means one test `remove_dir_all`s
another's tree mid-run — a flake that would present as a random unplanned-path
failure. See O3.

---

## 4. Validation sufficiency — can this tripwire pass when it should fail?

Three ways, in descending severity.

### 4.1 A case can exit 0 having done nothing (Critical)

`src/main.rs:215-227`:

```rust
if matches!(err.kind(), DisplayHelp | DisplayHelpOnMissingArgumentOrSubcommand | DisplayVersion) {
    print!("{err}");
    std::process::exit(0);
}
```

A case whose args are malformed such that clap reports a *missing argument* — the
positional input list swallowed by a preceding flag, a quoting slip on
`--fastqc_args "-o qc"`, a renamed flag after a future CLI change — prints help to
**stdout** and exits **0**. Then:

- §3.2 step 5's "require exit 0" is satisfied.
- `created` is empty → `created \ allowed ⊆ planned` is **vacuously true**.
- `planned` is empty → T1 still passes.
- T2 has nothing to iterate → passes.
- T3's hashes are unchanged → passes.

All three assertions pass on a case that never reached its dispatch arm. This is not
hypothetical: three of the plan's own probe cases died on fixture bugs, and the plan
explicitly inherits those as case requirements — a fourth one that lands on the
help path instead of an error path would be invisible.

The same vacuity arises from any run that writes entirely outside the case root
(absolute `--output_dir`, or a namer change that escapes it). `snapshot` is
root-scoped and cannot see it.

**Fix (5 lines):** for every case assert `!created.is_empty()`; for every
`Arm::Planned` case assert `!planned.is_empty()`. That single pair of assertions makes
every `⊆` in the suite non-vacuous by construction and kills both variants.

### 4.2 There is no permanent end-to-end negative control (Critical)

V2 proves the chain can go red, but it is manual, one-off, and performed on a scratch
worktree that is then reverted — nothing in CI preserves it. V3 (`compare` called with
synthetic sets) makes redness permanent for the **pure function only**; it does not
exercise hook → stderr → `parse_planned` → `compare`. So the wiring that V2 validated
by hand is unguarded from the moment V2 is reverted, which is precisely the
false-negative class the repo has been burned by (`feedback_harness_false_negatives`:
"prove the check can fail before trusting that it passed").

**Fix (≈5 lines, and it subsumes half of V1):** one extra `#[test]` that runs a known
arm with `DUMP_PLANNED_ENV` **unset**, asserts `parse_planned(stderr)` is empty, and
asserts `compare(root, created, /*planned*/ &[], inputs, false)` returns `Err`. That
drives the whole chain and proves it goes red, permanently, in both CI profiles. It
also gives V1's "inert when unset" half an implementable home.

Second, weaker gap: V2 proves redness for **one** site (927). The other 12 planned
arms' wiring is unproven, and the fix in 4.1 (`planned` non-empty per case) is what
generalises it.

### 4.3 The allowance has two ways to widen without anyone weakening it (Important)

§4's signature takes `fastqc: bool` as a parameter *separate from* the args, and §3.3
defines the rule in terms of *the args*. Two sources of truth for one predicate:

- case passes `--fastqc`, sets `fastqc: false` → spurious red → weakening pressure;
- case omits `--fastqc`, sets `fastqc: true` → **allowance silently active on an arm
  that produced no FastQC output**, and the "must be non-empty" guard fires… which is
  a red, so this direction self-corrects. But invert it once more: a case that passes
  `--fastqc_args` only, with `fastqc: true`, is *correct* under production semantics
  and *incorrect* under §3.3's prose.

And the non-empty half fires spuriously on the five arms where FastQC never runs
(§2 above): `--fastqc` on `hardtrim5/3`, `clock`, `implicon`, or the FASTQ-out
interleaved orphan produces zero allowed files, so "the allowance is masking a broken
case" would be asserted about a case that is fine. A red that a reader can prove is
wrong is the strongest possible argument for deleting the guard.

**Fix:** derive the predicate once, inside `compare`, from the args slice, as
`args.contains("--fastqc") || args.iter().any(|a| a == "--fastqc_args")` — mirroring
`main.rs:1458` — and scope the non-empty assertion to a named `FASTQC_CAPABLE_ARMS`
set. V4's negative control then becomes "created `_fastqc.zip` with no FastQC flag in
args must fail", which tests the same property with one source of truth.

### 4.4 Nothing ties a case to the site it is meant to exercise (Important)

The plan's central claim is "all 11 sites covered". Nothing in the test verifies which
pre-flight ran. Today the mapping holds (I checked all 11). But dispatch conditions in
`main()` are dense and change often — the `--clump_only --paired` N=1 bail at
`main.rs:598-604` was added exactly because a shape silently routed to the wrong arm.
A future condition change can quietly funnel two cases onto one site, leaving a site
unexercised while the suite stays green and the module doc keeps claiming coverage.

This is fixable **without touching any call site**, which keeps it inside settled
decision 3: add `#[track_caller]` to `preflight_output_collisions` and dump
`std::panic::Location::caller()` alongside the candidates. Then one extra `#[test]`
runs a representative arm per site, collects the dumped locations, and asserts **11
distinct values** — no hard-coded line numbers, so it is not brittle to unrelated
edits in `main.rs`. `#[track_caller]` embeds only the relative `src/main.rs` spelling
(the same strings every `unwrap()` already embeds), so the reproducible-build job's
"no absolute paths in the binary" rule is untouched.

### 4.5 Things I checked that are **not** false-negative paths

Stated plainly so the list above is not padded:

- **Interleaving on stderr.** The dump is at the top of the pre-flight, before any
  worker thread exists, and `eprintln!` locks. No split lines.
- **Prefix collision.** `TRIM_GALORE_PLANNED\t` cannot be produced by any other
  stderr writer, and FastQC's chatter is unaffected.
- **Env leakage between cases.** `run_arm` sets the variable on the child `Command`,
  not the test process, so nextest and threaded `cargo test` are both safe (with the
  one exception in C2).
- **Non-zero exit swallowing.** §3.2 step 5 makes it a hard failure, so a case that
  cannot run fails. Good — that is what closes most "case didn't exercise its arm"
  variants, and it is why 4.1 (exit **0** with no work) is the residual.
- **`Arm::Unplanned` asserting `planned` empty.** Correct, and it will genuinely go
  red the day either orphan gains a pre-flight. Nice design.
- **The dump firing before the collision check.** Correct placement: a run that bails
  still dumps, so the record exists for debugging even though the case then fails on
  exit status.

---

## 5. Alternatives (within the settled decisions)

1. **`#[track_caller]` site dump** — §4.4. Best value in this list: it turns the "11
   sites covered" claim from prose into a machine-checked assertion, costs one
   attribute plus one line in the dump, and edits no call site.
2. **Hand-enumerated exact set for the two orphan arms only.** The plan rejects
   hand-enumeration on principle ("it rots the same way"), and for the 13 pre-flighted
   arms that is right. But for arms 14/15 there *is* no expectation available from the
   binary — `Arm::Unplanned` skips T1 and T2 is near-vacuous (§1.4), so their write set
   is exercised but not checked. `assert_dir_holds_only`, which
   `integration_output_collision.rs` already provides, gives them T1-equivalent
   strength for two arms and 8 hand-written paths. The rot argument is weaker here
   because the follow-up issue will replace it with a real pre-flight.
3. **One "refused run creates nothing" case.** The tripwire tests the pre-flight's
   *coverage*; nothing tests its *purpose*. One case whose args provoke a genuine
   collision, asserting non-zero exit **and** `created.is_empty()`, pins the
   fail-before-any-write contract and gives O3's partial-output question a concrete
   anchor. Cheap, and it is the assertion #409 would most directly have wanted.
4. **Assert every planned path resolves inside the case root.** One line; closes the
   "writer escaped the snapshot" variant of §4.1 permanently rather than by convention.
5. **Split the hook for testability** — see C2. `fn dump_planned_lines(&[PlannedOutput])
   -> Vec<String>` (pure, unit-tested) plus an env-gated printer is strictly better
   than the plan's step-2 unit test and removes the `env::set_var` problem entirely.

---

## 6. Action items

### Critical

- **C1 — Assert `created` non-empty for every case, and `planned` non-empty for every
  `Arm::Planned` case.** `main.rs:219-226` exits **0** on
  `DisplayHelpOnMissingArgumentOrSubcommand`, so a malformed case prints help, creates
  nothing, and passes T1/T2/T3 vacuously. Same vacuity if a run writes outside the
  root. This is the one defect that can make the whole suite green while covering
  nothing. (§4.1)

- **C2 — Drop the in-process env-var unit test; split the hook instead.** Outline step
  2 (`dump_env_emits_one_line_per_candidate`, "with the var set…") requires
  `std::env::set_var`, which (a) is `unsafe` under edition 2024, (b) leaks to every
  concurrently-running test in the same binary under `cargo test --release`
  (`ci.yml:76`), and (c) directly contradicts the invariant recorded in
  `.github/workflows/ci.yml:55-56`: *"no `static mut` / `OnceCell` / `env::set_var`
  globals exist — the codebase is nextest-clean."* I confirmed `set_var` appears
  nowhere in `src/`, `tests/`, or `build.rs` today. Replace with a pure
  `dump_planned_lines(&[PlannedOutput]) -> Vec<String>` (unit-tested, no env) plus an
  env-gated printer; test the gate end-to-end per C3.

- **C3 — Add a permanent end-to-end negative control.** V2 is manual and reverted; V3
  covers only the pure `compare`. Add one `#[test]` that runs a known arm with
  `DUMP_PLANNED_ENV` unset, asserts zero prefixed lines, and asserts `compare` with an
  empty `planned` returns `Err`. Drives hook → stderr → parse → compare and proves the
  chain can go red, in both CI profiles, forever. (§4.2)

### Important

- **I1 — One source of truth for the FastQC predicate.** Derive it inside `compare`
  from the args, as `--fastqc || --fastqc_args`, matching production's gate
  (`main.rs:1458`, `:1777`, `:2241`, `:2398`, `:2521`). §4's `fastqc: bool` parameter
  plus §3.3's args-based prose is two sources of truth for the one predicate whose
  drift the plan itself identifies as fatal. (§4.3)

- **I2 — Scope the "allowed set must be non-empty" rule to FastQC-capable arms.**
  FastQC is never invoked from `specialty::hardtrim*`, `clock`, `implicon`, or
  `run_paired_ubam_single_file` (verified: no `fastqc` reference in
  `main.rs:1810-2004`). `--fastqc` there is a silent no-op, so the guard would go red
  on a legitimate case — a provably-wrong red is the strongest argument for deleting
  the guard. Name the capable arm set explicitly. Worth a separate product issue that
  `--fastqc` is silently ignored on five arms.

- **I3 — §2.3's arm-14 write set is missing `{stem}_unpaired_{1,2}`**
  (`main.rs:1837-1839`, legal under `--paired --retain_unpaired ip.bam` since
  `--retain_unpaired` is refused only for uBAM output). Add
  `orphan_paired_fastq_retain_unpaired`, and correct the §2.3 table — the follow-up
  pre-flight issue will be planned from it.

- **I4 — Add a paired `--output_dir` case.** Case 2's note claims it exercises "the
  `-o` branch of every namer"; it exercises the SE namers only. `pair_output_dir`'s
  `output_dir.map(...)` arm (`io.rs:91`) and the `Some(&pair_dir)` threading at
  `main.rs:921-922` versus the writer at `:1745-1751` can diverge **only** in that
  branch, and that is the exact shape of #398. No case in the matrix covers it.

- **I5 — Correct T2's stated strength.** `created = after \ before` excludes any path
  that already existed, so T2 cannot detect an input being overwritten — the #409
  harm — and on case-insensitive APFS it cannot detect a case alias either. Either
  redefine it over `after`, or stop describing it as "the assertion that gives the two
  orphan arms real coverage" and say plainly that T3 carries them. (§1.4)

- **I6 — Machine-check the site coverage** via `#[track_caller]` +
  `Location::caller()` in the dump, plus one test asserting 11 distinct sites. Without
  it, "all 11 sites covered" is prose that decays silently on the next dispatch-
  condition change. (§4.4)

- **I7 — `inputs` must mirror `guarded_inputs`** (`main.rs:74-83`): `cli.input` plus
  the `--passthrough` file plus the `--demux` barcode file. Cases 5, 6, 14 and 16
  otherwise skip T3 on exactly the two files the pre-flight was extended to guard.

### Optional

- **O1 — Case 4 must pre-create `qc/`.** `ensure_output_dir` only creates
  `--output_dir`; the probe did the `mkdir` itself
  (`spikes/probe_flag_variants.sh:58`). Without it the run errors and the case fails
  for the wrong reason.
- **O2 — Case 9 should assert the `clumpify:` notice on stderr.**
  `resolve_clump_layout` (`main.rs:172-204`) falls back to plain mode with a *warning*
  and exit 0 when `--memory` is below the floor. At `--cores 2` the floor is ≈536 MiB
  against the 1 G default so it passes today, but a cores bump or a change to
  `STATIC_OVERHEAD_BYTES` / `MIN_BIN_BYTES` silently turns "confirms no temp-file
  writers" into "confirms nothing".
- **O3 — Require unique case-root tags.** `tempdir` keys on `std::process::id()`
  (`integration_output_collision.rs:31-36`), identical for all 34 tests under `cargo
  test --release`; a duplicate tag means one test `remove_dir_all`s another's tree
  mid-run.
- **O4 — Record that §2.3's non-aliasing proof depends on `--basename` not being
  threaded into either orphan arm** (verified absent from `main.rs:1810-2004` and
  `:2413-2560`). With basename honoured, `--basename ip` plus an input named
  `ip_val_1.fq` aliases and the proof lapses. Put the dependency in the
  `Arm::Unplanned` doc comment where a future implementer will read it.
- **O5 — Add one "refused run creates nothing" case** (non-zero exit **and** empty
  `created`). Tests the pre-flight's purpose rather than its coverage; gives O3's
  partial-output question an anchor.
- **O6 — Assert every planned path resolves inside the case root.** One line, closes
  the escaped-writer blind spot permanently.
- **O7 — One `--no_report_file` case on a paired or `clump_only` arm.** Case 7 covers
  the SE gate only; `clump_report_candidates`' `no_report_file` early return
  (`main.rs:128-130`) and the writers' own gates are the other half.

---

## 7. Where the plan is simply right

Said plainly, so the action list is not read as a verdict on the whole document.

- Every line reference in the plan checks out against `3008b37` — all 11 sites, all
  six appendix consumers, `demux.rs:156`, both `ci.yml` refs, both `io.rs` anchors.
  That is unusual and it made this review much faster.
- A1's enumeration is exhaustive; I looked for a third orphan and there isn't one.
- The arm→site mapping, including 2670 serving three modes, is correct.
- The assertion direction, the `collision_key` comparison relation, and the decision
  to compare sets rather than sequences are all correctly reasoned, and A4's
  no-`cfg`-gate argument is exactly right for the reason it gives.
- T1 as specified would have caught all five prior instances — I checked each.
- The FastQC allowance is a rule rather than a list, and it is complete for the
  supported `--fastqc_args` subset (no `--extract`, so no directory case).
- §3.4's insistence that the failure message forbids the wrong fix is well judged, and
  §12's "how this tripwire could die" is the right thing to have written down.
- The evidence discipline — rebuilding from HEAD after finding `target/` stamped
  `38deb1e`, and converting three of its own probe failures into stated case
  requirements — is the reason the residual findings above are as narrow as they are.
