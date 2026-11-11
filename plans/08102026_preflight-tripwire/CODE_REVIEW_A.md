# Code Review A — #414 pre-flight tripwire

**Branch:** `feature/414-preflight-tripwire`
**Reviewer:** A (independent; no access to Reviewer B's work)
**Date:** 2026-08-10

Files under review:
- `src/io.rs` — modified, +68 lines (2 consts, `dump_planned_lines`, `#[track_caller]`, 7-line gated dump, 1 unit test)
- `tests/integration_preflight_tripwire.rs` — new, untracked; 1263 lines / 50 tests as reviewed, 1284 / 51 after my fixes

**Verdict: sound, and the mechanism works.** Three low-risk fixes applied (§4); four
recommendations, none blocking (§6).

---

## 0. The tree I reviewed, and an unplanned negative control

`git diff --stat dev` showed **`src/main.rs | 1 +`** alongside the expected `src/io.rs`:

```rust
// src/main.rs:1014
std::fs::write(input.with_extension("MUTANT"), b"x")?; // REVIEW-B MUTATION
```

A live probe belonging to the other reviewer, not part of the change. My first run of the
suite in the shared tree failed 11 tests because of it. It was reverted part-way through my
review, so the tree is clean again (`git diff --stat dev` → `src/io.rs` only), but the
window is why I did **all** verification in an isolated `git worktree` pinned at `dev` @
`3008b37`, with `git diff dev -- src/io.rs` applied and the test file copied in
(`grep -c MUTANT src/main.rs` → 0). Results below are from that worktree.

The accident is worth recording, because it is the strongest single piece of evidence in
this review: **an unplanned new writer, added by someone who was not coordinating with the
test, turned 11 cases red with exactly the intended message.**

```
created but never planned: ["a.MUTANT", "b.MUTANT"]
Add each to the candidate list at src/main.rs:1006. Do NOT add them to the FastQC
allowance and do NOT narrow the candidate list — those are how this tripwire dies.
The 6 candidates that arm did plan: ["a.fastq_trimming_report.json", ...]
```

That is the whole chain — hook → stderr → parse → compare → message — firing on a real
mutation, under conditions nobody arranged.

---

## 1. Summary

The hook is 7 lines in the one function all 11 call sites funnel through, touches no
candidate-assembly code, and is inert without the environment variable. `Location::file()`
yields the source-relative `src/main.rs:NNNN` spelling — observed directly in dump output —
so neither the reproducible-build nor the `--version` provenance invariant is at risk. 50/50
tests pass in debug and release; 640 tests total; clippy `-D warnings` and `cargo fmt` clean.

On the question that matters — *can this pass when it should fail?* — I found no path that
makes T1 vacuous, and I confirmed the harder half empirically rather than by reading:

- **Every one of 21 probed case shapes reaches exactly the pre-flight site the plan's table
  claims**, exactly one site per run, 11 distinct sites, both orphan arms dumping nothing (§2).
- **T3 fires.** I planted a clobber-the-input mutation and confirmed the T3 assertion goes
  red with the right message while T1 stays green — which is also a live demonstration of why
  the plan was right to delete v1's T2 (§5).
- **The marker system is fail-safe in both directions.** Mis-marking a planned arm
  `.unplanned()`, or a capable arm `.fastqc_incapable()`, or leaving a capable arm marked
  capable after FastQC silently stops running, each produces a red rather than a green. That
  is a genuinely good property and it was not obvious it would hold.

The real weaknesses are not in T1:

| | Finding | Priority |
|---|---|---|
| H1 | T3 has no *permanent* self-test, contrary to the file's own standard — and it is the sole substantive assertion on the three orphan cases | High |
| M1 | T5 cannot detect the scenario its comment says it guards; A10's "machine-checked" is not true as implemented | Medium |
| M2 | No case asserts which arm it reached, which the plan said the site field was for | Medium |
| L1–L5 | Allowance breadth, cosmetics, unprobed A2 residual | Low |

Two of these (H2, L2 in my original numbering) I fixed directly; see §4.

None is a live false negative. All are decay paths inside a test whose entire purpose is to
resist decay — the class the plan itself calls fatal — which is why I am naming them rather
than filing them as polish.

---

## 2. Where the code is sound

Stated plainly rather than manufacturing findings against it.

**The hook.** Placement (`src/io.rs:174-181`) is correct: before `input_keys`, so the dump
fires on a run that goes on to `bail!` — which `a_refused_run_writes_nothing` depends on.
Inert when unset: one `std::env::var_os` per run. `#[track_caller]` is sound and
semver-compatible on a `pub fn`; it embeds `Location` statics, but the same strings are
already in the binary from every `unwrap()`.

**Reproducibility and provenance, checked against a release binary:**

| Check | Result |
|---|---|
| dump lines with the variable unset | 0 |
| `--version` long form em-dash provenance line | present |
| `-V` | `trim_galore 2.3.0`, no provenance line |
| absolute paths referencing the repo/worktree source | **0** |
| absolute paths total | 424, **all** `~/.cargo/registry/...` (pre-existing dependency paths) |
| relative `src/main.rs` `Location` string | present |
| new strings introduced | exactly the two hook constants |

**No candidate-assembly code touched.** `git diff --stat` confirms `src/io.rs` only, so
candidate order — which decides which path a collision message names, and which existing
tests depend on — is preserved. This was the plan's central constraint and it was respected.

**`dump_planned_lines`** is pure, one line per candidate, site repeated per line so
`parse_planned` needs no state. Its unit test pins the wire format, the tab count, and that
an empty list dumps nothing. Avoiding `env::set_var` was right: `unsafe` under edition 2024,
and it would contradict `ci.yml`'s nextest-clean invariant.

**The FastQC predicate is correct and the contested claim holds.** I verified `cli.rs:337` is
`#[clap(long = "fastqc_args", allow_hyphen_values = true)]` with **no** `requires`, and that
all five production gates are `cli.fastqc || cli.fastqc_args.is_some()` (`main.rs:1458`,
`:1777`, `:2241`, `:2398`, `:2521`). So `fastqc_active` (`:211`) mirrors production, and the
plan was right to decline Reviewer B's O-g. Better than "verified in the plan":
`se_trim_fastqc_args_only` *proves* it in the suite, because the self-policing rule
`(true, true) if allowed.is_empty()` would fail that case if FastQC had not actually run.

**The `fastqc_incapable` markers are right, including the pair that is easy to conflate.**
No `cli.fastqc` reference exists in `run_paired_ubam_single_file` (`main.rs:1810-2003`), so
orphan arm 14 is correctly incapable — but `run_ubam_output_paired_single_file`
(`:2413-2558`) **does** gate at `:2521`, so orphan arm 15 is correctly left *capable*.
Marking both orphans alike was the natural mistake and the code avoids it. Nothing before
`:636` reads `cli.fastqc`, confirming hardtrim/clock/implicon are no-ops (#421).

**Case→site binding, verified independently.** 21 shapes probed against a clean debug binary:

| case shape | site | case shape | site |
|---|---|---|---|
| se_trim_two_inputs | 1006 | clock_pair | 2670 |
| se_trim_from_ubam_input | 1006 | implicon_pair | 2670 |
| se_trim_demux | 1006 | clump_only_paired_split_dirs | 2670 |
| se_trim_clumpify | 1006 | clump_only_se | 661 |
| paired_trim_split_dirs | 927 | clump_only_ubam_se | 809 |
| paired_trim_passthrough_cores | 927 | clump_only_ubam_paired_interleaved | 706 |
| hardtrim5_{fastq,ubam} | 466 | clump_only_ubam_paired_two_files | 752 |
| hardtrim3_{fastq,ubam} | 497 | ubam_out_se | 2109 |
| orphan_paired_fastq | *(none)* | ubam_out_paired_two_files | 2047 |
| orphan_paired_ubam | *(none)* | | |

Row-for-row match with the plan's §2.4 table; one site per run; 11 distinct.
`grep -n preflight_output_collisions src/main.rs` returns exactly those 11 lines.

**A10 holds today.** `main.rs:927` sits after the candidate loop closes at `:925`, outside
`cli.input.chunks(2)`; per-pair processing starts at `:939`. (See M1 for why that is *not*
machine-checked.)

**Tag uniqueness.** All 51 `case_root` tags are distinct — checked mechanically, not by eye.
This matters because `case_root` opens with `remove_dir_all` and every case shares a PID
under `cargo test --release`.

**`parse_planned` fails safe.** The A3 newline caveat cannot produce a false *pass*: a
newline truncates the dumped path, making the planned key wrong and driving T1 **red**; the
dropped continuation line does the same. Extra tabs in a path survive intact, because
`split_once` splits once and keeps the remainder. `eprintln!` takes the `Stderr` lock per
invocation and the dump runs on the main thread before any worker spawns, so lines cannot
interleave.

**`snapshot` (`:133`) does not miss writes under the root.** It recurses into directories
created during the run (pushed onto the stack), records symlinks as files rather than
following them — so no loop can hang it and a symlinked output still registers — and excludes
directories, which is what lets `ensure_output_dir`'s pre-created `out/` and the hand-created
`qc/` not count as writes, and what makes `a_refused_run_writes_nothing`'s
`created.is_empty()` meaningful.

**T4 is not decoration.** `created.is_empty()` genuinely closes the
`DisplayHelpOnMissingArgumentOrSubcommand` → stdout/exit-0 hole; `comparison_rejects_an_empty_run`
pins it. `hash_file` (`:152`) is content-sensitive and deterministic within a process
(`DefaultHasher::new()` is fixed-key SipHash-1-3; `Hash for Vec<u8>` length-prefixes then
writes the bytes).

**The marker system is symmetric.** Mis-marking fails red in every direction:
`.unplanned()` on a planned arm → "promote the case to `Arm::Planned`"; `.fastqc_incapable()`
on a capable arm with a flag → "drop the marker"; capable+active with no artefact → "the
allowance is masking a broken case". Three of the four ways a maintainer could quietly
mis-classify a new case are caught by construction.

---

## 3. Issues

### H1 — High: T3 has no permanent self-test, and it is the whole of three cases' coverage

`tests/integration_preflight_tripwire.rs:287-297`:

```rust
// T3 — no pre-existing file's bytes changed. This is #409's harm asserted
// directly, and it is what covers the `Arm::Unplanned` arms.
for (path, before_hash) in &hashes {
    assert_eq!(hash_file(path), *before_hash, ...);
}
```

**The mechanism works — I proved it.** Planting `std::fs::write(input, b"CLOBBERED")?` after
`run_single_file` (worktree `main.rs:1014`) makes `se_trim_two_inputs` fail with
`input …/a.fastq was modified by the run`, while T1 stays green. That also demonstrates
in passing why deleting v1's T2 was right: the clobbered path is in `before`, so it never
enters `created`, and no `created`-based assertion could ever see it.

The finding is that this proof is a **throwaway probe that gets reverted** — precisely the
arrangement the plan rejected for the rest of the file. §5 step 5 requires "five self-tests,
**all permanent**", and §13 deviation 1 claims "Every assertion inside it now has its own
proof it can fire". T3 escapes that claim on a technicality: it is not inside `compare`.

Every other assertion has a permanent self-test — T1
(`comparison_reports_an_unplanned_path`), T4 twice, T5, the allowance, over-planning. T3
has none, and across all 50 tests the loop never observes a changed hash, because no case
stages an input that any arm overwrites. So a future regression that empties `hashes` or
builds it against the wrong root would pass everywhere, silently.

This matters more than a missing test usually would: §3.3 says T3 "is what gives the two
orphan arms their coverage", so for `orphan_paired_fastq_from_interleaved_bam`,
`orphan_paired_fastq_retain_unpaired` and `orphan_paired_ubam_from_interleaved_bam` — where
T1 and T5 are skipped by construction — T3 is the only substantive assertion, and it is the
one assertion with no standing proof it can fail.

**Root cause is an unrecorded deviation.** The plan's §4 signature is
`compare(root, created, planned, guarded_inputs, args, arm)`. The implementation dropped
`guarded_inputs` from `compare` and moved T3 into `Tripwire::run`, which is what put it
outside the self-test net. §13's "Deviations from §5, all additive" does not mention it.

**Recommendation.** Extract the loop as a pure checker and self-test it — ~15 lines, and it
restores the property the rest of the file already has:

```rust
/// T3 — no pre-existing file's bytes changed. Pure, so `comparison_detects_a_clobbered_input`
/// can prove it fires without a planted mutation in `main.rs`.
fn unchanged(before: &BTreeMap<PathBuf, u64>) -> Result<(), String> {
    for (path, h) in before {
        if hash_file(path) != *h {
            return Err(format!("input {} was modified by the run", path.display()));
        }
    }
    Ok(())
}
```

### M1 — Medium: T5 cannot detect the scenario its own comment says it guards

`:391-401`:

```rust
// T5 — one arm per run. Every pre-flight call site sits outside its per-pair loop,
// so a run dumps exactly one arm's list. Moving one inside a loop would union two
// arms' lists and weaken T1 silently; this makes that a visible failure.
let site = &planned[0].0;
if let Some((other, path)) = planned.iter().find(|(s, _)| s != site) {
```

`find(|(s, _)| s != site)` detects two *distinct* call sites. Moving an existing pre-flight
inside its per-pair loop produces N invocations from the **same** `file:line`, so `find`
returns `None` and T5 passes — the exact scenario the comment names.

I checked whether the dump can distinguish the two shapes at all, and it cannot: for
`paired_trim_multipair` (2 pairs), one call carrying 12 candidates and two calls carrying 6
each emit byte-identical stderr. There is no invocation boundary in the format. So A10
("Exactly one pre-flight invocation per run … **T5 asserts it**"; table: "Fixed,
machine-checked") is not machine-checked.

Practical impact is limited — the union of per-pair lists equals the single list, so T1's
strength is unchanged for any path present in its own pair's list. The defect is that the
assertion, its comment, and the plan all claim a property the code does not deliver, inside
the one file whose value depends on its claims being exact.

**Recommendation**, in preference order:

1. Give the dump an invocation boundary and assert exactly one — one
   `TRIM_GALORE_PLANNED_BEGIN\t<site>` line per call, or an invocation counter appended to
   the site field. ~4 lines in `io.rs`, and it makes T5 true as advertised.
2. If that is more hook than the property is worth, correct the comment **and** A10 to what
   T5 actually asserts: *no run reaches two distinct pre-flight sites*. Do not leave the
   stronger claim standing.

I did not apply either: choosing between them is a design call, and silently rewriting the
comment would pre-empt option 1.

### M2 — Medium: no case asserts which arm it reached

The plan's §3.1 says the site field is "what binds a case to the arm it is meant to
exercise", and §5's table assigns a site to every case. Neither survived into the code:
`Tripwire` has no expected-site field, and `compare` uses `planned[0].0` only for T5 and the
failure text. A case whose args silently route to a different arm still passes — T1 holds
against whatever arm it hit — while the suite's claim to cover 13 pre-flighted arms quietly
becomes false.

Not hypothetical in kind: `main.rs:597`'s explicit `--clump_only --paired` N=1 bail exists
precisely because a shape routed to the wrong arm.

`every_site_is_reached_by_some_case` is a partial guard — 11 sites, required distinct — but
it cannot separate the three arms sharing site 2670 (`--clock`, `--implicon`,
`--clump_only --paired`), and only `--clock` is in its shapes list. If
`clump_only_paired_split_dirs` ever routed to 927, nothing would go red.

I verified by hand that all 21 shapes route correctly today (§2), so this is decay risk, not
a present defect.

**Recommendation.** Cheapest useful version: extend `every_site_is_reached_by_some_case`'s
shapes with `--implicon` and `--clump_only --paired` and assert those two *do* collide with
`--clock` on 2670 — that pins the three-way sharing as an asserted fact rather than prose,
without hard-coding a line number. A fuller version adds `expect_site: Option<&str>` to
`Tripwire`, populated from a sibling shape's captured site rather than a literal, so V5's
no-line-numbers property survives.

### L1 — Low: the FastQC allowance is broader than the rule it encodes

`is_fastqc_artifact` (`:201-204`) matches any file whose name ends `_fastqc.html` /
`_fastqc.zip`, anywhere under the root, with any stem. §2.2(b)'s rule is narrower:
`<primary-stem>_fastqc.{html,zip}` *beside the primary*. On a `--fastqc` run an unplanned
write named `*_fastqc.zip` would be swallowed — contrived today, since no other writer
produces that suffix, and the loose form is what lets the `-o qc` case work without
special-casing. My §4 fix removes the worst consequence of the breadth (an allowance that
takes the *whole* write set now fails), so this is a note rather than a defect. If tightening
is wanted, require the stem to match a created non-artifact file's stem.

### L2 — Low: `compare` keys every planned path twice

`:404-407` builds `planned_keys`, then `:412-420` recomputes `collision_key(&root.join(p))`
for the containment check. ~20 paths, so irrelevant to runtime; the second loop could iterate
`planned_keys` instead.

### L3 — Low: the `!active` stray-FastQC branch is a message, not a guard

`:362-369`. With no FastQC flag a `_fastqc.zip` is in `rest` and absent from `planned_keys`,
so T1 reports it regardless. The branch buys a clearer message and
`comparison_rejects_stray_fastqc_artifacts` pins that message, which justifies keeping it —
worth knowing it adds no coverage. `Vec<&&PathBuf>` at `:363` is a double reference that
`.copied()` would flatten.

### L4 — Low: case roots are never removed

`case_root` clears a stale root on entry, nothing clears it on exit, so a full run leaves 52
directories in `$TMPDIR`. `integration_output_collision.rs` has the same shape, so this
matches precedent; noting it only because it is a new increment of 52 × 2 profiles on CI
runners.

### L5 — Low: A2's residual is still unprobed

§8 A2 says `--fastqc_args "--svg"` (`fastqc.rs:81`) "was not probed; probe it once during
implementation". §13 does not record it being probed and no case passes `--svg`. If it writes
an `.svg` beside the primary, that is an unplanned write outside the allowance and the suite
would not see it. One case, or an explicit module-doc note that the allowance covers only the
two-suffix rule.

### L6 — Low: `compare`'s trailing bare `bool`

`compare(&root, &created, &planned, self.args, self.arm, self.fastqc_capable)` — `arm` is
properly an enum, but `fastqc_capable` is a bare `bool` in position 6. A two-variant enum
alongside `Arm` would read better at the call site. Cosmetic.

---

## 4. Fixes applied

All three are in `tests/integration_preflight_tripwire.rs`. Each is behaviour-neutral on the
current suite (verified: 50→51 tests, all green in both profiles, clippy and fmt clean).

**F1 — T3 now covers every pre-existing file, not only the declared `guarded` list**
(`:258-274`). Previously `hashes` was built from `self.guarded`, so T3 read as "no *declared*
input's bytes changed" while §3.3 states it as "no *pre-existing file's* bytes changed". The
two coincide in every case today — I checked all 41 — but the coincidence is maintained by
hand, per case: a future case that stages a pre-existing report or barcode file and forgets
to declare it would fall silently outside T3, which is the exact bug family this file
polices. `before` was already computed one line above, so hashing it costs nothing.
`guarded` is retained as the staging assertion.

```rust
let before = snapshot(&root);
for rel in self.guarded {
    assert!(before.contains(Path::new(rel)), "{}: guarded input {rel} was not staged", self.tag);
}
// T3 covers every pre-existing file, so a case that stages one without declaring
// it in `guarded` is still guarded.
let hashes: BTreeMap<PathBuf, u64> = before
    .iter()
    .map(|rel| { let p = root.join(rel); (p.clone(), hash_file(&p)) })
    .collect();
```

**F2 — the allowance can no longer swallow the entire write set** (`:355-361`). This is the
"a filter that removes everything" hole: `rest` is `created` minus FastQC artefacts, and
nothing required `rest` to be non-empty, so if every created file were taken as an artefact,
T1 would compare two empty sets and pass against *any* candidate list. Unreachable today
(every FastQC run also creates a primary), but it is the failure mode a future widening of
`is_fastqc_artifact` would open, and the plan's own anti-decay argument rests on the
allowance being self-policing. Now it is policed in that direction too.

```rust
// An allowance that swallowed the whole write set would leave T1 comparing empty sets.
if rest.is_empty() {
    return Err(format!(
        "the FastQC allowance took every created file ({allowed:?}), so T1 has nothing \
         left to check — narrow `is_fastqc_artifact`, do not widen it"
    ));
}
```

**F3 — new self-test `comparison_rejects_an_allowance_that_takes_everything`** (`:1125`),
because F2 is an assertion and this file's standard is that every assertion has a permanent
proof it can fire. Pure `compare` call, no binary run.

**F4 — the case-root containment check is now a path-boundary test** (`:411-414`).
`key.starts_with(&root_key)` was a raw string prefix, so `/tmp/tg_trip_a_1` is a prefix of
`/tmp/tg_trip_a_10/x`. Unreachable here (one PID per run, unique tags), but one character
from correct:

```rust
let root_prefix = format!("{root_key}{}", std::path::MAIN_SEPARATOR);
…
if !key.starts_with(&root_prefix) {
```

I did **not** touch: T5's comment or A10 (M1 — the fix is a design choice between two
options), the expected-site binding (M2 — restructuring), or `is_fastqc_artifact`'s breadth
(L1 — scope).

---

## 5. Verification I ran

All in an isolated worktree at `dev` @ `3008b37` + the `src/io.rs` patch + the test file, so
none of it is contaminated by the shared-tree mutation described in §0.

| | Check | Result |
|---|---|---|
| V1 | `cargo test --test integration_preflight_tripwire` (as reviewed) | 50/50, 10.55 s |
| V2 | `cargo test --release --test …` (as reviewed) | 50/50, 1.38 s |
| V3 | Full `cargo test` | **640 passed, 0 failed** — matches §13's claim |
| V4 | `cargo clippy --all-targets --release -- -D warnings` | exit 0, 12.5 s with dependency compilation, i.e. real work rather than a cache hit |
| V5 | `cargo fmt --all -- --check` | clean |
| V6 | Site coverage, independently | 21 shapes probed; every one hits the site the plan's table claims, one site per run, 11 distinct, orphans dump nothing |
| V7 | Release-binary hook inertness | 0 dump lines with the variable unset |
| V8 | `--version` / `-V` | em-dash provenance present on long form, absent from `-V` |
| V9 | No absolute source paths introduced | 0 references to the repo/worktree; all 424 absolute strings are pre-existing `~/.cargo/registry` dependency paths; `src/main.rs` `Location` string is relative |
| V10 | **T3 positive control** | Planted `std::fs::write(input, b"CLOBBERED")` after `run_single_file`; `se_trim_two_inputs` fails with `input …/a.fastq was modified by the run`, T1 green. Probe reverted, `git status` clean |
| V11 | **T1 positive control** (unplanned) | Reviewer B's `MUTANT` writer turned 11 cases red with the correct site and candidate list (§0) |
| V12 | After my fixes: fmt / debug / release / full suite / clippy | clean / 51 pass / 51 pass / **641 passed, 0 failed** / exit 0 |

Two notes on method, both from prior bites in this project:

- `cp` is aliased to `-i` here and silently blocked one copy into the worktree; I used
  `cat X > Y`.
- A clippy run finishing in under a second is a cache hit, not a pass. The runs above took
  12.5 s and 2.3–2.5 s with a `Compiling trim-galore` line, so they did real work.

---

## 6. Recommendations, by priority

**Critical** — none.

**High**

1. **H1 — give T3 a permanent self-test.** Extract the hash loop as a pure `unchanged(...)`
   and add `comparison_detects_a_clobbered_input`. ~15 lines. T3 is the sole substantive
   assertion on the three orphan cases and the only assertion in the file whose ability to
   fail rests on a probe that gets reverted. I verified it *does* fire (V10); this is about
   keeping that true.

**Medium**

2. **M1 — make T5 true, or make its claim true.** Either add an invocation boundary to the
   dump and assert exactly one call (~4 lines in `io.rs`), or narrow the comment and A10 to
   "no run reaches two distinct pre-flight sites". Leaving the stronger claim standing in a
   file that documents its own guarantees is the part I would not ship.
3. **M2 — pin the three arms that share site 2670.** Add `--implicon` and
   `--clump_only --paired` to `every_site_is_reached_by_some_case` and assert they collide
   with `--clock`, converting the three-way sharing from prose into an assertion. Optionally
   add `expect_site` to `Tripwire`, populated from a sibling shape rather than a literal.

**Low**

4. **L5 — probe `--fastqc_args "--svg"` once** and either add a case or record in the module
   doc that the allowance is the two-suffix rule only. This is an explicit A2 residual the
   plan asked to close during implementation and §13 does not record closing.
5. L1 (tighten the allowance to stem-matched artefacts), L2 (key each planned path once),
   L3 (`.copied()`), L4 (clean up case roots), L6 (enum instead of the trailing `bool`) — all
   optional.

**Housekeeping, not a review finding**

6. Confirm `src/main.rs` is clean before committing. It was reverted during my review, but
   the `REVIEW-B MUTATION` line was in the branch's working tree for part of it.
