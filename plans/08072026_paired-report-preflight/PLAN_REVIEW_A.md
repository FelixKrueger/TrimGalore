# Plan Review A — Add report paths to the paired pre-flight (#388)

**Plan:** `plans/08072026_paired-report-preflight/PLAN.md` (v1, 2026-08-07)
**Target:** `dev` @ `c4599f5`, clean tree
**Reviewer:** A (independent; no coordination with B)
**Method:** source read of `src/main.rs`, `src/io.rs`, `src/cli.rs` plus empirical runs of
`target/release/trim_galore` (2.3.0, `533582a`) in
`/private/tmp/claude-501/…/scratchpad/prev388a/`.

**Verdict:** the *diagnosis* and *Step 1* are correct and land in the right place. The
**validation section is broken**: three of the four proposed tests (V1, V2, V4) fail or
are vacuous *even with the fix correctly applied*, because two earlier guards preempt
the pre-flight. As written, the plan ships a correct 6-line fix with **no test that
exercises it**, and its §5 negative control cannot hold. The real reproducer is not in
the plan; it is given below, verified.

---

## 1. What checks out

- **Insertion point exists as described.** `src/main.rs:730-767` is the `if cli.paired`
  FASTQ pre-flight. `candidates` is declared at `:741`, the passthrough block ends at
  `:764`, `planned.extend(candidates)` is at `:765`. Inserting Step 1's block between
  `:764` and `:765` compiles: `output_dir` (used `:736`), `cli.no_report_file`, and the
  `naming` alias are all in scope, and `&PathBuf` → `&Path` coerces at the call.
- **Multi-pair really feeds ONE pre-flight call.** The `for chunk in cli.input.chunks(2)`
  loop accumulates into `planned`; the single `preflight_output_collisions` is at `:767`,
  after the loop. §6's cross-pair claim is structurally correct — and the cross-pair hole
  is live today (§3.4 below).
- **The writer's gate claim is exactly right, including JSON.** `main.rs:1558` is
  `if !cli.no_report_file {`, the `write_paired_reports` call at `:1585` is inside it, and
  the descriptors at `:1570-1576` use `report_name`/`json_report_name` for R1 and R2 —
  the same four names Step 1 pushes. `write_paired_reports` writes txt (`:2392`) and json
  (`:2417`) unconditionally *inside* that one gate; there is no separate JSON flag. A2 holds.
- **The harm is real and worse than the plan claims.** Verified silent data loss on
  current `dev` (§3.4).
- **Step 3's CHANGELOG anchors exist** — `#### Bug fixes` at `CHANGELOG.md:6`, #383 at
  `:64` and `:115`.

---

## 2. Critical findings

### C1 — V1 cannot pass. `Cli::validate` rejects `--paired ./r1.fq r1.fq` today, before the pre-flight.

`path_identity_key` (`io.rs:78`) calls `lexical_normalise`, which **absolutises**
(`io.rs:55-56`). So `./r1.fq` and `r1.fq` are already ONE key, and `cli.rs:537` bails.
The source comment at `cli.rs:535` states this outright:

> `#383 — keyed like the output-collision pre-flight, so ./r1.fq r1.fq cannot pass as a pair.`

Verified:

```
$ trim_galore --paired ./r1.fq r1.fq
Error: Read 1 and Read 2 appear to be the same file: ./r1.fq. Did you mean to pass distinct R1 and R2 files?
```

V1 asserts "non-zero exit **with the collision message**". The actual message is the
same-file message, so **V1 fails with the fix applied**, and §5's negative control
("comment the block → V1 must fail") is meaningless — V1 fails either way.

§2.1's parenthetical ("This also upgrades the D13-noted self-pair from silent to refused
via a second, spelling-based route") and V1's "Runs on any filesystem" are the root
error. The report-key arithmetic the plan asks to be traced *is* correct —
`report_name(./r1.fq, None)` → `./r1.fq_trimming_report.txt` (parent of `./r1.fq` is
`.`), `report_name(r1.fq, None)` → `r1.fq_trimming_report.txt` (parent is `""`, so the
`unwrap_or(".")` never fires), and `collision_key` absolutises both to the same string —
but it is unreachable: control never gets to `:767`.

### C2 — V2 cannot fail. The format detector rejects `*_trimming_report.txt` before the pre-flight.

`main.rs:212` maps `detect_input_format` over **all** inputs, well before `:767`.
V2's exact shape, verified:

```
$ trim_galore --paired a_R1.fq a_R2.fq x.fq a_R1.fq_trimming_report.txt
Error: Input 'a_R1.fq_trimming_report.txt' is not recognised as FASTQ (plain or gzipped) or unaligned BAM
```

Identical before and after the fix. V2 asserts "the alias wording", so it fails with the
fix applied; had it asserted only "refused", it would be vacuous. Either way it cannot
observe the new code.

**Corollary — §2.1 ¶3's claimed benefit is near-unreachable.** For a report path to hit
the output-vs-input branch, an *input* must be named `*_trimming_report.{txt,json}` **and**
carry valid FASTQ/uBAM content. A genuine prior report never does. So "a prior run's
report fed back in are refused instead of overwritten" is not the win the plan states;
the only reachable case is a FASTQ file with a report-shaped name. Keep the fix (it is
free), but drop the claim or restate it honestly.

### C3 — V4 cannot pass, for C1's reason.

`--no_report_file` does not reach `validate`'s same-file check. Verified:

```
$ trim_galore --paired --no_report_file ./r1.fq r1.fq
Error: Read 1 and Read 2 appear to be the same file: ./r1.fq. …
```

V4's *intent* (pin the gate) is sound and A2 deserves a test — it just needs inputs that
survive `validate`. With the reproducer from C4 it works today: verified exit 0, output
dir holds only `r1_val_1.fq` + `R1_val_2.fq`, no reports.

### C4 — No proposed test exercises the fix. Here is a verified reproducer.

Two **genuinely distinct** input files whose filenames fold equal, sharing an output dir.
`validate` passes (`path_identity_key` is case-sensitive), primaries differ
(`_val_1` vs `_val_2`), reports collide. On current `dev`:

```
$ trim_galore --paired a/r1.fq b/R1.fq -o out ; echo $?
0
$ ls out
r1_val_1.fq  R1_val_2.fq  r1.fq_trimming_report.json  r1.fq_trimming_report.txt
$ grep 'Input filename' out/r1.fq_trimming_report.txt
Input filename: R1.fq          # ← R2's report overwrote R1's. Two inputs, one report.
```

Exit 0, silent loss of R1's report. Step 1's code catches it: both candidates fold to
`<abs>/out/r1.fq_trimming_report.txt`. This is the test V1 should have been. It requires
a case-insensitive FS to be a *true* positive but is rejected unconditionally (see I1) —
so as a **rejection** test it passes on Linux CI too.

---

## 3. Important findings

### I1 — A1's "strictly finer" is wrong, and §2.1's "no new false positive is possible" is false.

Primary keys encode **position** (`_val_1` vs `_val_2`); report keys do not. In that
dimension the report key is *coarser*, which is precisely why C4's case exists — the
plan's own headline scenario contradicts its assumption. The correct statement of the
new rejection set is: *two inputs whose filenames fold equal in the same output dir*.

`collision_key` (`io.rs:84`) always case-folds; the pre-flight has no filesystem
awareness. So `--paired r1.fq R1.fq` with two genuinely distinct files **will be newly
refused on a case-sensitive filesystem** — a real false positive. That is acceptable and
already the documented design trade-off (`io.rs:45-47`: "on opt-in case-sensitive APFS
volumes this may false-positive, but the penalty is a loud early error rather than silent
data loss"), but the plan must say so instead of claiming immunity. Rewrite A1 accordingly.

Checked the other over-rejection routes the plan should have considered, all benign:
`--basename` (reports keep per-input names, but with one pair the primaries are always
distinct and with multiple pairs the primaries *already* collide → no new rejection);
`-o` absent (reports land next to each input, so cross-dir inputs stay distinct);
R1/R2 of one pair sharing a *stem* but differing in extension (`sample.fq` + `sample.fastq`
→ report names distinct). No over-rejection beyond I1's fold-equal family.

### I2 — §2.2's out-of-scope call leaves an identical, live bug shipping.

The paired-uBAM twin has the same hole, verified on current `dev`:

```
$ trim_galore --paired --output-format ubam a/r1.fq b/R1.fq -o out ; echo $?
0
$ ls out
r1_val.bam  r1.fq_trimming_report.json  r1.fq_trimming_report.txt   # one report pair, two inputs
```

Its pre-flight is `main.rs:1851-1860` and its report names are at `:2163-2169`
(`run_ubam_output_paired_two_files`), both using `report_name`/`json_report_name` on
`input_r1`/`input_r2`. Deferring is a defensible scope call, but "out of scope because
this mirrors the SE FASTQ precedent" is not a reason — the harm class §2.1 invokes is
equally reachable here. Either cover it (four more lines, same shape) or say plainly
that a known-live instance is being left for a follow-up issue.

Two inaccuracies in that paragraph: the `:1825`-area pointer is the *dispatcher*
(`run_ubam_output` starts at `:1829`; its paired pre-flight loop is `:1851-1859`), and
**"the same four lines" will not compile verbatim** — that loop pushes to `planned`,
not `candidates`.

### I3 — The search token is ambiguous: two hits.

`Pre-flight across pairs before any I/O` matches `main.rs:731` (intended) **and**
`main.rs:2458` (`run_specialty_paired`, for `--clock`/`--implicon`). Since §2.2 gives the
insertion point by token rather than line, pin it to `:730-767` / "the `if cli.paired`
FASTQ branch". The `:2458` site has no `candidates` binding, so a wrong landing fails to
compile rather than silently misbehaving — but it costs the implementer a cycle.

### I4 — The asymmetry is wider than §1 states: the uBAM SE path has the hole too.

`planned_secondary_outputs` is called **exactly once**, at `main.rs:833` — the FASTQ SE
loop. The uBAM SE pre-flight (`:1899-1904`) plans only BAM names, while
`run_ubam_output_single` writes reports at `:1929`/`:2017`. So #385 fixed FASTQ-SE only.
§1's "the SE path gained `planned_secondary_outputs`" should read "the SE **FASTQ** path".
(In practice the uBAM SE primaries catch the fold-equal case — verified `exit=1` for
`--output-format ubam a/r1.fq b/R1.fq -o out`, because `_trimmed.bam` names fold too —
so only the C2-class residual is exposed. Worth a sentence, not necessarily a fix here.)

### I5 — V3 is too weak to be A1's guard, and no test covers the cross-pair route.

"Two ordinary distinct pairs still run" cannot detect over-rejection: ordinary pairs have
nothing near the fold-equal boundary. To actually guard A1, V3 should pin the shapes
adjacent to it — `--basename` with one pair, and same-stem-different-extension R1/R2 —
and the plan should state that the case-differing-distinct-inputs shape is now
**rejected** (I1), not accepted.

Separately, the cross-pair route §6 claims to cover is untested and live. Verified on
current `dev` — four inputs, three report pairs, exit 0:

```
$ trim_galore --paired a/r1.fq a/r2.fq b/R2.fq b/x.fq -o out ; echo $?   # 0
$ ls out | grep trimming_report.txt
r1.fq_trimming_report.txt  r2.fq_trimming_report.txt  x.fq_trimming_report.txt
# b/R2.fq's report overwrote a/r2.fq's
```

---

## 4. Efficiency

No concern. Four extra `PathBuf`s per pair; `preflight_output_collisions` is one
`HashMap` pass over `planned`. `report_name`/`json_report_name` do no filesystem I/O
(`io.rs:485-514` are pure string joins), so the pre-flight stays I/O-free as documented.

## 5. Alternatives

- **Preferred:** have the paired branch call a shared helper (the existing
  `planned_secondary_outputs`, or a small `paired_secondary_outputs`) instead of
  open-coding the push. That is the #384 "two functions reading one property must fold
  together" lesson §6 invokes, applied properly — it makes writer/guard drift structural
  rather than test-enforced, and gives I2/I4 a single place to fix.
- Considered and rejected: adding `_val_1`/`_val_2`-style position to report names would
  remove I1's false positives, but breaks the MultiQC-visible report naming and Perl
  parity. Not worth it.

---

## 6. Action items

### Critical
1. **Rewrite V1** to C4's verified reproducer (`--paired a/r1.fq b/R1.fq -o out`,
   distinct files, fold-equal names), asserting the *collision* wording. The
   `./r1.fq r1.fq` scenario is already refused by `cli.rs:537` and must be dropped or
   re-labelled as a pre-existing `validate` guard.
2. **Drop or rebuild V2.** `main.rs:212` preempts any real `*_trimming_report.txt` input.
   If kept, the fixture must be FASTQ *content* under a report-shaped *name*, and §2.1 ¶3's
   benefit claim must be restated to match.
3. **Rewrite V4** onto C4's inputs + `--no_report_file` (verified exit 0 today).
4. **Re-derive §5's negative control** after 1–3; as written it cannot hold.

### Important
5. **Fix A1 and §2.1.** Report keys are coarser in the position dimension; new rejections
   are exactly "fold-equal filenames in one output dir"; false positives on case-sensitive
   filesystems **are** possible and are the accepted `io.rs:45-47` trade-off.
6. **Resolve §2.2.** Cover the paired-uBAM twin (`:1851-1859` pre-flight, `planned` not
   `candidates`) or state that a verified-live instance is being deferred, with an issue.
   Fix the `:1825` pointer and the "same four lines" wording.
7. **Add a cross-pair rejection test** (verified live hole) — it is what §6 claims.
8. **Pin the insertion point by line** (`:730-767`), not the two-hit token.
9. **Strengthen V3** per I5; correct §1 to "SE **FASTQ** path" per I4.

### Optional
10. Refactor toward a shared secondary-outputs helper (§5) rather than an open-coded push.
11. Note the uBAM SE residual (I4) in §2.2's scope paragraph.

---

## 7. Harness confidence

Every "pass" above was shown capable of failing: the same binary produced `exit=0` for
the four accepted shapes and `exit=1`/error text for the rejected ones, and the data-loss
claim is evidenced by `Input filename: R1.fq` inside a file named
`r1.fq_trimming_report.txt` — not by a file count alone. Case-sensitive-filesystem
behaviour (I1) is derived from `collision_key`/`norm_path` source, not executed; macOS
APFS here is case-insensitive.
