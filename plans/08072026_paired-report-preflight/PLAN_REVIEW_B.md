# Plan review B — #388 paired report pre-flight

**Reviewer:** B (independent; no coordination with Reviewer A)
**Target:** `plans/08072026_paired-report-preflight/PLAN.md` (v1, 2026-08-07)
**Repo state:** `dev` @ `c4599f5`, clean tree, plan not implemented
**Method:** read the cited source, then ran `target/release/trim_galore` (current with `dev`)
against each scenario the plan asserts. Scratch:
`/private/tmp/claude-501/.../scratchpad/prev388b/`. Filesystem here is
case-**insensitive** (APFS), confirmed by `ls R1.fq` resolving to `r1.fq`.

---

## Verdict

The **fix itself is right**: 6 lines in the right place, correct gate, and it closes a
harm I reproduced. The **validation section is not**: three of the four listed tests
cannot pass as specified, the negative control in §5 is void, and the one scenario the
fix actually exists for is not tested at all. Assumption A1 is false, and its failure
mode is a user-visible behaviour change the plan does not acknowledge.

Root cause of most of the below, stated once: **paired primaries carry a positional
discriminator (`_val_1` / `_val_2`) that report names do not.** #385's SE
"reports-are-finer-than-primaries" argument held because SE has one primary per input.
Inheriting it into paired mode is invalid. Every finding here flows from that.

---

## What checks out

| Plan claim | Verdict |
|---|---|
| Insertion point: per-chunk `candidates`, token `Pre-flight across pairs before any I/O` | **Correct** — `src/main.rs:731`; `candidates` built per chunk (741–764), `planned.extend` at 765 |
| Multi-pair feeds **one** `preflight_output_collisions` call | **Correct** — `planned` declared 732, extended in-loop, single call at `src/main.rs:767` |
| Writer gate is `!cli.no_report_file` at `main.rs:1558` | **Correct**, line-accurate, and it covers **both** `.txt` and `.json` (1570–1578) |
| SE precedent mirrors exactly | **Correct** — `main.rs:90-93` is `if !cli.no_report_file { report_name; json_report_name }`, identical shape |
| Harm in §2.1 (`--paired r1.fq R1.fq` self-pairs, one report silently lost) | **Reproduced** — see below |
| §2.2 "same four lines" for the uBAM twin | **Correct** — same four paths, same gate (verified in `run_ubam_output_paired_two_files`) |
| A3 (FastQC out of scope) | Reasonable, unchanged from #385 |

### The harm is real (§2.1 confirmed)

```
$ trim_galore --paired --dont_gzip r1.fq R1.fq     # only r1.fq on disk
exit=0
r1_val_1.fq  R1_val_2.fq  r1.fq  r1.fq_trimming_report.json  r1.fq_trimming_report.txt
```

Exit 0, self-paired, **one** report pair on disk — R1's report landed on r1's. Both
`_val_` files exist because the suffix differs. This is the #383 harm class surviving,
exactly as the plan says, and the fix does catch it (report keys fold to one
`collision_key`). The fix is worth making.

---

## Critical

### C1. V1 cannot pass, and cannot fail — `Cli::validate` gets there first

The plan's V1 (`--paired ./r1.fq r1.fq`) never reaches the pre-flight.
`cli.rs:537` compares `path_identity_key(chunk[0]) == path_identity_key(chunk[1])`, and
`path_identity_key` **absolutises** (`lexical_normalise`, `io.rs:78`), so `./r1.fq` and
`r1.fq` are already one key. Measured on today's binary, no fix applied:

```
$ trim_galore --paired ./r1.fq r1.fq
exit=1
Error: Read 1 and Read 2 appear to be the same file: ./r1.fq. ...
```

Consequences, all three fatal to V1 as written:

1. It **already exits non-zero** → with loose assertions the test cannot fail, so it
   guards nothing.
2. Written in house style it **cannot pass even after the fix**:
   `assert_rejected_cleanly` asserts `stderr.contains(PREFIX)` where `PREFIX =
   "Output path collision (case-insensitive, for APFS/NTFS safety)"`
   (`tests/integration_output_collision.rs:91,100-103`). Validate's message contains no
   such prefix and short-circuits before the pre-flight.
3. §5's negative control ("comment the block → V1 and V2 must fail") is therefore
   **void** — V1 passes or fails identically with and without the change.

§2.1's supporting claim is also wrong: "`./r1.fq_trimming_report.txt` and
`r1.fq_trimming_report.txt` absolutise to one key" is *true* (I traced it: `Path::new
("r1.fq").parent()` is `Some("")`, not `None`, so `report_name` joins onto empty and
`collision_key` absolutises both to `<cwd>/r1.fq_trimming_report.txt`) — but it is
unreachable, because the identical normalisation in `path_identity_key` rejects the
inputs first. **There is no filesystem-independent two-spelling input that reaches the
report route**: `path_identity_key` and `collision_key` differ *only* in case folding,
so any pair of spellings that folds together in report space is also caught by validate
unless they differ in case — and case-only aliasing is precisely the FS-dependent case.

**Fix — replace V1 with the harm case, made portable by writing both spellings:**

```rust
write_fastq(&dir.join("r1.fq"), "R1");
write_fastq(&dir.join("R1.fq"), "R2");   // one file on APFS, two on ext4
let (ok, stderr) = run_in(&dir, &["--paired", "--dont_gzip", "r1.fq", "R1.fq"]);
assert!(!ok); assert!(stderr.contains(PREFIX) && stderr.contains(DUP_MSG));
```

Portable because the pre-flight is purely lexical: on a case-insensitive FS the two
names are one file (the reproduced harm); on a case-sensitive FS they are two files
whose report keys still fold. Rejected on both after the fix; **accepted on both
before** it (I verified the case-insensitive half above), so it can genuinely fail.
Do not create only `r1.fq` and reference `R1.fq` — `main.rs:209-213` format-detects
*every* input, so on a case-sensitive FS that dies on a missing file, not a collision.

### C2. V4 is impossible as specified

V4 expects `--no_report_file` + V1's inputs to **succeed**. It cannot: the same
`cli.rs:537` check is gate-independent.

```
$ trim_galore --paired --no_report_file ./r1.fq r1.fq
exit=1
Error: Read 1 and Read 2 appear to be the same file: ./r1.fq. ...
```

Rebase V4 on C1's inputs instead: `--paired --no_report_file r1.fq R1.fq` must
**succeed** (reports are the only colliding keys, so removing them from the candidate
list must let it through). That version genuinely pins the gate — and note it pins it in
the direction that matters, since a gate-less implementation would reject it.

### C3. V2 is unreachable — the format guard precedes the pre-flight

`main.rs:209-213` runs `detect_input_format` over **all** inputs, ~550 lines before the
paired pre-flight. A fed-back report dies there:

```
$ trim_galore --paired a_R1.fq a_R2.fq x.fq a_R1.fq_trimming_report.txt
exit=1
Error: Input 'a_R1.fq_trimming_report.txt' is not recognised as FASTQ (plain or gzipped) or unaligned BAM
report md5 before == after   # byte-intact
```

So V2 asserting `ALIAS_MSG` can never pass (wrong message, before *and* after the fix),
and asserting only "refused + byte-intact" can never fail. **§2.1's third paragraph is
factually wrong**: these inputs are *already* refused and were *never* overwritten, so
the fix does not add the protection the plan credits it with. The output-vs-input half
of the change is effectively dead on the paired path — a report path can only equal an
input if that input is a valid FASTQ *named* `*_trimming_report.txt`, which the format
guard otherwise excludes. Harmless to include (it mirrors SE, cheap), but **drop V2 or
rewrite it around a valid-FASTQ file deliberately named `*_trimming_report.txt`**, and
correct the §2.1 claim either way.

### C4. A1 is false — a non-alias over-rejection exists, and it removes working runs

A1 claims report candidates "cannot reject a pair the primaries would have accepted —
except where two inputs alias one file or a report aliases an input". Counterexample,
two genuinely distinct files, no aliasing:

```
$ trim_galore --paired --dont_gzip -o out A/reads.fq B/reads.fq
exit=0
out/: reads_val_1.fq  reads_val_2.fq  reads.fq_trimming_report.{txt,json}
```

Primaries are distinct (`_val_1` vs `_val_2`); the two reports are **not** — one report
is silently lost today, and after the fix this invocation becomes a hard error. That is
arguably the *right* outcome (the project already rejects the SE analogue on purpose —
`se_trim_rejects_same_basename_across_dirs_with_output_dir`, test file line 175), but:

- A1 as written is wrong and must be restated. The exact residual class is: **one pair
  whose R1 and R2 share a filename in different directories, written to a common
  directory** (via `-o`, or one input already in the cwd).
- It **removes a currently-working invocation**, so it needs a CHANGELOG line and a
  deliberate sign-off, not silence.
- I checked the neighbouring cases the brief raised and they are all fine: multi-pair
  with repeated basenames + `-o` is *already* rejected on primaries
  (`out/r1_val_1.fq and out/r1_val_1.fq`); `--basename` with >1 pair is rejected by
  validate ("ambiguous output naming"); single pair + `--basename` without `-o` still
  succeeds; same-basename in different dirs *without* `-o` still succeeds (reports land
  beside their own inputs). So the class above is the only one.
- **V3 does not detect it.** "Two ordinary distinct pairs" have distinct report names
  and pass with or without the fix. Add this case as its own test (expecting rejection).

---

## Important

### I1. "No new false positive is possible" is wrong

§2.1 asserts this from D13 test H. On a **case-sensitive** FS, `sample_R1.fq` and
`sample_r1.fq` are two real files with two real distinct reports, yet `collision_key`
folds case and the run is now rejected. That is a new false positive on Linux — the
documented #385 trade-off (`io.rs:45-47`: "on opt-in case-sensitive volumes this may
false-positive"), which is fine, but the plan claims *immunity*, and that claim is what
justifies not testing the boundary. Restate as: residual false-positive class is inputs
differing only in ASCII case on a case-sensitive FS, accepted per #385. (C1's test
doubles as documentation of it.)

### I2. The uBAM-out twin has the identical live bug — the scope call is incoherent

§2.2 defers it. But §1 frames the whole plan as *closing* the asymmetry #385 left, and
deferring here opens a new one of the same shape. Reproduced on today's binary:

```
$ trim_galore --paired --output-format ubam r1.fq R1.fq     # FASTQ in, uBAM out
exit=0
r1_val.bam  r1.fq  r1.fq_trimming_report.{txt,json}     # one report pair, one lost
```

It is *more* exposed than the FASTQ path, not less: its pre-flight plans only **one**
primary per pair (`paired_bam_output_name`, `main.rs:1849-1860`), so fewer keys guard
it. The plan's own "same four lines" is accurate (same paths, same `!cli.no_report_file`
gate) and the gate is identical, so there is no technical reason to defer.
**Recommend folding it in** — 4 lines at `main.rs:1849-1860`, one extra test.

Two corrections to §2.2 while there: (a) the pointer `:1825` is the `run_ubam_output`
function header, not the insertion site — that is the per-chunk loop at 1849–1860;
(b) two-file *BAM input* is rejected upstream by `reject_bam_format_mismatch_in_pair`,
so this path is reachable only with FASTQ input + `--output-format ubam` — worth saying,
since it changes how the test is written. The single-file interleaved path is genuinely
safe, but for a reason the plan does not give: it has one input, so its `.txt`/`.json`
names cannot collide with each other or with the input.

### I3. V3 is weaker than the file's own standard

`tests/integration_output_collision.rs:11-14` requires acceptance cases to "assert
content or filename rather than mere existence". V3 says "run to completion". Use
`count_reads_from` on both `_val_` files and assert all four report files exist with
distinct names — otherwise a candidate list built from the wrong namer still passes.

### I4. No test covers the harm the fix exists for

Following C1–C3, the plan's post-fix test set contains **zero** tests that fail without
the change. Minimum viable set: C1's case-alias rejection, C2's rebased gate test, C4's
same-basename-with-`-o` rejection, a strengthened V3, and (if I2 is accepted) the
uBAM-out twin. Re-run the §5 negative control against *that* set.

---

## Optional

- **O1.** CHANGELOG (Step 3): add the removed-invocation note from C4. The #383 entry
  already sets this precedent; a bare "collisions are now refused" undersells it.
- **O2.** The new comment says a report can "overwrite an input" — per C3 that is
  unreachable on this path. Trim to the reachable half: reports are outputs too, and
  two inputs differing only in case name one report.
- **O3.** §6's self-review claims the multi-pair accumulation was checked and it was
  right — worth keeping. But it also says "Risk: none identified beyond over-rejection,
  which A1/V3 covers", and A1/V3 demonstrably do not cover it (C4). Update.
- **O4.** Consider whether report names *should* carry the `_1`/`_2` discriminator, the
  way primaries do. That would fix C4's over-rejection class at the source rather than
  converting it into an error — but it breaks MultiQC's expected filenames, so I am
  **not** recommending it; noting it so the trade-off is on the record.

---

## Summary of action items

| # | Severity | Item |
|---|---|---|
| C1 | Critical | V1 unreachable (`validate` cli.rs:537 fires first); replace with case-alias harm test; §5 negative control void |
| C2 | Critical | V4 impossible; rebase on `r1.fq R1.fq` inputs |
| C3 | Critical | V2 unreachable (format guard `main.rs:209-213`); correct §2.1 ¶3's false claim |
| C4 | Critical | A1 false; same-basename R1/R2 + `-o` is a non-alias over-rejection that removes working runs; V3 misses it |
| I1 | Important | "No new false positive is possible" is wrong on case-sensitive filesystems |
| I2 | Important | uBAM-out twin has the identical live bug (reproduced); fold in the 4 lines |
| I3 | Important | V3 must assert content/filename per the test file's own standard |
| I4 | Important | As specified, no test fails without the change |
| O1–O4 | Optional | CHANGELOG note; comment wording; §6 update; report-naming trade-off on record |
