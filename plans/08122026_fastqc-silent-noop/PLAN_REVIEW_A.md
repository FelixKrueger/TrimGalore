# PLAN_REVIEW_A — #421 `--fastqc` silent no-op (Reviewer A)

**Reviewed:** `plans/08122026_fastqc-silent-noop/PLAN.md` v1 (2026-08-12)
**Repo:** `/Users/fkrueger/Github/TrimGalore`, `dev` @ `3ffcef4`, clean tree
**Method:** dispatch read of `src/main.rs` (all early-return arms), `src/cli.rs::validate`,
`src/clump_only.rs`, `tests/integration_preflight_tripwire.rs`, `.github/workflows/ci.yml`, `docs/`,
plus **empirical probes** with the committed release binary (stamped `1fc00f6` — two commits behind
`3ffcef4`; `4e4a762` touched only `src/bam.rs` + two uBAM test files and `3ffcef4` is a docs npm bump,
so neither can affect FastQC gating). No file in the main tree was modified; probes ran in `$TMPDIR`.

**Verdict: the plan is sound and its highest-risk claim is correct.** One factual error in §2.1/A4
matters (a **sixth** silent no-op exists that the plan asserts does not), the tripwire deletion
inventory is incomplete and unreplaced, and the docs/help scope is under-counted.

---

## 1. Logic review

### 1.1 Arm 5 — verified independently, from the dispatch, and empirically. The plan is right.

The claimed condition `paired && input.len() == 1 && output_format != UBam && !clump_only` is correct,
and the `output_format` term is genuinely load-bearing. Three independent witnesses:

1. **Dispatch order.** `main.rs:845` returns into `run_ubam_output` for `OutputFormat::UBam` *before*
   `main.rs:853`'s `if cli.paired && cli.input.len() == 1` → `run_paired_ubam_single_file`
   (`1810-2003`, no `fastqc::run` anywhere in it). Inside `run_ubam_output`, `:2006` routes the same
   shape to `run_ubam_output_paired_single_file`, which runs FastQC at `:2521-2522`. So one term —
   `output_format` — is the entire difference between the incapable arm and its capable twin.
2. **CI.** `ci.yml:1165-1182` runs `--paired --output-format ubam --fastqc test_files/ubam_paired_test.bam`
   and asserts `ubam_paired_test_val_fastqc.{zip,html}` exist **and** that exactly one zip appears.
   Dropping the `output_format` term reddens that step.
3. **The tripwire already encodes the distinction.** `orphan_paired_fastq_from_interleaved_bam`
   (`:1021`, `--paired ip.bam`) is `.fastqc_incapable()`; `orphan_paired_ubam_from_interleaved_bam`
   (`:1048`, `--paired --output-format ubam ip.bam`) is **not**. Two adjacent cases, one flag apart.

Empirical confirmation (probes in `$TMPDIR`, `ubam_paired_test.bam` copied in as `ip.bam`):

| Invocation | exit | files created | `grep -ci fastqc stderr` |
|---|---|---|---|
| `--paired --fastqc ip.bam` | 0 | `ip_val_{1,2}.fq` + 4 reports, **no `_fastqc.*`** | **0** |
| `--paired --output-format ubam --fastqc ip.bam` | 0 | `ip_val.bam` + **`ip_val_fastqc.{html,zip}`** | — |
| `--hardtrim5 20 --fastqc x.fastq` | 0 | `x.20bp_5prime.fq` only | **0** |
| `--clock --paired --fastqc r1.gz r2.gz` | 0 | two `*clock_UMI*` files only | **0** |
| `--implicon=8 --paired --fastqc r1.gz r2.gz` | 0 | two `*8bp_UMI*` files only | **0** |

**Over/under-inclusion of the four-term condition.** I could not construct an
accepted-today invocation that it refuses, other than `--paired --fastqc <single FASTQ>` — which
already fails today (`Error: --paired with a single input file is only legal if that file is a uBAM…`,
verified), so nothing working is lost; only the message degrades (→ §1.5). Nor could I find a
silently-ignoring invocation in the `paired && len == 1` family that it misses. The `!clump_only`
exclusion behaves as described: `--clump_only --paired --fastqc ip.bam` still emits
`main.rs:597`'s "requires two FASTQ input files… add `--output-format ubam`" (verified).

### 1.2 A1 — the five-arm set is complete for *reachability*, but the plan misses a sixth silent no-op

`main()`'s dispatch early-returns are exactly: `:493` hardtrim5, `:524` hardtrim3, `:546` clock,
`:578` implicon, `:828` clump_only, `:846` uBAM output, `:864` paired-single-file-FASTQ, then the
paired loop and the SE loop. That is the whole surface; the plan's five arms cover every arm from
which `fastqc::run` is unreachable. Multi-pair specialty input creates **no** extra arm — the guards
key on the flag, not the file count — and `--clumpify` is not an arm at all (it routes inside
`run_single_file` / `run_paired` / `run_paired_ubam_single_file`, `main.rs:1316/1548/1847`).

**But `--clump_only` gates on `if fastqc` alone, not on the disjunction.** `clump_only.rs:373`,
`:553`, `:844`, `:1094` all read `if fastqc`, and `main.rs` passes the bare `cli.fastqc` positionally
at `:636`, `:671`, `:715`, `:778`, `:819`. So `--fastqc_args` **without** `--fastqc` is silently
ignored on all four `--clump_only` drivers. Verified:

```
$ trim_galore --clump_only --fastqc_args "--quiet" x.fastq   # rc=0
x_clumped.fq  x.fastq_clumping_report.txt                     # no _fastqc.*, 0 mentions of "fastqc"
$ trim_galore --clump_only --fastqc x.fastq                   # control: artefacts appear
$ trim_galore --fastqc_args "--quiet" x.fastq                 # control: SE trim DOES produce them
```

This is not a nitpick about §2.1's bookkeeping — it is **§2.3's own argument applied to a case the
plan declares safe**: "a refusal keyed on `cli.fastqc` alone would leave `--fastqc_args` silently
ignored — the same defect, one flag over." That defect is live in `clump_only.rs` today. Consequences
for the plan as written:

- **A4 is false as stated.** The disjunction is the gate at the **5 gate expressions in `main.rs`**
  (covering 7 of the 12 `fastqc::run` calls), not "all 12 call sites".
- §2.1's table lists 11 sites under a claimed 12; the missing one is **`clump_only.rs:555`**, the R2
  call inside `clump_only_paired`.
- The plan's acceptance twin `--clump_only --fastqc a.fq` is fine, but if the implementer adds the
  symmetric `--clump_only --fastqc_args …` acceptance test (the natural thing to do given §2.3), that
  test **pins a silent no-op**.
- `cli.rs:196` — `--clump_only`'s own help text — explicitly lists `--fastqc` among the flags it
  "composes with", and `cli.rs:333` says `--fastqc_args` "Implies --fastqc". Both are in-repo promises
  that the clump_only path breaks.

The fix is one line and it *simplifies* the plan: add `Cli::fastqc_requested()` (`self.fastqc ||
self.fastqc_args.is_some()`), use it for the five new guards **and** at the five `main.rs`
pass-through sites. §4's "no helper" decision is about a *message* helper and is right; an accessor is
a different object and is what makes the disjunction one fact instead of nine copies.

### 1.3 A6 / the tripwire deletion — both halves verified; the inventory and the replacement are not

Both halves of A6 hold. All eight `.fastqc_incapable()` cases (`:810`, `:821`, `:832`, `:843`, `:856`,
`:879`, `:1024`, `:1039`) pass **no** FastQC flag, so `active` is false and the `(true, false)` branch
cannot fire today; and they map onto exactly the plan's five arms (hardtrim5 ×2, hardtrim3 ×2, clock,
implicon, `--paired ip.bam` ×2 both FASTQ-output). No arm outside #421's scope is marked. There is
also no `comparison_*` self-test driving that branch, so nothing else depends on it.

**Does deleting it weaken the surviving rule?** No. The partition is `active && is_fastqc_artifact(p)`
and the surviving `(true, true)` rule ("a FastQC flag is set but no artefact was written") is
unchanged for every capable case; `fastqc_capable` only ever *relaxed* the rule.

Two problems remain:

- **The deletion inventory is incomplete.** `compare()` takes `fastqc_capable` as its **6th
  positional parameter**, and six self-tests call `compare()` directly with a trailing `true`:
  `:1104-1111`, `:1120`, `:1182`, `:1197`, `:1210-1217`, `:1255` — plus the pass-through at `:311`.
  Dropping the parameter touches all seven. These are compile errors, not silent bugs, but the plan
  enumerates its deletion set precisely ("the field, the builder and its eight call sites, and the
  `(true, false)` match arm") and an implementer working to that list will be surprised by the diff.
- **Nothing replaces it.** The marker is the only in-repo record that those five arms cannot run
  FastQC. Post-fix, a case that pairs a FastQC flag with one of them fails on `run()`'s exit-0 assert
  (`:291`) — a *reactive* failure that only exists if someone happens to write such a case. The
  positive form costs almost nothing and is strictly stronger than what is being deleted: five cases
  asserting non-zero exit + the message + **nothing written**, following
  `a_refused_run_writes_nothing` (`:1350`), which is already the file's idiom for refusals. That
  converts a guard that never fired into one that fires on every run.

### 1.4 Guard placement and message precedence

Placing the guards after the §3.4a block (≈`cli.rs:635`) is right: the uBAM-specific messages for
`--clock`/`--implicon` + `--output-format ubam` (`:608`, `:614`) keep winning, which is the better
diagnosis. Two knock-on precedence shifts, both benign and worth one line in §3's edge-case table:

- `--clock --fastqc r1.fq` (odd count) now reports the FastQC refusal instead of
  `validate_paired_input`'s "requires an even number of input files".
- `--hardtrim5 20 --clumpify --fastqc a.fq` now reports FastQC instead of `:737`'s clumpify message.

### 1.5 O1 — the assumption is the weaker call, and a wording fix beats relocation

The plan accepts that `--paired --fastqc single.fastq` reports the FastQC refusal before the
"only legal if that file is a uBAM" error. The regression is real and lands on what is very likely the
**more common** mistake — forgetting R2 is more frequent than feeding an interleaved uBAM — and the
proposed message ("…when paired FASTQ output is produced from a single interleaved uBAM…") talks about
a file format the user never mentioned. Today's message is exactly right for them (verified above).

I agree with rejecting the §3.4b relocation (it splits the fix across two files for a message
ordering). The cheaper fix is to word arm 5's message so it is true and useful under **both**
readings, e.g.:

> `--fastqc is not supported with --paired and a single input file: that shape means one interleaved
> uBAM producing FASTQ output, which writes no QC report. Add --output-format ubam (which does produce
> one), pass R1 and R2 as two files, or run FastQC separately.`

One sentence plus remedies, covers both users, no extra file, no ordering change.

### 1.6 The disjunction (§2.3) — confirmed, with the caveat above

`cli.rs:337` carries no clap `requires`, so `cli.fastqc` is false on an args-only run; the
`--fastqc_args`-only SE probe above proves FastQC still runs. The tripwire independently derives the
same predicate (`fastqc_active`, `:223`) and documents the reason. Correct — for `main.rs`. Not for
`clump_only.rs` (§1.2).

### 1.7 Nothing in the repo currently promises `--fastqc` on the five arms

I swept `README.md`, `docs/src/content/docs/**`, `.github/workflows/**`, `tests/**`, `scripts/`,
`examples/`, `Dockerfile`, and the clap doc comments. No text promises `--fastqc` on
hardtrim/clock/implicon/arm-5, and **no existing test or CI step pairs a FastQC flag with any of
them** — so nothing in-repo breaks. But the plan's docs scope (one sentence) is short by four:

| File:line | Text | Why it needs a touch |
|---|---|---|
| `docs/.../guide/outputs.md:145` | "PE FASTQ produces two reports (one per mate)" | the plan's one item; correct to fix |
| `docs/.../guide/outputs.md:143` | "Works on both FASTQ output (default) and uBAM output" | now false for the interleaved-uBAM → FASTQ shape |
| `docs/.../quickstart.md:53` | "Works on both FASTQ and uBAM output paths." | same claim, higher-traffic page |
| `docs/.../guide/flags.md` §"Cross-flag interactions" | — | the canonical home for refusals; has no FastQC entry |
| `docs/.../modes/clump-only.md:77` | "`--fastqc` — runs on the reordered output" | where §1.2's args-only gap lives |

`docs/.../index.mdx:65` ("`--fastqc` is the only flag required") is marketing copy — optional.
The specialty pages (`modes/hardtrim.md`, `clock.md`, `implicon.md`) mention FastQC nowhere, which is
where a user who hits the new error will look — one line each, optional.

**Help text (step 5):** the plan edits `--fastqc` (`cli.rs:327-331`) but not `--fastqc_args`
(`:333-338`), whose "Implies --fastqc" is both the load-bearing premise of the refusal and currently
false under `--clump_only`. And `cli.rs:196` lists `--fastqc` as composing with `--clump_only`.

### 1.8 Integration-test home (step 4) — wrong file suggested

`tests/integration_no_args_help.rs` is narrowly about the bare-invocation help path (its whole module
doc is about routing clap's help to stdout with exit 0); it is not a general home for CLI refusals.
**`tests/integration_paired_format_guard.rs` is the right file**: it already owns `--paired`
single-file refusals, has a `run(args) -> (bool, String)` helper and a `nonexistent_out()` idiom for
asserting *no side effects*, and its module-doc property 2 is literally V7's concern —
"`--clump_only` on the FASTQ-output path keeps its own aux-tag diagnosis… where the first draft of the
fix regressed". Arm 5 and V7 belong there. Only the hardtrim case has no natural home.

---

## 2. Assumptions

| | Plan's status | Reviewer A |
|---|---|---|
| A1 | Verified | **Verified for reachability** (all nine dispatch arms traced; five confirmed empirically). Incomplete as a statement about the *defect class* — see the clump_only args-only gap. |
| A2 | Verified | **Verified** from `main.rs:845/853` + `validate_paired_input:536-544` + probes. |
| A3 | Verified | **Verified** three ways: `main.rs:2006/2521`, `ci.yml:1173`, tripwire `:1048` vs `:1021`; and reproduced. |
| A4 | Verified | **FALSE as stated.** The disjunction is the gate at the 5 `main.rs` gate expressions (7 calls); the 4 `clump_only.rs` gates (5 calls) read `if fastqc`. The clap half ("no `requires`") is correct. |
| A5 | Verified | **Verified for `--fastqc`. False for `--fastqc_args`-only** — the four drivers accept it and do nothing. |
| A6 | Verified | **Verified** (all 8 marked cases pass no FastQC flag; the branch has no self-test). Deletion is safe; the *inventory* is short by 7 call sites and the guard is unreplaced. |
| A7 | Decided by Felix | Not re-argued. One factual note for the CHANGELOG framing: nf-core's `trimgalore` module takes `--fastqc` only from `task.ext.args` (it is not in the module command), and nf-core/rnaseq runs a standalone `FASTQC` module rather than the flag, so the flagship consumer does not break. Residual risk is any wrapper that sets `--fastqc` globally while exposing specialty modes through pass-through args — that wrapper goes from silently QC-less to hard-failing. |

**Unstated assumption worth surfacing:** that a *hard error* is the right severity for a flag whose
absence costs the user nothing but a missing report. The plan's remaining-risks §2 acknowledges the
consequence; the CHANGELOG entry should carry a `**Breaking:**` prefix under the house convention
(one sentence naming what broke and what replaces it), since previously-succeeding invocations now
exit non-zero.

---

## 3. Efficiency

Nothing to analyse and the plan says so correctly: five boolean tests plus one `Option::is_some()`
before any I/O, on the error path only for allocation. If the `fastqc_requested()` accessor route is
taken it is one inlined method call. No memory or scalability dimension.

---

## 4. Validation sufficiency

The ten V-rows cover the failure mode the plan fears most (over-rejection) well: V3 is exactly the
test that catches a dropped `output_format` term, V5 pins the capable arms, V7 pins the `!clump_only`
exclusion. Gaps, in order of how much they matter:

1. **No row exercises arm 5's refusal end-to-end.** V6 drives the binary for `--hardtrim5` only —
   the arm with the *simplest* condition. Arm 5 has the four-term condition and the only real
   look-alike, so it is the one that most needs a binary-level case with the real
   `test_files/ubam_paired_test.bam`. Add it (in `integration_paired_format_guard.rs`).
2. **No row asserts a refused run writes nothing.** The house contract for refusals is
   `a_refused_run_writes_nothing` (tripwire `:1350`) and the "no side effects" property in
   `integration_paired_format_guard.rs`. A CLI-layer refusal trivially satisfies it — assert it
   anyway, because that is what makes the refusal's *position* (before `ensure_output_dir`) a pinned
   fact rather than an accident.
3. **Nothing covers the clump_only args-only case** (§1.2) in either direction — neither a fix nor a
   documented gap.
4. **No row replaces the deleted anti-decay marker** (§1.3). V8 only asserts the file is still green,
   which it would be if the deletion removed real coverage.
5. **No row sweeps for existing consumers** of the now-refused combinations. I ran it and the answer
   is "none" (`tests/`, `.github/workflows/`, `docs/`, `scripts/`, `examples/`, `README.md`) — record
   it in the plan as a completed check so the next reader does not have to redo it.
6. **V10 is not runnable as written.** Content-wise it is fine — all four steps use committed fixtures
   (`illumina_10K.fastq.gz`, `ubam_test.bam`, `ubam_paired_test.bam`, `BS-seq_10K_R{1,2}.fastq.gz`)
   and need only `unzip`. But they `rm -rf` and write into `/tmp/rust_fastqc`, `/tmp/ubam_out_*`,
   which the project sandbox does not permit (`$TMPDIR` and `/tmp/claude` only, and `rm -rf` is
   blocked); they also need a `--release` build at the change's commit. Restate V10 as "re-run the
   four steps with `-o` retargeted to `$TMPDIR/<fresh>`", or state plainly that the real gate is CI
   on the PR — an unrunnable row that everyone skips is worse than an honest "CI covers this".
   Related: V9 lists `cargo nextest run`; if nextest is not installed that row degrades silently.
7. V8's "53 cases pass" is **correct today** (53 `#[test]` in the file) and the deletion removes no
   test, so the number holds — unless recommendation 4 is taken, which changes it.

**Where could this silently produce a wrong result?** Only one place: a guard that is subtly
*narrower* than intended still leaves a silent no-op, and no proposed row detects that in general —
V2 checks the five arms the plan already knows about. The generalisable check is the one the tripwire
was reaching for: enumerate the arms, run each with a FastQC flag, and require every arm to either
produce artefacts or exit non-zero. Recommendation 4 is that check, restricted to the five known
arms; extending it to *all* arms (the capable ones already have CI/tripwire cases) would close the
class rather than the instance, and would have caught the clump_only gap.

---

## 5. Alternatives

| Alternative | Trade-off | Recommendation |
|---|---|---|
| **`Cli::fastqc_requested()` accessor + fix the five `main.rs` clump_only pass-sites** | +1 method, −9 copies of one predicate; makes `--fastqc_args`' "Implies --fastqc" true everywhere; closes §1.2 in the same PR | **Take it.** Smallest change that stops A4 from being wrong. |
| Fix clump_only in a separate issue | Keeps this PR to one behaviour change; but leaves the plan asserting A4/A5 as verified when they are not, and ships a `--fastqc_args`-shaped hole one PR after refusing `--fastqc_args` elsewhere | Acceptable **only** if §2.1/A4/A5 are corrected and the follow-up issue is filed in the same breath |
| Reword arm 5's message (§1.5) vs relocating the check to `main.rs` post-detection | Wording: zero extra files, covers both user populations. Relocation: perfect ordering, two files, splits the family | **Reword.** Keep O1's decision, fix its symptom. |
| Table-driven guards `[(cond, flag, remedy)]` instead of five flat `if`s | Marginally shorter; departs from the §3.4a precedent the file uses eight times | Plan's flat-sequence choice is right — no change. |
| Warn-and-continue / wire FastQC into the five arms | Settled by Felix (A7); no new evidence against it surfaced — the nf-core check in §2 argues *for* refusal being safe | No change. |

---

## 6. Action items

### Critical

1. **Correct §2.1, A4 and A5, and decide the clump_only args-only case.** `clump_only.rs:373/553/844/1094`
   gate on `if fastqc`, not the disjunction; `--clump_only --fastqc_args "…"` is a sixth silent no-op
   (reproduced). Also add the missing 12th site, `clump_only.rs:555`. Preferred resolution: introduce
   `Cli::fastqc_requested()` and use it both for the new guards and at `main.rs:636/671/715/778/819`;
   add a V-row asserting `--clump_only --fastqc_args "--quiet" a.fq` now produces
   `*_clumped_fastqc.{html,zip}`. If it is deferred instead, file the follow-up and delete the
   "Verified" status from A4/A5 — as written the plan certifies a fact that is false.

### Important

2. **Complete the tripwire deletion inventory.** `fastqc_capable` is also the 6th positional parameter
   of `compare()`, called directly at `:1104-1111`, `:1120`, `:1182`, `:1197`, `:1210-1217`, `:1255`,
   plus the pass-through at `:311`. Name all seven in step 3.
3. **Replace the deleted marker with five active refusal cases** (non-zero exit + message + nothing
   written), following `a_refused_run_writes_nothing` (`:1350`). Deleting a latent guard is fine;
   deleting it without leaving an assertion behind removes the only in-repo record that those arms
   cannot run FastQC, in a file whose stated purpose is anti-decay.
4. **Reword arm 5's message** so it is correct for the forgot-R2 user as well as the interleaved-uBAM
   user (§1.5 has draft wording). Today's message for that shape is precisely right and the plan
   replaces it with one that names a format the user never mentioned.
5. **Move the integration case to `tests/integration_paired_format_guard.rs`**, not
   `integration_no_args_help.rs` (§1.8), and put arm 5 + V7 there — the file already owns `--paired`
   single-file refusals, the `--clump_only` message-precedence property, and a no-side-effects idiom.
6. **Widen the docs scope** to `quickstart.md:53`, `guide/outputs.md:143`, and a `guide/flags.md`
   "Cross-flag interactions" bullet (§1.7); and extend step 5 to `--fastqc_args`' help text
   (`cli.rs:333-338`) plus `--clump_only`'s (`cli.rs:196`) if item 1 is deferred.
7. **Make V10 runnable or honest** (§4.6): retarget the four CI steps' `-o` to `$TMPDIR/<fresh>`, or
   state that CI on the PR is the gate. Add the end-to-end arm-5 refusal row and the
   "refused run writes nothing" row.

### Optional

8. Record the "no existing consumer" sweep in the plan (`tests/`, CI, `docs/`, `scripts/`,
   `examples/`, `README.md` — verified clean), so it is not re-derived.
9. Note the two benign precedence shifts in §3's edge-case table: `--clock --fastqc r1.fq` (odd
   count) and `--hardtrim5 --clumpify --fastqc` now report the FastQC refusal first.
10. CHANGELOG: use the `**Breaking:**` prefix — previously-succeeding invocations now exit non-zero.
11. One line on `modes/hardtrim.md`, `modes/clock.md`, `modes/implicon.md` — those pages never mention
    FastQC, and they are where a user who hits the new error will look.
12. While `compare()` is being edited, consider a `comparison_*` self-test for the surviving
    `active && allowed.is_empty()` branch; it has no standing proof either, and the file's stated
    standard (`comparison_detects_a_modified_input`, `:1128`) is that assertions must be shown to fail.
13. Out of scope but adjacent, worth an issue: `--fastqc` with `--demux` reports on the *pre*-demux
    trimmed file (`main.rs:1458` runs before `:1468`), and PE `--retain_unpaired` produces no report
    for the `*_unpaired_*` files. Neither is a silent no-op (artefacts do appear), so neither belongs
    in this fix — but both are the same "which outputs get QC'd?" question, undocumented.
