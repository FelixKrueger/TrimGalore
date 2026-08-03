# CODE review — Reviewer B — `-a2` honoured with presets (#369)

**Target:** working tree on `fix/a2-honoured-with-presets`, based on `dev` @ `4276810`. Nothing committed.
**Reviewed:** `src/main.rs`, `src/cli.rs`, `tests/integration_adapter2.rs` (new), `.github/workflows/ci.yml`, `CHANGELOG.md`, `docs/src/content/docs/guide/flags.md`, plus `PLAN.md` §13 as a claim set to verify.
**Read-only:** no repo file was modified by me except this report. A throwaway git worktree at `4276810` was built for pre/post comparison and removed again.

## ⚠ The tree changed under me during this review

When I started, `git status` showed five modified files (`ci.yml`, `CHANGELOG.md`, `flags.md`, `cli.rs`, `main.rs`) plus the untracked test file. Partway through, `src/adapter.rs` appeared as modified and `src/cli.rs` / `src/main.rs` changed content — someone else edited the tree concurrently, applying fixes for three of my findings while I was writing them up.

The report below is organised accordingly: findings are stated as I verified them, each annotated with its status **as of my last check of the live tree**. Every "already fixed" annotation was re-verified by running the rebuilt binary, not inferred from the diff.

One consequence of the concurrent edit is a **current CI-breaking failure** (F0 below).

## Verdict

The design is right and the load-bearing claims hold. I tried to refute the two that matter most — the suppression sentinel and "no accepted path moved" — and could not; the second I verified by measurement, not by reading.

- **F0 (Critical, live):** `cargo fmt --all -- --check` **fails right now** on `src/adapter.rs:375`, from the concurrent edit. CI's lint job gates on this.
- **B1, B2, B3 (High/Medium):** found and verified against the tree as I read it; **since fixed** in the live tree, and I confirmed the fixes work.
- **B4-B8 (Medium/Low):** still open as of my last check — documentation accuracy, `--help` hedging, test-coverage gaps, comment length, register.

---

## What I verified, and how

| Claim | Method | Result |
|---|---|---|
| `cargo test --release` | ran, on both tree states | 454 pass, 0 fail, both times |
| `cargo clippy --all-targets --release -- -D warnings` | ran, on both tree states | clean both times |
| `cargo fmt --all -- --check` | ran, on both tree states | clean before the concurrent edit; **fails after** (F0) |
| §13 dev-1: suppression is the only source of an empty R1 adapter | read every producer + ran the clap conflict cases | **confirmed** (§A) |
| No accepted path moved | **built `dev` @ `4276810`, byte-compared 7 invocations** | **confirmed** (§B) |
| Read 1 guarantee is structural | read the diff | confirmed (§C) |
| `--clump_only` error still wins over the new warning | read line order + ran it | confirmed (§D) |
| No duplicated error message for malformed `-a2` | ran it | confirmed (§D) |
| CI step behaves as claimed **in bash** | reproduced the loop under `bash -c` with an arg printer | confirmed (§E) |
| CI grep cannot pass vacuously | inspected a real R2 report | confirmed (§E) |
| Perl 0.6.11 accepts `--illumina -a2 SEQ` and trims R2 with it | read `git show 0.6.11:trim_galore` | confirmed (§E) |
| Perl md5 half of the new CI step | **not run** — no conda/cutadapt here | unverified (§E) |

---

## F0 — [Critical, currently live] `cargo fmt --check` fails

Introduced by the concurrent edit, not by the original implementation.

```
$ cargo fmt --all -- --check
Diff in /Users/fkrueger/Github/TrimGalore/src/adapter.rs:375:
-fn parse_adapter_specs_inner(
-    specs: &[String],
-    announce: bool,
-) -> Result<Vec<(String, String)>> {
+fn parse_adapter_specs_inner(specs: &[String], announce: bool) -> Result<Vec<(String, String)>> {
```

The signature fits on one line, so rustfmt wants it there. CI's `lint` job runs `cargo fmt --all -- --check` and fails on any diff. `cargo fmt --all` fixes it. (Verified by running the check twice: clean at the tree state I first read, failing now.)

---

## Findings verified against the tree, with current status

### B1 — [High] The override NOTE was not gated on `cli.paired`; single-end output contradicted itself — **since fixed**

As I first read it, `src/main.rs` printed the NOTE unconditionally while gating the `Adapter 2 (Read 2):` line on `cli.paired`, and `apply_adapter2_override` ran regardless of `cli.paired`. Reproduced by running the binary:

```
$ trim_galore --small_rna -a2 AAATCAAAAAAAC test_files/BS-seq_10K_R1.fastq.gz -o …
WARNING: -a2/--adapter2 was given but is not used in this mode (… single-end run). Ignoring.
Adapter: smallRNA (TGGAATTCTCGG)
NOTE: Read 2 adapter taken from -a2; the smallRNA default (GATCGTCGGACT) is not used.
```

Two statements in one run, one of which must be wrong to the reader. This is the exact self-contradiction V20 was written to close, and the V20 test cannot see it: `single_end_does_not_announce_a_read2_adapter` (`tests/integration_adapter2.rs:410-418`) uses `--illumina`, which has no R2 default, so `displaced_r2` is `None` and no NOTE can fire. The contradiction needs a preset *with* an R2 default — the very row §3.2's truth table marks as the NOTE case.

I confirmed the mutation was otherwise harmless in single-end mode: `run_single_file` never receives `adapters_r2`, so SE text and JSON reports were unaffected (checked the JSON: `"adapters_r2": []`); `trim_read` only consults `config.adapters_r2` when `is_r2` (`src/trimmer.rs:112`), and `r2_adapter_count` is called only from the paired functions. Stderr-only — but on the most-read surface.

**Status: fixed in the live tree** with `if cli.adapter2.is_empty() || !cli.paired { return Ok(None); }` plus the one-line comment *"Not paired: `Cli::validate` has already warned, and Read 2 does not exist."* Re-ran the reproduction: only the WARNING now appears, no NOTE.

**Still worth doing:** extend `single_end_does_not_announce_a_read2_adapter` to a `--small_rna` case asserting `!stderr.contains("NOTE: Read 2 adapter taken from -a2")`. Nothing in the suite currently guards the fix.

### B2 — [Medium] Malformed `-a2` was fatal in exactly the modes that declare it unused — **since fixed**

As I first read it, `src/cli.rs` parsed `-a2` *before* computing `unusable_reason`. All four of these ran to completion on `dev` @ `4276810` and exited 1 (verified by running both binaries):

| Invocation | pre-fix | as first reviewed |
|---|---|---|
| `--paired --clock -a2 ZZZQQQ …` | exit 0 | exit 1 — `Adapter sequence must contain only DNA characters …` |
| `--paired --clock -a2 file:missing.fa …` | exit 0 | exit 1 — `Cannot open adapter FASTA file: missing.fa` |
| `--paired --hardtrim5 30 -a2 '' …` | exit 0 | exit 1 — `Empty adapter sequence` |
| `--illumina -a2 ZZZQQQ single.fq.gz` | exit 0 | exit 1 |

PLAN §3.6 and §7 item 5 scope the new non-zero exit to "preset and auto-detect paths"; the CHANGELOG says only that `-a2` given where it cannot apply "now warns instead of being silently dropped". Neither covered "and hard-fails if the value it will not use is malformed". The `--clock`/`--hardtrim` rows were the awkward ones: those modes never look at an adapter, so the run failed on a flag the tool was about to declare irrelevant. A pipeline reusing one flag template across paired and single-end samples with `-a2 file:r2_adapters.fa` would break on the single-end sample if that FASTA is absent.

A related cosmetic consequence of the same ordering: `--paired --clock -a2 'A{10}'` printed the expansion notice and *then* the "is not used in this mode" warning.

**Status: fixed in the live tree** — `unusable_reason` is now computed first and the parse only runs in the `None` arm (`src/cli.rs:941-960`). Re-verified: `--paired --clock -a2 ZZZQQQ` now exits 0 with only the WARNING. `--paired --illumina -a2 ZZZQQQ` still exits 1, so V17 and `malformed_a2_is_rejected_on_preset_paths` are unaffected (suite still green).

Two notes on the fix as applied: (a) its comment is **three lines** (`src/cli.rs:938-940`), over CLAUDE.md's "one line, two maximum" limit — the third line, *"A value the run is about to ignore is not worth failing on"*, is the rationale and belongs in the commit message; (b) no test covers the new lenient behaviour, so nothing stops it regressing.

### B3 — [Medium] `-a2 A{N}` announced its expansion twice — **since fixed**

`-a2` was parsed in `Cli::validate` and again in `apply_adapter2_override`, and `parse_adapter_spec` has an `eprintln!` side effect. Measured on both binaries:

| Invocation | pre-fix | as first reviewed |
|---|---|---|
| `--paired -a SEQ -a2 'A{10}' R1 R2` | 1 message | 2 |
| same, 2 pairs | 2 | 3 |

i.e. `1 + pairs`. The `-a` + `-a2` path is the one PLAN §3.5 / V4 calls "unchanged" — output files and both report formats *are* unchanged (verified byte-identical, §B), but stderr was not, and that path was already correct before the fix. The double call also meant a `-a2 file:adapters.fa` FASTA was opened once in `validate` plus once per input file.

**Status: fixed in the live tree** via `parse_adapter_spec_quiet` / `parse_adapter_specs_quiet` in `src/adapter.rs`, called from `Cli::validate`. Re-verified: the count is back to 1. This is the edit that introduced F0.

**Still worth doing:** the `file:` FASTA is still read twice per run when the mode is usable (`validate` + once per input file). Harmless, but §6's per-file re-read note now undercounts by one.

### B4 — [Medium, open] Gating the display on `cli.paired` also changes runs that pass no `-a2` at all

`src/main.rs:837`. `resolve_adapter`'s `--small_rna` / `--bgiseq` branches return `to_r2_vec()` regardless of `cli.paired`, and the old display was `if !adapters_r2.is_empty()`. So single-end preset runs printed the line before and do not now. Verified against the pre-fix binary:

```
SE --small_rna : pre-fix "Adapter 2 (Read 2): GATCGTCGGACT" → post-fix absent
SE --bgiseq    : pre-fix 1 occurrence                       → post-fix 0
```

The new behaviour is better. But PLAN §3.5 and §7 both assert "No change for … any invocation without `-a2`", and the CHANGELOG files the SE display change under the `-a2`-given bullet ("single-end runs no longer print an `Adapter 2 (Read 2):` line for an adapter they will not use"). Documentation-only: say that single-end `--small_rna`/`--bgiseq` stop printing it with **or without** `-a2`.

---

## Verified-clean claims worth recording

### §A — The suppression sentinel: I could not refute it

`src/main.rs` infers `--consider_already_trimmed` suppression from `adapters_r1.len() == 1 && adapters_r1[0].1.is_empty()`. §13 deviation 1 claims that is the only way an empty R1 sequence reaches that point. **Confirmed** — every producer bails on empty:

- single sequence: `anyhow::bail!("Empty adapter sequence")` in `parse_adapter_spec_inner`, and `A{0}` bails immediately after;
- embedded `" -a "` form: empty parts are skipped and an all-empty spec bails with `"No valid adapter sequences found in multi-adapter specification"`;
- `file:` FASTA: empty per-entry sequences bail (both the flush-previous and flush-last branches), and an entry-less file bails with `"No adapter sequences found in FASTA file"`;
- all five preset consts have a non-empty `seq` (`src/adapter.rs:36-64`);
- the synthetic suppressed preset (`src/adapter.rs:219-227`) is the sole `seq: ""`.

Two reachability facts tighten it further, both verified by running the binary:

- `--consider_already_trimmed` is `conflicts_with_all` the five presets (`src/cli.rs:292-293`), so `--illumina --consider_already_trimmed 10000` is rejected by clap before any of this — suppression can only arise on the auto-detect branch. (Confirmed: `error: the argument '--illumina' cannot be used with '--consider_already_trimmed'`.)
- `-a` is **not** in that conflict list, but the `-a` branch returns before detection runs, so its R1 sequence is always a validated non-empty one.

Caveat worth knowing: `validate_adapter_sequence` would *accept* an empty string — `bytes().all()` is vacuously true — so the guarantee rests entirely on the callers' explicit emptiness bails, not on the validator. That is a sentinel coupling with nothing enforcing it. A `debug_assert!`, or a pointer to `adapter.rs:219` in the existing one-line comment, is cheap insurance. Not a defect today.

### §B — "No accepted path moved": verified by measurement, not reading

Built `dev` @ `4276810` in a worktree and ran both binaries into the *same* output path (sequentially) so no path string could differ. Compared decompressed `.fq.gz` and the `.txt` reports:

| Invocation | FASTQ + `.txt` report |
|---|---|
| `--paired -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC` | identical |
| `--paired` (auto-detect, no `-a2`) | identical |
| `--paired --small_rna` | identical |
| `--paired --bgiseq` | identical |
| `--paired --illumina` | identical |
| single-end default | identical |
| `--paired --rrbs` | identical |

`*_trimming_report.json` differed in exactly one field, `"command_line"`, which embeds `argv[0]` — an artefact of the two binaries living in different directories, not a behaviour change (diff inspected; it is the only differing line). A5 and V13 hold, and the Perl-parity md5 matrix is safe: no invocation without `-a2` produces different bytes. This also confirms B3 was stderr-only.

### §C — The Read 1 guarantee is structural

Read the whole diff of `resolve_adapter`. The only change inside it is in the `-a` branch: `parse_adapter_specs(&cli.adapter2)` deleted and the R2 slot of the returned tuple changed from that value to `Vec::new()`. No R1 expression, no branch condition and no branch ordering was touched; the other six branches are byte-identical. `adapters_r1` is passed to `apply_adapter2_override` by shared reference, so it cannot be mutated there. Sound.

### §D — Double validation and `--clump_only` ordering are both clean

- **Ordering.** The `--clump_only` bail for `-a2` is at `src/cli.rs:716`, inside the `if self.clump_only` block opened at `:695`; the new block is at `:938`. Ran it: `--clump_only --cores 2 -a2 AAATCAAAAAAAC …` prints only `Error: --clump_only does not trim; -a2/--adapter2 is not compatible` — no warning first.
- **No double error.** `Cli::validate` is called once (`src/main.rs:166`) and bails before the override's parse can run. Ran `--paired --illumina -a2 ZZZQQQ`: the DNA-characters error appears exactly once, exit 1. The two parses only both execute on *valid* input — that was B3, not an error-path problem.
- Minor layering note: `Cli::validate` now does filesystem I/O (`file:` FASTA). It already does `path.exists()` checks a few lines earlier, so not new in kind.

### §E — The new CI step

**Bash semantics — correct.** Reproduced the loop with an argument printer:

```
$ bash -c 'for mode in "--illumina" ""; do printargs --paired $mode -a2 AAATCAAAAAAAC …; done'
argc=8  [1]=<--paired> [2]=<--illumina> [3]=<-a2> [4]=<AAATCAAAAAAAC> …
argc=7  [1]=<--paired> [2]=<-a2> [3]=<AAATCAAAAAAAC> …
```

The empty `$mode` iteration elides the word entirely, which is the intended auto-detect invocation. (The zsh/bash divergence from iteration-log #1 only bites multi-word values; neither value here is multi-word, so both shells agree.)

**Vacuous-pass analysis.** The job sets `defaults.run.shell: bash -l {0}` (`.github/workflows/ci.yml:249-251`), which *replaces* GitHub's default `bash -e {0}` — no `errexit`, no `pipefail`, so a failing `conda run` would not abort the step. If both binaries failed, both md5s would be the md5 of an empty stream and the comparison would pass. The step is saved by its last command: `grep -q` on a missing R2 report exits 2 and triggers the `||` block. So it cannot pass vacuously. (Several pre-existing steps in this job have the same structure *without* a trailing grep — out of scope, but worth knowing.)

**The grep target is right and discriminating.** Inspected a real R2 report from the fixed binary. Three lines mention the sequence:

```
 9: Optional adapter 2 sequence (Read 2): 'AAATCAAAAAAAC'
18: Command line parameters: -j 1 -e 0.1 -q 20 -O 1 -a AAATCAAAAAAAC BS-seq_10K_R2.fastq.gz
33: Sequence: AAATCAAAAAAAC; Type: regular 3'; Length: 13; Trimmed: 5444 times.
```

Only line 18 matches `a AAATCAAAAAAAC` (line 9 has `'` before the sequence, line 33 has `: `). The filename matches `naming::report_name`. Pre-fix that line reads `-a AGATCGGAAGAGC`, so §13's V15 result stands. Line 9 is written from `config.adapters_r2`, so it too is absent pre-fix — either line would have worked as the assertion.

**Perl-side premise — confirmed from the tagged source** (`git show 0.6.11:trim_galore`):

- the adapter-vs-preset `die`s are all inside `if (defined $adapter)` (lines 3129-3143) and none tests `$adapter2`, so `--illumina -a2 SEQ` is accepted;
- `-a2` without `--paired` dies (3147-3149) — the documented divergence;
- with `$a2` set, Read 2 is trimmed with `-a $a2` in the same single pass (1355-1358); with `$a2` empty both reads use `-a $adapter` (1362-1365). Exactly our post-fix semantics, so the md5 comparison should hold.

**What I could not verify:** the Perl md5 half itself — no conda/cutadapt in this environment, so that assertion is first exercised on the branch push. Everything checkable about it points the right way; flagged only so "CI step reviewed" is not read as "CI step executed".

### Efficiency

Nothing in a per-read path; one extra `parse_adapter_specs` per input file over a vector that is empty or length 1. §6's note that a `-a2` over 64 bp skips the Myers prefilter is accurate — `src/alignment.rs:19-20` and `:271` both cap at 64 bp — and the fix does make `-a2 T{150}` reachable from the preset and auto-detect paths for the first time. Inherent to the fix; correctly documented; no action.

### Errors

No unhandled conditions found. `apply_adapter2_override` propagates its parse error with `?`; `displaced` is captured before the replacement, so the NOTE cannot report the new value as the displaced one; `adapters_r2.first()` cannot lose information because every preset R2 default is single-entry.

---

## Remaining open findings

### B5 — [Low] `--help` was not hedged, only the docs page

`docs/.../guide/flags.md:12` gained "unless the user provides their own `--adapter2` value" — good, and it copies the `--clip_R2` phrasing on the next line as §4 Step 8 intended. The same claim is still unconditional in three `--help` strings, the more-consulted surface:

- `src/cli.rs:58` — `-a2`: "Auto-set by --small_rna and --bgiseq presets."
- `src/cli.rs:73` — `--small_rna`: "…and sets --adapter2 (GATCGTCGGACT, Illumina small RNA 5')."
- `src/cli.rs:82` — `--bgiseq`: "Sets --adapter2 for Read 2."

Suggest the same hedge on at least `:73` and `:82`.

### B6 — [Low] Test-coverage gaps against the plan's own V-rows

`tests/integration_adapter2.rs`. None of these is a vacuous pass — I checked each assertion for the `split_whitespace` trap that bit `suppression_wins_over_a2`, and the rest read report lines or stderr substrings that cannot collapse. These are missing guards, not false ones:

- `:369-405` `unusable_a2_warns_with_a_mode_accurate_reason` covers `--clock`, `--hardtrim5`, single-end. `--hardtrim3` and `--implicon` each have their own arm in the `src/cli.rs` cascade and are untested; §13 says both were checked by hand.
- `:425-443` `malformed_a2_is_rejected_on_preset_paths` asserts only `!ok`, so it would pass if the run failed for an unrelated reason. Add a per-case `stderr.contains(…)`: `"only DNA characters"`, `"Empty adapter sequence"`, `"Cannot open adapter FASTA"`.
- `:448-467` `repeated_a2_applies_on_a_preset_path` asserts the stderr display line, whereas V22's claim is about the R2 report's per-adapter breakdown and `r2_adapter_count`. The display reads the same vector that becomes `TrimConfig.adapters_r2`, so it is a fair proxy, but it does not prove the stats array was sized from `adapters_r2.len()`.
- `:213-257` and `:284-295` discard `run()`'s success flag. They fail via a helper panic rather than silently, but the message would be misleading.
- No positive test for `-a2 file:adapters.fa` on a preset path, which §3.1/A3 lists alongside `A{N}` and repeats. Only the missing-file negative is covered.
- Nothing guards the B1 or B2 fixes (see those entries).

### B7 — [Low] Comment length

The original inline comments comply with CLAUDE.md's one-line rule (`src/main.rs` override call site, the paired-display note, the sentinel note, and the two-line `#369` note in `src/cli.rs` — one wrapped sentence with a stated reason).

Two deviations:

- the doc comment on `apply_adapter2_override` spends its last two lines arguing *why* suppression wins ("trimming R2 while R1 is left alone would be asymmetric, and the mode announces that only quality trimming will happen") — that reasoning is already in §3.6 and belongs in the commit message;
- the concurrent B2 fix's comment is three lines (`src/cli.rs:938-940`), one over the limit.

For the record, the new test file is **not** over-commented by this repo's standards: 60 comment lines / 517 total = 11%, against 12% for `tests/integration_clump_only.rs` and 16% for `tests/integration_ubam.rs`.

### B8 — [Low] Register nit on the suppression warning

`src/main.rs`: "WARNING: -a2/--adapter2 **not applied** — adapter trimming is suppressed for this library (--consider_already_trimmed). **Ignoring.**" says the same thing at both ends. Either drop the trailing "Ignoring." or open with the house-standard "was given but is not used…". The single-end/specialty warning reads well and matches the existing deprecation-warning pattern; em-dash usage is consistent with house style.

---

## Recommendations, prioritised

**Critical**
1. **F0** — run `cargo fmt --all`. `src/adapter.rs:375` currently fails `cargo fmt --check`, which gates CI's lint job.

**High**
2. **B1 follow-through** — the fix is in; add the guard. Extend `single_end_does_not_announce_a_read2_adapter` to `--small_rna` and assert the NOTE is absent (`tests/integration_adapter2.rs:410-418`).

**Medium**
3. **B2 follow-through** — the fix is in; add a test for the new lenient path (`--paired --clock -a2 ZZZQQQ` exits 0 with the warning) and trim the three-line comment at `src/cli.rs:938-940` to two.
4. **B4** — CHANGELOG/plan wording: single-end `--small_rna`/`--bgiseq` stop printing `Adapter 2 (Read 2):` with **or without** `-a2`, so "no change for any invocation without `-a2`" needs qualifying.
5. Re-run `cargo test --release` and `cargo clippy --all-targets --release -- -D warnings` after the fmt fix. Both were clean at my last check, but the tree has been edited more than once during this review and my results are only as current as that check.

**Low**
6. **B5** — hedge the `--help` text for `--small_rna` / `--bgiseq` / `-a2` to match the docs page (`src/cli.rs:58`, `:73`, `:82`).
7. **B6** — add `--hardtrim3` and `--implicon` to the unusable-mode test; assert the error message in the malformed-`-a2` test; add a `-a2 file:` positive case.
8. **B7** — move the suppression rationale out of `apply_adapter2_override`'s doc comment into the commit message.
9. **B8** — de-duplicate the suppression warning's phrasing.
10. Optional: `debug_assert` or a source pointer next to the suppression sentinel so the invariant in §A is enforced rather than merely true.
