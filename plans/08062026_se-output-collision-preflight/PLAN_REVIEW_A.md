# Plan Review A — Output-collision pre-flight for the four uncovered dispatch paths (#383)

**Plan:** `plans/08062026_se-output-collision-preflight/PLAN.md` (v1, 2026-08-06)
**Reviewed against:** `dev` @ `72624c7`, binary `target/release/trim_galore` (built 2026-08-06 16:41)
**Reviewer:** A (independent; no shared state with Reviewer B)
**Verdict:** Sound plan, correct diagnosis, accurate anchors. One Critical gap: the fix
does not close #383 for path-aliased spellings of the same directory, which I reproduced
against the current binary. Five Important items, six Optional.

---

## 1. Verification of the plan's factual claims

I checked every load-bearing claim rather than taking it on trust. Unless noted, the claim
is **confirmed**.

### 1.1 Anchors

| Plan claim | Verified |
|---|---|
| `preflight_collision_bam` at `main.rs:64-80`, path-generic `&[PathBuf]` | ✅ exact |
| Doc-comment at `:55-59` documents `resolve_clump_layout`, misattached to the pre-flight | ✅ exact |
| Three call sites `main.rs:533`, `:566`, `:614` | ✅ exact |
| `--paired` FASTQ pre-flight `main.rs:672-727` | ✅ (`if cli.paired {` at 672, candidate loop ends 727) |
| `--paired` uBAM pre-flight `main.rs:1799-1821` | ✅ (block 1800-1820) |
| `run_specialty_paired` pre-flight `main.rs:2410-2430` | ✅ (block 2412-2429) |
| `--clump_only` SE FASTQ pre-flight `main.rs:481-498` | ✅ exact |
| `--clump_only` SE uBAM pre-flight `main.rs:609-614` | ✅ exact |
| Insertion points `:344` (hardtrim5), `:369` (hardtrim3), `:781` (SE `} else {`), `:1857` (`// Single-end loop.`) | ✅ all four exact |
| `hardtrim5` branch returns at `main.rs:367` | ✅ exact |
| `gzip` derived from `cli.input[0]` at `main.rs:295`; `ensure_output_dir` at `:300` | ✅ exact |
| `norm_path` at `io.rs:44-46` | ✅ exact |
| `specialty.rs:435/:456/:471/:504` are the `None => PathBuf::from(filename)` arms | ✅ all four exact |
| `hardtrim_output_name` at `:423`, `hardtrim_bam_output_name` at `:446`, both private | ✅ exact |
| `report_name` / `json_report_name` key on full input filename, `io.rs:401-430` | ✅ exact |
| `strip_fastq_extensions` at `io.rs:444` | ✅ exact |
| `--basename` + >1 SE input rejected at `cli.rs:620` | ✅ exact, and outside the UBam block so it covers both formats |
| `--paired` rejects duplicate pairs / R1==R2, `cli.rs:~535/~549` | ✅ (R1==R2 at 535, duplicate-pair loop 543-556) |
| CI case-alias guard at `ci.yml:617-640`; `test -z "$(ls -A …)"` at `:615`/`:640` after the `-o` runs at `:604`/`:628` | ✅ (plan's `:604`/`:628` name the invocations, not the assertions — harmless off-by-a-few) |
| CHANGELOG has `### Unreleased` → `#### Bug fixes` | ✅ (lines 4/6) |

### 1.2 The hole map is exhaustive, not just illustrative

I enumerated every loop over `cli.input` in `main.rs`: lines 345, 370, 484/499, 558/572,
611/615, 683/737, 784, 1803/1823, 1858, 2415/2432. Exactly four processing loops lack a
preceding pre-flight — 345, 370, 784, 1858 — which is precisely the plan's set. **No fifth
hole exists.** This is the strongest part of the plan.

### 1.3 Measured behaviour (all four reproductions re-run)

Tiny 40-read fixtures with `ALPHA_`/`BETA_` read-ID prefixes, in
`…/scratchpad/review-a/`. Negative control run first: single-input `alpha.fastq.gz` gives
`ALPHA=40 BETA=0`, so the attribution grep is capable of reporting a non-zero ALPHA count.

| Repro | Result |
|---|---|
| (a) `sample.fastq.gz` + `sample.fq.gz` | exit 0, one `sample_trimmed.fq.gz`, `ALPHA=0 BETA=40`, **two** reports |
| (c) `--hardtrim5 20 dirA/same.fastq.gz dirB/same.fastq.gz` | exit 0, one `./same.20bp_5prime.fq.gz` in CWD, `ALPHA=0 BETA=40` |
| (d) `--hardtrim5 20 --output-format ubam` colliding stems | exit 0, "Writing hard-trimmed" ×2, one `.bam` |
| SE trim uBAM colliding stems | exit 0, `Output: sample_trimmed.bam` printed twice, one `.bam` |
| §2.3 accept case (SE trim, two dirs, no `-o`) | exit 0, two outputs, `dirA`=ALPHA only, `dirB`=BETA only — **correct attribution**, not merely two files |
| A6 premise (same file twice, today) | exit 0, one output, 40 reads, correct content |
| Open 4 (`--hardtrim5 20 --hardtrim3 15`) | exit 0, only `x.20bp_5prime.fq.gz` written — hardtrim3 silently ignored |

Every behavioural claim in §2.1/§2.2/§2.3/§8-A3/§10-Open-4 holds.

### 1.4 Downstream claim (validation matrix)

I read every `trim_galore` invocation in `ci.yml` (48 of them). **No step passes more than
one single-end input**, and the multi-pair steps use `A_`/`B_` distinct stems. The plan's
"validation matrix unaffected" is correct.

### 1.5 One factual error

§2.4 says "**five** hand-rolled copies … (`main.rs:482`, `:680`, `:1801`, `:2414`, and the
`--clump_only` SE FASTQ one at `:486`)". `grep -n "Output path collision"` returns five
sites, but one is the helper itself (`:70`). There are **four** hand-rolled copies —
`:482`, `:680`, `:1800`, `:2412` — and `:482`/`:486` are the same block counted twice
(`:482` is `let mut out_paths`, `:486` is the `if let Some(existing)` inside it). This
misstates the size of the out-of-scope consolidation commit (§10-3). Cosmetic, but the
count appears twice in the plan.

### 1.6 Signatures compile as written

`main.rs` imports `trim_galore::io as naming`, `trim_galore::specialty`, `std::path::Path`;
at `:344`/`:369` `output_dir: Option<&Path>` and `gzip: bool` are both in scope, and `cli`
is owned so `&cli` type-checks. §4.1's body uses `norm_path` unqualified (correct inside
`io.rs`), fully-qualified `std::collections::HashMap`, and `anyhow::bail!` — `io.rs`
already imports `anyhow::{Context, Result}` and `std::path::{Path, PathBuf}`, so no new
`use` lines. `preflight_output_collisions` being `pub` in the library means no
`dead_code` warning even though only the binary calls it. §4.3's five-arg helper is under
clippy's `too_many_arguments` threshold. `preflight_collision_bam(&[planned])` at `:533`
coerces `&[PathBuf; 1]` → `&[PathBuf]` unchanged. Nothing here will fail
`clippy -D warnings`.

---

## 2. Logic review

### 2.1 CRITICAL — the fix does not close #383 for path-aliased spellings

`norm_path` is a pure string lowercase. The candidate paths handed to it come from
`single_end_output_name`, which does `input.parent().join(filename)` — so the *directory
portion of the key is whatever the user typed*. Two spellings of the same directory produce
two different keys for the same physical file, and the pre-flight accepts.

Reproduced against the current binary (the fix cannot change these outcomes, because the
candidate list Step 7 builds is exactly the `Output:` line the binary already prints):

```
$ cd <dir>; trim_galore <abs>/sample.fastq.gz sample.fq.gz     # ALPHA + BETA
exit=0
Output:   /private/tmp/.../t2/sample_trimmed.fq.gz
Output:   sample_trimmed.fq.gz
ALPHA=0 BETA=40
$ ls -i <abs>/t2/sample_trimmed.fq.gz sample_trimmed.fq.gz
193763599 …/t2/sample_trimmed.fq.gz
193763599 sample_trimmed.fq.gz          ← same inode
```

```
$ trim_galore ./sample.fastq.gz sample.fq.gz
exit=0
Output:   ./sample_trimmed.fq.gz
Output:   sample_trimmed.fq.gz
ALPHA=0 BETA=40                          ← same inode, 193763651
```

Lowercased, `"/private/…/t2/sample_trimmed.fq.gz"` ≠ `"sample_trimmed.fq.gz"` and
`"./sample_trimmed.fq.gz"` ≠ `"sample_trimmed.fq.gz"`, so `seen.insert` returns `None`
twice and the helper returns `Ok(())`. **After the fix, #383 reproduces verbatim — same
exit 0, same 40 reads lost, same two reports — if either argument carries a `./` prefix or
an absolute path while the other does not.**

This is not exotic. `find`/`xargs` lists, `ls -d ./*` idioms, Snakemake and Nextflow input
staging, and hand-assembled file lists all mix absolute and relative spellings within one
argument vector. A uniform glob (`*.fastq.gz` or `./*.fastq.gz`) is safe because the prefix
is consistent; a *mixed* list is not.

§1 claims the change makes it "impossible for one invocation to silently overwrite one
input's output with another's". That claim is false as written, and no assumption in §8
records the limitation.

**Fix — cheap and strictly safe.** Key on the lexically absolutised path. I verified
`std::path::absolute` (stable, no filesystem access, well inside the 1.88 floor):

```
sample_trimmed.fq.gz       -> /…/review-a/sample_trimmed.fq.gz
./sample_trimmed.fq.gz     -> /…/review-a/sample_trimmed.fq.gz   ← `.` removed
a//sample_trimmed.fq.gz    -> /…/review-a/a/sample_trimmed.fq.gz  ← separators collapsed
a/../sample_trimmed.fq.gz  -> /…/review-a/a/../sample_trimmed.fq.gz  ← `..` kept (POSIX-correct)
```

So absolutising the **key only** (keep the message printing the user's original paths, so
the CI greps and the two `test -z "$(ls -A …)"` assertions stay valid) closes both aliases
above at the cost of one `getcwd`. It can only ever *increase* rejections, never decrease
them, and two paths that absolutise equal are genuinely the same file — so no new false
positives beyond the case-fold one #216 already accepted. It also strictly improves the
existing paired / clock / implicon / clump_only guards once the consolidation commit lands.
`..` and symlinks still alias; that residual should become an explicit assumption rather
than going unstated (the plan already argues correctly against `fs::canonicalize`, and that
argument is untouched — `absolute` is lexical).

Note the namers' `parent().join()` already collapses duplicate separators, which is why
`a//x` vs `a/x` is *not* a hole today. It is specifically `.` and absolute-vs-relative that
leak.

### 2.2 Assumption A2 is false in one documented flag combination

A2 says the primary output path is a sufficient collision key for "reports, JSON, FastQC,
demux". I verified the four secondaries individually:

- **reports / JSON** — ✅ A2's reasoning is right, and I checked the implication in both
  directions. `report_name` is injective in `(parent-string, file_name-string)`, so equal
  report names force equal input path strings, hence equal stems, hence equal primaries.
  Conversely distinct primaries in the same directory require distinct stems, which require
  distinct file names, so distinct reports. The claim holds.
- **`--demux`** — ✅ safe, and for a reason the plan does not state: `demux::demultiplex`
  resolves its output directory as `output_dir.unwrap_or(trimmed_file.parent())`
  (`demux.rs:142-147`) and derives `base_name` from the trimmed file's name
  (`demux.rs:129-140`). Both keys are functions of the primary path, so the resolution is
  parallel to the primary's by construction. Worth one sentence in §2.3 — the current
  argument ("`--demux` is handed `&output_path`") is not sufficient on its own, because
  being handed the path says nothing about how the directory is resolved.
- **FastQC, default** — ✅ safe. `fastqc-rust 1.0.1` `runner.rs:224-233` falls back to
  `group.files.first().parent()`, i.e. the primary's parent. Parallel resolution again.
- **FastQC with `--fastqc_args "-o DIR"`** — ❌ **counterexample.**
  `fastqc.rs:91-96` (`apply_fastqc_args`) overwrites `config.output_dir` with an arbitrary
  directory that the pre-flight never sees. Verified:

```
$ trim_galore --fastqc_args "-o <qc>" dirA/same.fastq.gz dirB/same.fastq.gz
exit=0
Output:   dirA/same_trimmed.fq.gz      ← distinct primaries: pre-flight ACCEPTS
Output:   dirB/same_trimmed.fq.gz
$ ls <qc>
same_trimmed_fastqc.html
same_trimmed_fastqc.zip               ← ONE pair of artifacts, not two
$ unzip -p <qc>/same_trimmed_fastqc.zip …/fastqc_data.txt | grep Filename
Filename	same_trimmed.fq.gz          ← dirA's QC report silently clobbered by dirB's
```

Consequence is a lost QC artifact, not lost reads, and the hole is pre-existing and
identical on the `--paired` path (which has had the pre-flight since #216) — so this is not
a regression the plan introduces. But A2 is stated as an unqualified universal and flagged
as "load-bearing", so it needs narrowing: either add the `--fastqc_args`-derived output
directory to the candidate list, or scope A2 to "given FastQC's *default* output-directory
resolution" and record `--fastqc_args -o` as a known residual. V6 as designed cannot detect
this (it only compares `report_name` with `single_end_output_name`).

### 2.3 A3, A4, A5, A6, A7

- **A3** ✅ verified two ways: `main.rs:367` returns, and empirically
  `--hardtrim5 20 --hardtrim3 15` writes only the 5′ output. Per-block checks suffice.
  Also confirmed `5prime`/`3prime` keeps the two modes' names distinct, so even a
  hypothetical combined run could not self-collide.
- **A4** ✅ verified. `gzip` is one run-wide flag (`:295`) threaded into both the candidate
  namer and the writer, so prospective and actual paths cannot disagree. I checked the
  awkward orderings: `.fastq.bgz` first → `gzip=false` → both candidates `sample_trimmed.fq`
  → still collides; `.fastq.gz` first → both `.fq.gz` → still collides. Order-independent.
- **A5** ✅ `cli.rs:620` is outside the UBam block, so it guards the SE uBAM arm too. Zero
  inputs is unreachable (`required = true`).
- **A6** ✅ premise verified: `trim_galore dup.fastq.gz dup.fastq.gz` exits 0 today with a
  correct 40-read output. So this genuinely is a behaviour change. **I agree with the
  decision to reject** — it matches `--paired`, avoids canonicalisation, and a duplicated
  argument is almost always a mistake. But see Important-3: the message the user will get is
  actively misleading.
- **A7** ✅ `ensure_output_dir` at `:300` precedes all four new insertion points, so a
  rejected run leaves an empty `-o` directory, matching the tested paired contract.

### 2.4 Ordering side effects the plan does not mention (benign, but worth a line)

Inserting the check before the SE loop at `:784` moves the collision error *ahead of*
`sanity_check_any` for inputs 1..n (that call lives inside the loop at `:787`). So
`trim_galore good.fastq.gz truncated.fq.gz` with colliding stems now reports the collision
rather than the truncation. That is the correct precedence — fail before any read — but it
is an observable change in which of two errors surfaces, and someone will eventually notice.

---

## 3. Validation sufficiency

### 3.1 The §9 pairing property does not hold as claimed

§9 opens: "Every rejection test is paired with an acceptance test **on the same dispatch
path**. Without that pairing a helper that rejected unconditionally would pass the whole
suite." Auditing the proposed list:

| Dispatch path | Rejection test | Acceptance test |
|---|---|---|
| SE trim FASTQ | V2, V4, V5 | V3 (two dirs; single input) ✅ |
| `--hardtrim5` FASTQ | V2, V5 | V3 (`-o`, distinct stems) ✅ |
| `--hardtrim3` FASTQ | V2 | **none** ❌ |
| SE trim uBAM | V2 | **none** ❌ |
| `--hardtrim5` uBAM | V2 | **none** ❌ |

The strong version of the property (an unconditionally-rejecting *helper* is caught) does
hold, because two paths do have acceptance partners. The stated version — per dispatch path
— does not, for three of five. A mistake confined to one arm (e.g. Step 8 accidentally
handing the pre-flight a constant, or `planned_hardtrim_outputs` mis-`match`ing) would pass
the whole suite. Two cheap additions close it: one SE-uBAM acceptance case (two distinct
inputs → exit 0, two `.bam`) and one `--hardtrim3` acceptance case.

Mitigating note the plan could make: several plausible slips are *behaviourally inert* for
the collision decision, because every candidate namer is stem-keyed. Passing `"5prime"` in
the hardtrim3 block, or the FASTQ namer on the uBAM arm, yields candidate paths that are
isomorphic to the real ones — colliding exactly when the real outputs collide. That is why
the acceptance gap is Important rather than Critical.

### 3.2 IMPORTANT — the one slip nothing tests: `output_dir` omitted from Step 7

`single_end_output_name` is the only candidate namer whose directory resolution depends on
`output_dir` in a *behaviour-changing* way (`None` → `input.parent()`, i.e. a different
directory per input). If Step 7's candidate builder were written with `None` instead of
`output_dir`:

- inputs `dirA/same.fastq.gz` + `dirB/same.fastq.gz` with `-o /shared` → candidates
  `dirA/same_trimmed.fq.gz` and `dirB/same_trimmed.fq.gz`, distinct → **accepted**, while
  both writes land on `/shared/same_trimmed.fq.gz`. #383 survives, with `-o` set.
- Every proposed test still passes: V2, V4 and V5's SE step all use two inputs in the *same*
  directory, so their candidates collide with or without `output_dir`; V3's acceptance cases
  are unaffected.

§3.3 lists "Same basename, different dirs, **with** `-o` → reject" as an edge case but no
V-item exercises it. **This is the single highest-value test to add** — it is the exact
shape of #383's sibling in a pipeline that stages inputs per-sample and writes to a shared
results directory, and it is the only case where a natural implementation slip in the
plan's own Step 7 goes undetected. (The corresponding hardtrim slip is inert: both hardtrim
namers key on a bare filename either way.)

### 3.3 IMPORTANT — V6 as specified asserts nothing

V6: "unit test in `io.rs` asserting that for a set of paths where
`report_name(a) == report_name(b)`, `single_end_output_name(a) == single_end_output_name(b)`
also holds."

`report_name(a) == report_name(b)` requires equal parent strings and equal file-name
strings — i.e. `a` and `b` are the same path (or trivially separator-equivalent, e.g.
`a/x.fq.gz` vs `a//x.fq.gz`). So any hand-built witness set either contains identical inputs,
making the assertion a tautology, or separator variants, for which both sides collide anyway.
The test will be green on day one and stay green under a naming change that *breaks* A2,
because the hypothesis can never be satisfied by an interesting pair. It is exactly the
class of self-satisfying check the plan's own §9 preamble warns about.

The constructible, meaningful version is the contrapositive with witnesses: assert that for
`sample.fastq.gz` / `sample.fq.gz` / `sample.fastq.bgz` in one directory, the **report names
are all distinct** *while* the **primary names are all equal** — i.e. demonstrate directly
that the primary key is strictly coarser than the report key. That fails loudly if a future
change makes `strip_fastq_extensions` finer-grained than `file_name`, which is the actual
risk A2 carries. Consider adding the same shape for `single_end_bam_output_name`.

### 3.4 V1–V5, V7 — otherwise well-shaped

V1's case-fold assertion genuinely is unreachable from an integration test on APFS, which
justifies Open 2's decision to move the helper into `io.rs`. I agree with that call. V3's
upgrade from existence to per-file read attribution is the right instinct and I confirmed
it discriminates (dirA=ALPHA only / dirB=BETA only). V2's `current_dir(tempdir)` pin is
necessary — I re-confirmed hardtrim writes to CWD. V4 is a worthwhile #382 regression pin.
V5 mirrors the existing guards' shell idiom correctly.

Small gaps: V2 asserts only that the *output* file is absent; the #383 signature was "one
data file, **two reports**, both claiming 100%", so V2 should also assert no
`*_trimming_report.txt`/`.json` exists (the CI step gets this free via `ls -A`, the
integration tests do not). And V7's repro (c) needs the same scratch-CWD discipline V2
gets, or `cargo`-adjacent runs will litter the repo root — the plan's own §11 records this
happening during investigation.

---

## 4. Efficiency

Nothing to object to. O(n) time and memory over a command line's worth of paths; one
`String` key and one `PathBuf` clone per input; no I/O. Running before the ≤1 M-record
adapter scan means colliding invocations fail faster than today, which the plan correctly
notes. `planned_hardtrim_outputs`'s `Vec` is dropped immediately.

Two notes if Critical-1 is adopted: §6's "no syscalls" becomes "one `getcwd`" (or zero, if
the CWD is fetched once and reused), and `HashMap::with_capacity(paths.len())` is a free
micro-win nobody will measure.

---

## 5. Alternatives considered

1. **One hoisted pre-flight above dispatch instead of four insertion points.** Rejected —
   correctly, though the plan does not discuss it. A unified check needs a mode → namer
   `match` that duplicates the dispatch logic it sits above, so it trades four small
   omissions for one large one. The plan's per-branch shape is right for a data-loss fix.
   The residual risk is structural: a fifth dispatch path added later will silently lack the
   guard again, exactly as these four did. Cheapest mitigation is a comment at the helper
   naming all call sites (the `norm_path` doc-comment already uses this pattern), which Step 1
   half-does; extend it to list all guarded paths so an omission is visible from the helper.
2. **Fold duplicate inputs before hashing instead of rejecting (A6).** I prefer the plan's
   rejection, but with a dedicated message (Important-3).
3. **`fs::canonicalize` for the key.** The plan's argument against it is sound (syscall per
   input, symlink and bind-mount semantics). `std::path::absolute` gets most of the benefit
   with none of those costs — see Critical-1.
4. **Refusing to clobber files already on disk** (§10-2) — agreed, separate feature.
5. **Consolidating the four hand-rolled copies** (§10-3) — agreed, separate commit. Splitting
   it out of the data-loss fix is the right instinct.

---

## 6. Action items

### Critical

**C1. The pre-flight does not close #383 for path-aliased directory spellings.** Reproduced
against the current binary: `trim_galore <abs>/sample.fastq.gz sample.fq.gz` (cwd =
`<abs>`) and `trim_galore ./sample.fastq.gz sample.fq.gz` both exit 0 and lose 40 of 80
reads, and the candidate keys the plan hashes are exactly the two textually-different
`Output:` paths the binary prints for the same inode. Because `norm_path` is a pure string
lowercase, the planned helper accepts both invocations. Either (a) key on
`std::path::absolute(p)` — verified to strip `./` and prepend the CWD, lexical, no
filesystem access, message keeps printing the user's original paths so the CI greps stay
valid — and record `..`/symlink aliasing as the remaining limitation; or (b) if that is
deferred, strike the "impossible" framing in §1, add it as an explicit assumption, and say
so in the CHANGELOG, because the issue as filed will still reproduce in one common spelling.
Add a test: two inputs in the same directory, one absolute and one relative, sharing a stem
→ must reject.

### Important

**I1. Add the one test that catches an omitted `output_dir` in Step 7.** SE trim,
`dirA/same.fastq.gz` + `dirB/same.fastq.gz`, **with** `-o <third dir>` → must reject. §3.3
lists this case but no V-item covers it, and it is the only proposed-code slip that
silently under-rejects: every existing V2/V4/V5 SE case uses same-directory inputs whose
candidates collide with or without `output_dir`. See §3.2.

**I2. Narrow A2 — it is false for `--fastqc_args "-o DIR"`.** Verified: two distinct
primaries are accepted while both FastQC artifact sets are written to the overridden
directory under one name, silently clobbering the first. `fastqc.rs:91-96` sets
`config.output_dir` from `--fastqc_args`, and the pre-flight never sees it. Pre-existing and
identical on the `--paired` path, and the loss is a regenerable QC artifact — so scope A2 to
FastQC's *default* directory resolution and record the override as a known residual, or add
the parsed `-o` directory to the candidate list.

**I3. Give the A6 duplicate-input rejection its own message.** With the generic message, the
user who typed `trim_galore dup.fastq.gz dup.fastq.gz` gets "…: `dup_trimmed.fq.gz` and
`dup_trimmed.fq.gz` would be written to the same file. Check that inputs produce distinct
output paths (e.g., different source directories or `--output_dir`)" — two byte-identical
paths and two remedies that both fail to fix it. Mirror `cli.rs:535`'s dedicated wording
("Read 1 and Read 2 appear to be the same file"). Detecting exact-duplicate inputs in
`Cli::validate` is unit-testable, runs before `ensure_output_dir` (so no empty `-o`
directory is left behind), and keeps the pre-flight message meaning what it says.

**I4. Rewrite V6 — as specified it cannot fail.** `report_name(a) == report_name(b)` implies
`a` and `b` are the same path string, so the assertion is a tautology on any constructible
witness set. Invert it: assert that `sample.fastq.gz` / `sample.fq.gz` / `sample.fastq.bgz`
in one directory produce three **distinct** report names and one **shared** primary name,
demonstrating the primary key is strictly coarser. See §3.3.

**I5. §9's stated pairing property holds for only two of five rejection paths.** `--hardtrim3`,
SE trim uBAM and hardtrim uBAM get a rejection test and no acceptance test. Add one
acceptance case each (or weaken the claim and say which arms are rejection-only, and why the
stem-keyed namers make the residual risk small — see §3.1).

### Optional

**O1.** §2.4/§10-3 say "five hand-rolled copies"; there are four (`:482`, `:680`, `:1800`,
`:2412`). `:482` and `:486` are the same block. The count appears twice.

**O2.** V2 should also assert no `*_trimming_report.txt`/`.json` was written. Two reports
next to one data file *is* the #383 signature; existence-of-output-absent alone misses a
regression that writes reports before bailing.

**O3.** V7's repro (c) needs the scratch-CWD discipline V2 gets, or a manual re-run writes
`same.20bp_5prime.fq.gz` into the repo root — §11 records this having happened once already.

**O4.** V2's hardtrim cases use `Command::current_dir(tempdir)`, which breaks the
`fixture()` helper's relative `test_files/…` path used elsewhere in `tests/`. Copy fixtures
into the tempdir or pass absolute paths, and say which.

**O5.** §2.3's `--demux` safety argument ("is handed `&output_path`") is incomplete — being
handed the path says nothing about directory resolution. The actual reason it is safe is
that `demux.rs:142-147` resolves `output_dir.unwrap_or(trimmed_file.parent())`, so its key
is a function of the primary. One sentence, and it makes A2 auditable rather than
plausible.

**O6.** Note in §3.4 that the new SE check moves the collision error ahead of
`sanity_check_any` for inputs 1..n (`main.rs:787`), so a run with both a collision and a
truncated second input now reports the collision. Correct precedence, observable change.

**O7.** A6 is a behaviour change; the CHANGELOG's `### Unreleased` has both
`#### Bug fixes` and `#### Changes` (line 63). Step 12 sends everything to `#### Bug fixes`
— the duplicate-input refusal probably belongs under `#### Changes`.

**O8.** Open 4 (`--hardtrim5` + `--hardtrim3` → 3′ silently ignored) is confirmed
reproducible. Since a reviewer will ask, file it as its own issue before this lands; the
fix is a clap `conflicts_with`.

---

## 7. Summary

The diagnosis is correct and unusually well evidenced — I re-ran every reproduction and
confirmed the hole map is exhaustive by enumerating all ten `cli.input` loops in `main.rs`.
The design (reuse the existing generic helper, move it to `io.rs` for unit-testability,
four minimal per-branch insertions, defer consolidation) is the right shape for a data-loss
fix, and the anchors are accurate enough to implement from. A6's rejection is the right
call.

The one thing I would not ship without addressing is **C1**: `norm_path`'s pure-string key
means a `./` prefix or a mixed absolute/relative argument list walks straight through the
new guard and reproduces #383 exactly, which I demonstrated on the current binary against
the very paths the plan hashes. A one-call lexical absolutisation of the key closes it at
essentially zero cost and strictly improves the four existing guards too. Beyond that, the
validation section is one test short of catching the most plausible slip in its own Step 7
(**I1**), one of its two assumption-tests cannot fail (**I4**), and A2 is a universal claim
with a real counterexample (**I2**).
