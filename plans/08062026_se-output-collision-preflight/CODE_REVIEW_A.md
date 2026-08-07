# Code Review A — output-collision pre-flight (#383)

**Reviewer:** A (independent; a Reviewer B may have run concurrently — no shared state)
**Branch:** `fix/383-output-collision-preflight` off `dev` @ `72624c7`, uncommitted
**Scratch:** `…/scratchpad/creview-a3/` (`p1.sh`–`p6.sh`) and `…/scratchpad/creview-a/`

## No fixes applied

**Nothing was modified** — not source, tests, CI, or the plan. The change is uncommitted in a
shared working tree with a second reviewer possibly active, and a concurrent edit would corrupt
the diff both reviewers read. Every item below is a recommendation, including the one-liners.
The only file this review wrote is this report.

### Provenance note

An earlier Reviewer-A pass wrote a report to this path before it was interrupted. This file
merges it with a second independent pass. Items it raised are kept **only where I reproduced
them myself in this pass** — marked *(reproduced)*. C1 in particular turned out worse than the
earlier pass recorded: it exits **0**, not non-zero. Items neither pass could confirm were
dropped.

---

## Summary

The change does what it sets out to do. Eleven dispatch paths now route through
`io::preflight_output_collisions`; I verified assumption A9 site by site (candidate expression vs.
what the writer actually uses) and found **no mismatch and no fifth hole**. The four hand-rolled
conversions dropped nothing — an 18-case acceptance battery and a 10-case rejection battery both
behave. Defect 1 (`./x` vs `x`) and defect 2 (output aliases an input) are genuinely fixed and
genuinely tested, and there is an unclaimed win: a single-input self-overwrite
(`--basename s s_trimmed.fq`) is now refused where the writer previously truncated its own input.

Four things I would not merge without addressing:

- **C1** — defect 2 is closed for *primary* outputs only. A `--demux` **secondary** output still
  overwrites a file named as an input, destroying its reads **at exit 0**. §11 Out-of-scope 2
  explicitly puts this in scope ("there the user has named the file as an input, so writing it is
  unambiguously wrong"). Reproduced.
- **H1** — the lexical key still lets #383 reproduce verbatim through a `..` component. Disclosed
  as A8, so this is a scope call, not a surprise — but the CHANGELOG states the residual exists
  without stating that it costs the same silent read loss.
- **H2** — `--clock`/`--implicon` name output into the CWD *exactly* like hardtrim, so the generic
  remediation is provably wrong there too, and only hardtrim got the hint built for that purpose.
- **H3** — the hardtrim hint is factually wrong whenever `-o` is given, and contradicts the clause
  it is appended to. The CHANGELOG says it replaces the generic advice; the code appends.

Counts: **1 Critical, 3 High, 4 Medium, 7 Low.**

---

## Verified correct

Recorded so the next reader does not redo it.

**A9 — candidate expression vs. writer, all eleven sites.** Argument for argument, including
`basename` and `gzip`:

| Site | Candidate | Writer |
|---|---|---|
| SE trim FASTQ | `main.rs:796` | `main.rs:1103` (character-identical) |
| SE trim uBAM | `main.rs:1867` | `main.rs:1893` (character-identical) |
| `--paired` FASTQ | `:701` / `:710` / `:725` | `run_paired` `:1293` / `:1313` / `:1305` |
| `--paired` uBAM | `:1818` | `:2028` |
| `--clock` / `--implicon` | closures `:426` / `:440` | same closure reaches the writer arm |
| `--hardtrim5/3` × 2 formats | `planned_hardtrim_outputs` | `specialty.rs:45`/`:88`/`:139`/`:191`, enum-enforced |
| `--clump_only` SE FASTQ / SE uBAM / PE A / PE B | `:511` / `:631` / `:578` / `:542` | `clump_only.rs:278` etc. |

**No fifth hole.** Every writing loop over `cli.input` — `:365`, `:396`, `:514`, `:591`, `:634`,
`:744`, `:800`, `:1828`, `:1870`, `:2433` — has a `preflight_output_collisions` immediately in
front of it. `main.rs:679` (`run_paired_ubam_single_file`) is the one unguarded path; it is safe
because `:1610-1620` builds names from `input.file_stem()` and ignores `--basename`, so
`<stem>_val_{1,2}` can never equal the single BAM input (see L2).

**A2 holds for output-vs-output.** `report_name` (`io.rs:447`), `json_report_name` and demux
(`demux.rs:129-147`) all resolve their directory by the same `output_dir`-else-`input.parent()`
rule and key on `input.file_name()`, so secondary-collide ⇒ primary-collide by construction:
equal dir + equal file name forces equal stem. It does **not** hold for output-vs-input — that is
C1.

**No false rejections.** 18-case acceptance battery (`p2.sh`): `--paired` × {2 pairs,
`--retain_unpaired`, `--passthrough`, `--basename`}, `--clock` × {1 pair, 2 pairs},
`--clump_only` × {SE, PE, SE uBAM}, SE trim × {FASTQ, uBAM}, paired-uBAM out, hardtrim5/3,
`-o .`, and re-trimming a previous output — all exit 0.

**Conversions kept their rejections and gained the alias check.** 10-case rejection battery
(`p3.sh`, `p4.sh`): every duplicate-output case still refuses; new alias refusals fire on SE trim,
`--paired`, `--clump_only` SE and `--clock`; the output directory is empty after each.

**`guarded_inputs()` is the right set.** `cli.input` + `--passthrough`. `passthrough_output_name`
appends `_passthrough`, so it cannot alias its own input, and the candidate at `:725` uses the
same `basename`/`gzip` as the writer at `:1305`.

**`collision_key`'s `Err` fallback cannot weaken the guard.** `std::path::absolute` fails only on
an empty path (clap rejects empty positionals first) or `getcwd` failure, which degrades *every*
key uniformly to the pre-change raw string — so it cannot create a mixed keyspace.

**Ordering (alias before duplicate) is right.** Where both hold, the alias message is the
actionable one and its advice is the correct one.

**D1's predicate is complete.** I looked for a fourth paired-shaped mode that does not set
`paired` and found none: `--clump_only --paired` sets it, hardtrim is genuinely per-file.
Duplicates `validate_paired_input` misses (`--clock a b a c` — neither within-pair nor a duplicate
pair) are still caught downstream by `run_specialty_paired`'s pre-flight.

**Unclaimed win** *(reproduced)*: `trim_galore --basename s s_trimmed.fq` is now refused with
nothing written. Previously `run_single_file` opened the reader (`:1123`) then truncated the same
file with `FastqWriter::create` (`:1124`). Worth a CHANGELOG sentence — it is sharper than the
failure #383 reports.

---

## Critical

### C1 — a `--demux` secondary output still overwrites a named input, destroying reads at exit 0

`main.rs:792-799` builds SE candidates from `single_end_output_name` only. Demux outputs
(`demux.rs:152`, `<primary-stem>_<sample>.fq`) are never compared against the input list.

Reproduced (`p6.sh`, claim 2) — 40 reads destroyed, **exit 0, no warning**:

```bash
printf 'sampleA\tATCG\n' > bc.txt
# m.fastq: 40 reads tagged SRC, each carrying the 3' barcode ATCG
# m_trimmed_sampleA.fq: 40 pre-existing reads tagged PRECIOUS
trim_galore --demux bc.txt m.fastq m_trimmed_sampleA.fq
# rc=0
# PRECIOUS reads in m_trimmed_sampleA.fq: 40 → 0
# SRC      reads in m_trimmed_sampleA.fq:  0 → 40
```

Input 1's demultiplexing pass writes `m_trimmed_sampleA.fq`, which **is input 2**; input 2 is then
trimmed from the file that was just overwritten, so `m_trimmed_sampleA_trimmed.fq` also holds the
wrong sample. That is #383's exact signature — silent loss plus a wrong-sample file — on the very
dispatch path this change was written to protect.

A2 does not cover this. A2 reasons about output-vs-**output** distinctness (defect 1's domain);
defect 2 is output-vs-**input**, and primary distinctness says nothing about whether a *secondary*
equals an input (`m_trimmed.fq` ≠ `m_trimmed_sampleA.fq`). Nor is it §11 Out-of-scope 2, whose own
wording carves it in: the user named the file as an input. And unlike A2's `--fastqc_args -o`
exception, the cost here is reads, not a regenerable QC artifact.

Reachable by a plain re-run glob: demux outputs carry the input extension, and `m.fastq` sorts
before `m_trimmed_sampleA.fq` under both C and en_US collation (`.` = 0x2E < `_` = 0x5F), so
`trim_galore --demux bc.txt *.f*q` over a directory holding a previous demux run lands in the
destructive order.

Cheapest containment: `--demux` is single-end-only and the barcode names are known before any
write, so add the prospective demux names to the candidate list at `main.rs:792` when
`cli.demux.is_some()`. If that is deferred, the CHANGELOG must stop saying defect 2 "is now
refused on every path", and A2 must state that output-vs-input is checked for primaries only.

---

## High

### H1 — `..` defeats the key; #383 reproduces verbatim *(reproduced)*

```bash
# cwd = work/, data/ is a sibling
trim_galore ../data/s.fastq /abs/…/data/s.fq
# Output: ../data/s_trimmed.fq   and   /abs/…/data/s_trimmed.fq
# rc=0; ALPHA reads surviving: 0, BETA: 40; both misleading reports present
```

This is A8 / §11 Out-of-scope 3 and the CHANGELOG does disclose it. My objection is to the scope
decision, not the honesty: the plan closes `./x`-vs-`x` and abs-vs-rel *because* mixed argument
lists are "the common shape" (§2.2(e)), and `../data/x` mixed with an absolute path is that same
shape — it is what a workflow manager's relative-to-workdir input looks like.

Two fixes, both consistent with A1's "loud error beats silent loss" trade:

- *Minimal, still lexical:* after `absolute()`, fold `Component::ParentDir` against the accumulated
  component stack. Can only increase rejections; the only false positive is a symlinked
  intermediate directory — the same class `norm_path`'s case-fold already accepts.
- *Thorough:* `fs::canonicalize` the path's **parent** and join the file name. Closes symlink
  aliasing too, one syscall per path, and retires A8 rather than documenting it.

If A8 stays, the CHANGELOG bullet should say what the residual *costs* (silent read loss identical
to #383), not merely that it exists.

### H2 — `--clock`/`--implicon` need the hint as much as hardtrim does *(reproduced)*

`clock_output_name` (`specialty.rs:488-491`) and `implicon_output_name` (`:521-524`) end
`None => PathBuf::from(filename)` / `Some(dir) => dir.join(filename)` with a **stem-only**
filename — identical to both hardtrim namers, as §2.3 itself records. So the collision class is
identically wide and `--output_dir` identically fails to rescue it, yet `run_specialty_paired`
passes `hint: None` (`main.rs:2431`):

```bash
trim_galore --clock -o out dirA/same_R1.fastq.gz dirA/same_R2.fastq.gz \
                           dirB/same_R1.fastq.gz dirB/same_R2.fastq.gz
# Error: … out/same_R1.clock_UMI.R1.fq.gz and out/same_R1.clock_UMI.R1.fq.gz would be written
# to the same file. Check that inputs produce distinct output paths (e.g., different source
# directories or `--output_dir`).      ← the user already did both of these
```

Multi-pair clock is a supported, CI-tested shape, so `--clock */*_R{1,2}.fq.gz` over a per-sample
layout now hard-fails with advice that provably does not work — the outcome §3.4 exists to
prevent. Fix: generalise the hint to a CWD-output hint (the sentence is already mode-neutral apart
from the words "Hard-trimmed") and thread it through `run_specialty_paired`, which must keep
`None` for `--clump_only --paired` (an `input.parent()` namer).

### H3 — the hardtrim hint is wrong with `-o`, and contradicts the clause it is appended to

`main.rs:55-59`, appended at `io.rs:81`. Rendered message with `-o` — the shape the new test *and*
the new CI step both use:

> … `out/same.20bp_5prime.fq` and `out/same.20bp_5prime.fq` would be written to the same file.
> Check that inputs produce distinct output paths (e.g., **different source directories or
> `--output_dir`**). Hard-trimmed output is written to the **current working directory**, so
> inputs sharing a basename collide **whatever `--output_dir` is set to** — run one invocation per
> input, or give the inputs distinct basenames.

Two defects in one sentence. The paths it just printed are inside `out/`, so "written to the
current working directory" is visibly false — `hardtrim_output_name` (`specialty.rs:449-452`)
honours `-o`; what it ignores is the **input's parent**. And the generic clause it is appended to
recommends exactly the two remedies §3.4 established do not work, so the user reads "use different
source directories or `--output_dir`" immediately followed by "that will not help".

`CHANGELOG.md` says these modes "add their own remediation to the error **rather than** the generic
advice". The code concatenates. Smallest fix that matches the CHANGELOG: build the advice sentence
as `hint.unwrap_or(GENERIC_ADVICE)` instead of appending. Repro: `p1.sh`, or
`trim_galore --hardtrim5 20 -o out dirA/same.fastq dirB/same.fastq`.

Neither the test nor the CI step can catch this: both assert only that `one invocation per input`
is present, never that the rest of the sentence is true.

---

## Medium

### M1 — the A2 table test does not test A2's direction, and drops most of V5

`io.rs`, `primary_output_key_is_coarser_than_secondary_keys`. A2 needs *primary distinct ⇒
secondary distinct*. The test instead picks three variants whose **primaries are equal** and
asserts their **secondaries are distinct** — a witness that the *converse* fails, which is
reassuring but is not evidence for A2. All three variants also live in one directory and have
pairwise-distinct `file_name()`, which every secondary namer embeds verbatim, so the `assert_ne!`
loop is close to a tautology for any implementation of that shape; the real failure mode it should
detect — a secondary namer that ignores the directory — is invisible to it.

It also departs from V5 with no §13 entry: V5 asked for pairs whose **primaries differ**, covering
the demux base name and the default-configuration `<stem>_fastqc.zip`, across `--basename` /
`--dont_gzip` / `-o` on and off. The implemented table covers none of those. The property does
hold today (see "Verified correct"), so this is rigour, not a live bug — but V5 existed precisely
because v1's equivalent could only ever pass. Keep the first half (all three primaries equal —
that pins `.bgz` stem stripping); replace the `assert_ne!` loop with same-basename inputs in
**different** directories, asserting secondaries differ without `-o` and collide with it.

### M2 — the duplicate message names neither input, and the alias message names only one

When the two colliding paths are textually identical — the common case on every mode — the message
repeats itself and identifies nothing:

```
out/same.20bp_5prime.fq and out/same.20bp_5prime.fq would be written to the same file
out/same_R1_val_1.fq and out/same_R1_val_1.fq would be written to the same file
```

The user is told to make the inputs produce distinct paths, but not which two inputs to look at.
The alias branch (`io.rs:66-73`) has the mirror problem: it prints the aliased path but not the
input that would clobber it, where the duplicate branch names both participants. On a 40-file
command line neither is actionable. Pass `(candidate, source_input)` pairs so both messages can
name inputs; or minimally special-case `existing == p` to "two inputs produce the same output
path: X". Pre-existing wording inherited from all four hand-rolled copies — but it is now the one
shared message, so this is the moment.

### M3 — `Cli::validate`'s duplicate-input check is byte-equality, so two spellings get the wrong message

`cli.rs:631` compares `other == path`:

```bash
trim_galore dup.fq dup.fq     # → "Input file dup.fq was given more than once …"   (intended)
trim_galore ./dup.fq dup.fq   # → "Output path collision (case-insensitive, for APFS/NTFS …"
```

§4.4's whole rationale is that a duplicated input should get a precise message "rather than the
case-insensitive output-collision pre-flight's APFS/NTFS message", and A6 explicitly worries about
"an inconsistency worse than either clean choice". Both shapes are rejected, so this is diagnostic
quality only. Keying that loop with the same normalisation (`collision_key` as `pub(crate)`) makes
the two consistent.

### M4 — `--hardtrim3`'s hint is untested, one argument from the gap the plan designed around

`tests/integration_output_collision.rs:360-379` asserts the prefix and `DUP_MSG` but not
`one invocation per input`; only the `--hardtrim5` test (`:354`) and the `--hardtrim5` CI step check
the hint. §4.2 introduced `HardtrimEnd` precisely because a copy-paste slip in the `--hardtrim3`
block "would produce a pre-flight that hashes paths the run never writes, while every rejection
test still passed" — the enum closes that for the end discriminator and nothing closes it for the
adjacent `hint` argument. I confirmed `main.rs:391` does pass it, so this is coverage, not a
defect. One line closes it.

---

## Low

**L1 — "(arguments N and M)" indexes positionals, not argv.** `cli.rs:632-639`;
`trim_galore -q 20 a.fq a.fq` reports "arguments 1 and 2" for argv items 4 and 5. Suggest "input
files 1 and 2".

**L2 — the helper's doc-comment overstates and miscounts.** `io.rs:51-54` says "all four output
formats"; there are two (`Fastq`, `UBam`) — it presumably means the four `--clump_only` arms. And
"every `main.rs` dispatch path that writes more than one file" is untrue: `run_paired_ubam_single_file`
writes 2 (4 with `--retain_unpaired`) and does not call it. It is safe, but §4.1's stated purpose
is to make a fifth omission visible *from the helper*, which only works if the list is exact. Name
that path as deliberately exempt and why.

**L3 — `HARDTRIM_HINT` carries its own leading space** (`main.rs:56`) and the helper appends it raw
via `hint.unwrap_or("")`. Correct today (verified in the rendered message), but a future caller
that forgets the space glues two sentences together. Prefer the helper owning the separator.

**L4 — stale line references in code this change touched.** `main.rs:2422` says "the `--paired`
pre-flight at lines 82–125" (now `:699-734`); `ci.yml` still points at `src/main.rs:101–111` for the
`--retain_unpaired` branch (now `:709-719`). Both predate the change.

**L5 — `input_keys: HashMap<String, &PathBuf>` collapses two inputs sharing a key** to whichever
was listed last (`.collect()` keeps the last value), so the alias message can name a spelling the
user did not type. Detection is unaffected — the surviving entry has the same key. `HashSet` for
membership plus `entry().or_insert()` if the `&PathBuf` is wanted makes it deterministic.

**L6 — D5's YAML `#` truncation is worth fixing while the steps are new.** The three new step names
lose everything from ` #383` onward, so the Actions UI shows "… pre-flight (issue". Consistency
with the existing `(issue #216)` step is defensible, but quoting three new names is a one-line
change and the consistency argument then applies in the other direction.

**L7 — one fixture passes for a reason its comment does not give.**
`se_trim_rejects_gz_and_bgz_sharing_a_stem` (`:196`) writes plain text into `wide.fastq.gz` /
`.bgz`; the comment says "Content is irrelevant: the pre-flight runs before any read", true only
because `detect_input_format` peeks the first byte and classifies `@…` as plain FASTQ while `gzip`
is derived from the *filename*. It does exercise the intended path today, and
`assert_rejected_cleanly` requires the collision prefix so it cannot pass for an unrelated
failure — but a clause naming the filename-vs-content split would stop the next reader trusting
the wrong reason.

---

## Vacuous-pass audit (§13 D3's failure class)

| Test | Can it fail? | Note |
|---|---|---|
| `preflight_rejects_dot_slash_alias`, `…absolute_versus_relative_alias` | yes | fail without `absolute()` |
| `preflight_rejects_output_that_aliases_an_input` | yes | also asserts *absence* of the duplicate wording |
| `preflight_accepts_dotdot_alias_known_limitation` | yes | inverts if the key gains `canonicalize` |
| `primary_output_key_is_coarser_than_secondary_keys` | yes, but not on A2 | **M1** |
| `se_trim_rejects_dot_slash_alias`, `…absolute…` | yes | D3's fix (drop `-o`) is right; `assert_dir_holds_only` pins it |
| `se_trim_rejects_same_basename_across_dirs_with_output_dir` | yes | the `output_dir`-vs-`None` slip detector; its acceptance twin at `:241` is what makes it non-trivial |
| `hardtrim3_rejects_same_basename_across_dirs` | partly | hint unasserted — **M4** |
| hardtrim5/3 acceptance tests | yes | filename assertions catch a wrong `HardtrimEnd` |
| `paired_rejects_output_that_aliases_an_input` | yes | all four planned outputs are distinct, so only the input-key set can reject |
| `duplicate_input_gets_a_precise_message` | yes | asserts absence of `APFS/NTFS` |
| CI steps 1 and 3 | yes | `rc != 0` + message + residue; step 3 also md5-pins the input |
| CI step 2 | presence yes, accuracy no | **H3** |

`tempdir()`'s `canonicalize` (D4) is load-bearing and correctly commented: without it the
absolute-vs-relative case compares `/tmp/…` against the child's `/private/tmp/…` `getcwd` and
passes vacuously. I confirmed the symlink is real on this machine, so D4 is a genuine catch.

Acceptance coverage note: the new file has no acceptance case for `--clump_only` or
`run_specialty_paired`. Both are covered elsewhere (`tests/integration_clump_only*.rs`, and the
multi-pair `--clock`/`--implicon` CI steps), so this is a completeness note — worth a sentence in
the module doc-comment so a reader does not conclude those two conversions are untested.

## Efficiency

Nothing to report. One `getcwd` per path via `std::path::absolute`, one `String` per key, one
`PathBuf` clone per candidate, `with_capacity` on the planned map — unmeasurable at command-line
n. `guarded_inputs()` clones the input vector per call, but exactly one site executes per run. The
guard now runs before adapter auto-detection's ≤1 M-record scan, so colliding invocations fail in
microseconds instead of after a full detection pass — a real improvement on the four new paths.
