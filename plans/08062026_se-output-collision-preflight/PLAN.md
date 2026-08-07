# Plan — Output-collision pre-flight: close the four uncovered paths and the two key defects (#383)

**Issue:** [#383](https://github.com/FelixKrueger/TrimGalore/issues/383)
**Branch:** `dev` @ `72624c7`
**Revision:** v2 (2026-08-06) — see §12

---

## 1. Goal

Make every prospective output path of a run checked, against every other prospective
output path **and** against every input path, before any file is opened. Fix #383 and the
three sibling holes on the same dispatch pattern.

Two defects in v1's design, each found independently by one reviewer and reproduced by the
orchestrator, are in scope:

- **The key was a raw string.** `./sample.fastq.gz` and `sample.fastq.gz` hash differently
  while naming one inode, so #383 reproduced verbatim through v1's own guard.
- **The check was output-vs-output only.** A planned output can be an *input* of the same
  invocation; the run then destroys that input's reads and emits a second output holding
  the wrong sample's data. `--paired` has this hole today, with the pre-flight in place.

Scope note: closing the second defect requires the shared helper to see the input list, and
`--paired` needs it as much as the new paths do. The four existing hand-rolled copies
therefore **convert to the shared helper in this commit** — what v1 deferred as a cosmetic
"commit 2" is now load-bearing, and it also corrects two error messages that advise a flag
that does not parse (§2.5).

---

## 2. Context

### 2.1 The verified hole map

Measured against `target/release/trim_galore` at `72624c7`. Every row was executed; read
counts come from grepping read-ID prefixes in the surviving output.

| Dispatch path | Pre-flight today | Measured on a colliding invocation |
|---|---|---|
| `--paired` trim FASTQ | ✅ `main.rs:672-727` | exit 1, message, nothing written |
| `--paired` trim uBAM (2-file) | ✅ `main.rs:1800-1820` | — |
| `--paired` trim uBAM (1-file interleaved) | n/a — one output | — |
| `--clock` / `--implicon` | ✅ `run_specialty_paired`, `main.rs:2412-2429` | exit 1, message |
| `--clump_only` SE FASTQ | ✅ `main.rs:481-498` | exit 1, message |
| `--clump_only` SE uBAM | ✅ `main.rs:609-614` | exit 1, message |
| `--clump_only` PE (both formats) | ✅ 2 sites | exit 1, message |
| **SE trim FASTQ** | ❌ `main.rs:781-792` | **exit 0, 40 of 80 reads lost** |
| **SE trim uBAM** | ❌ `main.rs:1857-1865` | **exit 0**, one `.bam`, two reports |
| **`--hardtrim5/3` FASTQ** | ❌ `main.rs:344-393` | **exit 0, 40 of 80 reads lost** |
| **`--hardtrim5/3` uBAM** | ❌ `main.rs:344-393` | **exit 0**, one `.bam` written twice |

Both reviewers independently enumerated all ten `cli.input` loops in `main()` and confirmed
exactly these four lack a preceding pre-flight. **No fifth hole exists.**

### 2.2 Reproductions (all re-run at v2)

```bash
# (a) #383 as filed — same stem, different FASTQ extension
trim_galore sample.fastq.gz sample.fq.gz
# → exit 0; one sample_trimmed.fq.gz; ALPHA 0 / BETA 40; TWO reports, each "40 (100.0%)"

# (b) the .bgz widening from #382
trim_galore wide.fastq.gz wide.fastq.bgz          # → exit 0; GZALPHA 0 / BGZBETA 40

# (c) hardtrim, same basename in DIFFERENT directories (outputs land in CWD)
trim_galore --hardtrim5 20 dirA/same.fastq.gz dirB/same.fastq.gz
# → exit 0; one ./same.20bp_5prime.fq.gz; DIRA 0 / DIRB 40

# (d) hardtrim uBAM arm
trim_galore --hardtrim5 20 --output-format ubam s_R1.fastq.gz s_R1.fq.gz
# → exit 0; "Writing hard-trimmed" twice; one .bam

# (e) DEFECT 1 — v1's guard accepts these; both lose 40 of 80 reads at exit 0
cd d && trim_galore /abs/d/sample.fastq.gz sample.fq.gz
#   Output: /abs/d/sample_trimmed.fq.gz  +  Output: sample_trimmed.fq.gz   (one inode)
cd d && trim_galore ./sample.fastq.gz sample.fq.gz
#   Output: ./sample_trimmed.fq.gz       +  Output: sample_trimmed.fq.gz   (one inode)

# (f) DEFECT 2 — a planned output IS an input; --paired is affected too
trim_galore s.fastq.gz s_trimmed.fq.gz
# → exit 0; s_trimmed.fq.gz overwritten with input 1's reads (prior 40 destroyed);
#   s_trimmed_trimmed.fq.gz ALSO holds input 1's reads — a file of the wrong sample
trim_galore --paired a_R1.fq.gz a_R2.fq.gz a_R1_val_1.fq.gz a_R2_val_2.fq.gz
# → exit 0; a_R1_val_1.fq.gz's 40 prior reads destroyed; a_R1_val_1_val_1.fq.gz created
```

Reachability, stated precisely so the CHANGELOG can be too:

- (e) needs an argument list that mixes spellings — `find`/`xargs` output, staged pipeline
  inputs, hand-assembled lists. A *uniform* glob (`*.fastq.gz` or `./*.fastq.gz`) is safe,
  because the prefix is then consistent. Mixed lists are the common shape.
- (f) needs a previous run's output in the input list, with the source file ordered
  **first**. A bare `*.gz` glob is safe under the default shell collation on this machine —
  both zsh and bash expand to `s_trimmed.fq.gz s.fastq.gz`, the harmless order. The
  destructive order arises from C-locale sorting (`LC_ALL=C`), an `ls | sort` pipeline, or an
  explicit list. Reviewer B's report claims a bare glob suffices; that is **not** reproducible
  here, and the plan does not rely on it.

### 2.3 Two naming asymmetries the plan depends on

**Output directory.** With no `--output_dir`, `io.rs` trim namers end
`None => input.parent()` (`io.rs:101`, `:121`, `:189`), while all four `specialty.rs` namers
end `None => PathBuf::from(filename)` (`specialty.rs:435`, `:456`, `:471`, `:504`) — a bare
filename, hence **CWD**. Measured consequence: SE trim on two same-basename inputs from
different directories does not collide (2 outputs, correct attribution), but `--hardtrim5`
on the same two does. Hardtrim's collision class is strictly wider and includes the
ordinary `--hardtrim5 20 */*.fastq.gz`. Neither naming rule changes here (§11 Out of scope 1).

**Secondary outputs.** See A2 in §9 for the narrowed claim and its one exception.

### 2.4 What already exists

`preflight_collision_bam` (`main.rs:64-80`) is path-generic — `&[PathBuf]`, hashes
`naming::norm_path`, bails on the first duplicate. Only its name and doc-comment say BAM.
Beside it sit **four** hand-rolled copies of the same HashMap-and-bail: `main.rs:482-498`
(`--clump_only` SE FASTQ), `:680-727` (paired FASTQ), `:1800-1820` (paired uBAM),
`:2412-2429` (`run_specialty_paired`). Three call sites use the helper: `:533`, `:566`, `:614`.

### 2.5 `--output-dir` is not a valid flag

`cli.rs:152` declares `short = 'o', long = "output_dir"` only. Verified:
`trim_galore --output-dir /tmp …` → `error: unexpected argument '--output-dir' found`.
Two of the four hand-rolled copies (`:721` paired FASTQ, `:1815` paired uBAM) advise
`--output-dir` in their error text. Converting them to the shared helper replaces that with
`--output_dir`, so this commit fixes a user-visible wrong instruction as a side effect.

---

## 3. Behavior

### 3.1 The collision key

```
key(p) = norm_path(std::path::absolute(p).unwrap_or_else(|_| p.to_path_buf()))
```

`std::path::absolute` is lexical — no filesystem access, no symlink resolution — and
stabilised in Rust 1.79, below the crate's 1.88 floor. Verified behaviour:

| input | absolute() |
|---|---|
| `sample_trimmed.fq.gz` | `<cwd>/sample_trimmed.fq.gz` |
| `./sample_trimmed.fq.gz` | `<cwd>/sample_trimmed.fq.gz` — `.` removed |
| `a//sample_trimmed.fq.gz` | `<cwd>/a/sample_trimmed.fq.gz` — separators collapsed |
| `a/../sample_trimmed.fq.gz` | `<cwd>/a/../sample_trimmed.fq.gz` — `..` **kept** |
| `/tmp/x/sample_trimmed.fq.gz` | unchanged |

Absolutising can only ever *increase* rejections: two paths that absolutise equal are the
same file. `..` and symlinks still alias — a stated residual (A8), not a silent one. The
`Err` fallback (empty path, or `getcwd` failure) degrades to v1's raw-string key rather
than skipping the check.

**Messages print the user's original paths**, never the absolutised key, so the existing CI
greps and every "nothing written" assertion stay valid.

### 3.2 Two distinct rejections

1. **Duplicate output.** Two prospective outputs share a key. Existing wording, canonical
   `--output_dir`, plus a mode-specific hint where the generic advice is wrong (§3.4).
2. **Output aliases an input.** A prospective output shares a key with a file this run
   reads. Its own message — "would be written to the same file" is untrue and its advice
   inapplicable. Names the input, says it looks like a previous run's output, and offers
   `--output_dir` (which does help in this case, unlike case 1 on hardtrim).

The input set is `cli.input` plus `--passthrough` when set. The `--demux` barcode file is
excluded by argument: every output name carries a `_trimmed` / `_val_` / `_unpaired_` /
`.Nbp_` / `_UMI_` segment plus a FASTQ or BAM extension, so it cannot equal a
user-supplied barcode text file except by deliberate contrivance.

### 3.3 Pre-flight contract

1. Runs before any reader is opened **on the guarded path** — hence before adapter
   auto-detection, the trimming banner, and any output file. Not before all I/O:
   `sanity_check_any(&cli.input[0])` (`main.rs:153`) and `detect_input_format` over every
   input (`:159`) precede dispatch, so **format and sanity errors take precedence over
   collision errors**. That ordering is correct and now stated.
2. Builds the full prospective-output list for **every** input, not per file.
3. On the first violation: `anyhow::bail!`, exit 1, nothing written.
4. Single input with no alias, or all-distinct keys: no-op.

### 3.4 Mode-specific remediation on hardtrim

The generic advice — "different source directories or `--output_dir`" — is wrong on
`--hardtrim5/3`, because the namer ignores the input's parent and `-o` joins a stem-only
filename. Verified: `--hardtrim5 20 -o out dirA/same.fastq.gz dirB/same.fastq.gz` still
collides inside `out`. Both suggestions lead back to the same error.

This matters because the mode's multi-file use is documented:
`docs/src/content/docs/modes/hardtrim.md` and the carried v0.6.x changelog text
(`docs/…/reference/changelog.md:955`, `:997`) both say it processes "one or more files",
and neither mentions CWD output. So `--hardtrim5 N */*.fastq.gz` is a documented pattern
that becomes a hard failure — correctly, since it loses data today, but it must fail with
usable advice. The two hardtrim call sites pass a hint naming the real remedies: one
invocation per input, or distinct input basenames.

### 3.5 Edge cases

| Case | After change | Why |
|---|---|---|
| One input, no alias | no-op | nothing to compare |
| `sample.fastq.gz` + `sample.fq.gz` | **reject** | #383 as filed |
| `sample.fastq.gz` + `sample.fastq.bgz` | **reject** | #382 widening |
| `./x.fastq.gz` + `x.fq.gz` | **reject** (new in v2) | absolutised key |
| `/abs/x.fastq.gz` + `x.fq.gz`, cwd `/abs` | **reject** (new in v2) | absolutised key |
| `a/../x.fastq.gz` + `x.fq.gz` | **accept** — known residual | `absolute` keeps `..` (A8) |
| Output path equals an input path | **reject** (new in v2) | defect 2; own message |
| Same file listed twice | **reject** in `Cli::validate` | A6; own message, before `ensure_output_dir` |
| Case-only path variants | **reject** | `norm_path` case-fold, #216 |
| Same basename, different dirs, no `-o` — SE trim | **accept** | outputs genuinely distinct |
| Same basename, different dirs, no `-o` — hardtrim | **reject** | outputs identical in CWD |
| Same basename, different dirs, **with** `-o` | **reject** | both land in `-o`; V2.3 pins this |
| `--basename` with >1 SE input | unchanged, already rejected | `cli.rs:620` |
| Zero inputs | unreachable | clap `required` |
| `--hardtrim5 N --hardtrim3 M` | unchanged | `hardtrim5` returns at `:367`; §11 Open 3 |

---

## 4. Signatures

### 4.1 `src/io.rs` — the shared helper

```rust
/// Case-folded, lexically-absolutised collision key (issues #216, #383).
fn collision_key(p: &Path) -> String {
    norm_path(&std::path::absolute(p).unwrap_or_else(|_| p.to_path_buf()))
}

/// Reject duplicate outputs, or an output that would overwrite an input, before
/// any file is opened. Called by every dispatch path in `main.rs` that writes
/// more than one file: SE trim (FASTQ + uBAM), `--paired` trim (FASTQ + uBAM),
/// `--hardtrim5/3`, `--clock`/`--implicon`, and all four `--clump_only` arms.
pub fn preflight_output_collisions(
    planned: &[PathBuf],
    inputs: &[PathBuf],
    hint: Option<&str>,
) -> Result<()>
```

Behaviour: build `HashSet<String>` of input keys; then walk `planned`, checking each against
the input set first and against a `HashMap<String, PathBuf>` of prior planned paths second.
Both messages print original paths. `hint`, when `Some`, is appended to the duplicate-output
message.

The doc-comment lists all call sites deliberately: this plan fixes the fourth instance of
one omission, and an enumeration at the helper makes a fifth omission visible from here.
(Alternative considered: an `OutputPlan { planned, inputs, hint }` struct. Rejected as
ceremony for eleven call sites, nine of which pass `None`.)

### 4.2 `src/specialty.rs` — visibility and a typed discriminator

```rust
pub enum HardtrimEnd { Five, Three }   // `.as_str()` → "5prime" / "3prime"

pub fn hardtrim_output_name(input: &Path, keep: usize, end: HardtrimEnd, output_dir: Option<&Path>, gzip: bool) -> PathBuf
pub fn hardtrim_bam_output_name(input: &Path, keep: usize, end: HardtrimEnd, output_dir: Option<&Path>) -> PathBuf
```

`fn` → `pub fn` follows the precedent of `clock_output_name` / `implicon_output_name`, which
are already `pub` so `main()` can build collision candidates. The `&str` → enum change exists
because after this plan the discriminator is written in two places — the candidate builder and
the writer — and a copy-paste slip passing `"5prime"` in the `--hardtrim3` block would produce
a pre-flight that hashes paths the run never writes, while every rejection test still passed.
The enum makes that slip a compile error instead of a test gap.

### 4.3 `src/main.rs`

```rust
fn planned_hardtrim_outputs(
    cli: &Cli,
    keep: usize,
    end: specialty::HardtrimEnd,
    output_dir: Option<&Path>,
    gzip: bool,
) -> Vec<std::path::PathBuf>
```

`std::path::PathBuf` spelled in full: `main.rs:5` imports `Path` only, as the existing helper
does. Same for every `Vec<PathBuf>` in §5.

### 4.4 `src/cli.rs` — A6 moves here

A single-end analogue of the duplicate-pair check, beside `cli.rs:534-558`, with its own
message in the style of `:535` ("Read 1 and Read 2 appear to be the same file"). Placing it
in `Cli::validate` gives a precise message, runs before `ensure_output_dir` so no empty `-o`
directory is left behind, covers every mode in one place, and is unit-testable. `cli.rs:504-516`'s
own doc-comment already states that duplicate inputs should get a precise error "rather than
the case-insensitive output-collision pre-flight's APFS/NTFS message" — v1 did the thing that
comment argues against.

---

## 5. Implementation outline

Anchors are against `72624c7` and shift as steps apply; each carries a searchable token.

**Step 1 — `src/io.rs`: add `collision_key` and `preflight_output_collisions`** (§4.1), after
`norm_path` (`io.rs:44`). Extend `norm_path`'s doc-comment to name the new caller.

**Step 2 — `src/main.rs`: delete `preflight_collision_bam`** (`:60-80`) and re-point its three
calls (`:533`, `:566`, `:614`) to `naming::preflight_output_collisions(&planned, &cli.input, None)`.
Incidental fix that falls out: the doc-comment at `:55-59` describes `resolve_clump_layout`'s
clumpify-budget fallback but is attached to the deleted function; removing it rejoins the comment
to `resolve_clump_layout` (`:82`).

**Step 3 — convert the four hand-rolled copies** to the shared helper, preserving each site's
existing semantics:
- `:482-498` `--clump_only` SE FASTQ — one candidate per input.
- `:680-727` paired FASTQ — accumulate candidates across **all** chunks into one `Vec` before
  the single call, so cross-pair collisions are still caught; keep the `--retain_unpaired` and
  `--passthrough` extra candidates.
- `:1800-1820` paired uBAM — one candidate per pair.
- `:2412-2429` `run_specialty_paired` — collect both names per pair across all pairs, then call once.
Every site passes `&cli.input` (plus `--passthrough` where applicable) as `inputs`, so all four
gain defect-2 protection. Note the two `--output-dir` message sites disappear here (§2.5).

**Step 4 — `src/specialty.rs`: `HardtrimEnd` enum + `pub`** on both namers (§4.2); update the
four internal call sites in `hardtrim5`/`hardtrim3`/`hardtrim5_to_bam`/`hardtrim3_to_bam` and
the `bgz_stem_reaches_specialty_output_names` test at `:822`.

**Step 5 — `src/main.rs`: add `planned_hardtrim_outputs`** (§4.3), matching on `cli.output_format`.

**Step 6 — guard `--hardtrim5`** at `:344`, first statement in the block:

```rust
// #383 — collide across all inputs before the first write.
naming::preflight_output_collisions(
    &planned_hardtrim_outputs(&cli, n, specialty::HardtrimEnd::Five, output_dir, gzip),
    &cli.input,
    Some(HARDTRIM_HINT),
)?;
```

with `HARDTRIM_HINT` a `const &str` naming the real remedies (§3.4).

**Step 7 — guard `--hardtrim3`** at `:369`, identically with `HardtrimEnd::Three`.

**Step 8 — guard SE trim FASTQ** at `:781`, in the `} else {` arm, before the loop:

```rust
// #383 — SE was the only trim path without the #216 pre-flight.
let planned: Vec<std::path::PathBuf> = cli
    .input
    .iter()
    .map(|input| naming::single_end_output_name(input, output_dir, cli.basename.as_deref(), gzip))
    .collect();
naming::preflight_output_collisions(&planned, &cli.input, None)?;
```

**Step 9 — guard SE trim uBAM** at `:1857`, same shape with `single_end_bam_output_name`.

**Step 10 — `src/cli.rs`: reject duplicate SE inputs** (§4.4) with its own message.

**Step 11 — tests.** Unit in `io.rs` (V1) and `cli.rs` (V4); new
`tests/integration_output_collision.rs` (V2, V3) following
`tests/integration_paired_format_guard.rs` conventions — `env!("CARGO_BIN_EXE_trim_galore")`,
a `tempdir(tag)` keyed on `std::process::id()`, `(bool, String)` from `Command::output()`.
**Any case that sets `Command::current_dir` must pass absolutised fixture paths** — `test_files/…`
is relative to the crate root, so otherwise the run fails with "Input file not found" and the
`!ok` half of the assertion passes for the wrong reason.

**Step 12 — CI**: three validation steps (V6).

**Step 13 — CHANGELOG**: `#### Bug fixes` for #383 and the two defects; `#### Changes` (the
section already exists at `:63`) for the A6 duplicate-input refusal and the
`--output-dir` → `--output_dir` message correction.

**Comment discipline:** one line per new comment. Detail belongs in the commit message; the
surrounding `main.rs` comments are far longer and that divergence is deliberate.

---

## 6. Efficiency

O(n) over a command line's worth of paths: one `HashSet<String>` of input keys, one
`HashMap<String, PathBuf>` of planned keys, one `String` per key, one `PathBuf` clone per
planned path. One `getcwd` per `absolute` call — n + m calls, or one if the CWD is fetched
once and reused; either is unmeasurable at these n.

The check runs before adapter auto-detection's ≤1 M-record scan, so colliding invocations
fail in microseconds rather than after a full detection pass. `HashMap::with_capacity` is a
free micro-win nobody will measure.

---

## 7. Integration

**Reads:** `cli.input`, `cli.passthrough`, `cli.output_format`, `cli.basename`, `output_dir`,
`gzip`. Nothing new from disk except `getcwd`.

**Writes:** nothing. The change only prevents writes.

**Order:** after `Cli::validate()` (which now also rejects duplicate SE inputs), after the
`--phred64`/BAM-format guards, after `ensure_output_dir`, before every reader open on the
guarded path.

**Downstream:**

- **Validation matrix — unaffected.** Both reviewers independently read every `trim_galore`
  invocation in `ci.yml`: no step passes more than one single-end input, and multi-pair steps
  use distinct `A_`/`B_`/`C_` stems. The two existing collision guards still reject
  (distinct dirs into a shared `-o`), and their `test -z "$(ls -A …)"` assertions (`:615`,
  `:640`) still hold. No byte-identity path changes.
- **`docs/`** — `modes/hardtrim.md` advertises "one or more files" without noting CWD output.
  A doc sentence is warranted, but it is a docs change, not part of this fix; recorded as
  §11 Open 4 so it is not lost.
- **Existing tests** — no test asserts that a colliding SE or hardtrim invocation succeeds.
  The four converted sites keep their behaviour, so `--paired`/`--clock`/`--implicon`/
  `--clump_only` tests should pass unchanged; V7's `cargo test` confirms rather than assumes.
- **nf-core / MultiQC** — a pipeline that today produces one trimmed file and two reports from
  a colliding SE invocation now exits non-zero. Newly loud, not newly broken.

---

## 8. Signatures compile

Checked, because v1 asserted this and was wrong: `main.rs:5` is `use std::path::Path;` with no
`PathBuf`, so every `PathBuf` in `main.rs` is spelled `std::path::PathBuf` (§4.3). `io.rs`
already imports `anyhow::{Context, Result}` and `std::path::{Path, PathBuf}`, so §4.1 needs no
new `use`. `preflight_output_collisions` being `pub` in the library avoids a `dead_code`
warning. The five-argument `planned_hardtrim_outputs` is under clippy's `too_many_arguments`
threshold. There is no crate-level `missing_docs` lint and no `[lints]` table, so `fn` → `pub fn`
on an undocumented namer will not trip `-D warnings`.

---

## 9. Assumptions

- **A1 (fixed).** `norm_path`'s ASCII case-fold is the right case treatment. Inherited from
  #216; may false-positive on opt-in case-sensitive APFS volumes, and a loud error is the
  accepted trade.
- **A2 (fixed, narrowed in v2).** *For every secondary output whose directory is resolved by the
  same `output_dir`-else-parent rule as the primary, primary-path distinctness implies
  secondary-path distinctness.* Verified for `report_name`, `json_report_name`,
  `single_end_bam_output_name`, and `--demux` (`demux.rs:142-147` resolves
  `output_dir` else `trimmed_file.parent()`, base name from the trimmed file's own name — so its
  key is a function of the primary). **Named exception:** `--fastqc_args "-o DIR"` /
  `--outdir` (`fastqc.rs:91-96`) overwrites `config.output_dir` with a directory the pre-flight
  never sees, while the artifact filename still derives from the trimmed output. Both reviewers
  reproduced two distinct primaries yielding one `same_trimmed_fastqc.{zip,html}`. Not widened
  into the candidate list: it costs a regenerable QC artifact rather than reads, it is
  pre-existing on every path including `--paired`, and parsing `--fastqc_args` twice would still
  miss the next `FastQCConfig` field. Recorded as §11 Out of scope 5.
- **A3 (fixed).** `--hardtrim5` and `--hardtrim3` never both execute: the `hardtrim5` branch
  returns at `main.rs:367`. Verified empirically. Per-block checks suffice.
- **A4 (fixed).** `gzip` is decided once from `cli.input[0]` (`main.rs:295`) and threaded to both
  the candidate namer and the writer, so prospective and actual paths cannot disagree — and
  argument order cannot produce a false accept (`.gz` first → both `.fq.gz`; `.bgz` first →
  both `.fq`; either way they collide).
- **A5 (fixed).** `cli.input` is non-empty (clap `required`), and `cli.rs:620` rejects
  `--basename` with multiple SE inputs on both output formats.
- **A6 (decision).** The same file listed twice is **rejected**, in `Cli::validate` (§4.4), not
  silently de-duplicated. It exits 0 with correct output today, so this is a behaviour change and
  is logged under `#### Changes`. Rationale: it matches `--paired`; de-duplicating by byte-equal
  path would accept `x.fq.gz x.fq.gz` while still rejecting `x.fq.gz ./x.fq.gz`, an
  inconsistency worse than either clean choice; and de-duplicating properly needs
  `fs::canonicalize`, a syscall per input with different answers across symlinks and bind mounts.
- **A7 (fixed).** Rejection by the pre-flight leaves an empty `--output_dir` behind, matching the
  paired contract that `ci.yml:615`/`:640` already assert. A6's rejection precedes
  `ensure_output_dir`, so it leaves nothing.
- **A8 (fixed, new in v2).** The key is *lexical*. `a/../x` vs `x`, and any symlink alias, still
  hash differently and are still accepted. Pinned by a test (V1.7) so the limitation is visible
  rather than discovered.
- **A9 (fixed, new in v2).** The candidate expression and the writer's expression must stay in
  sync. True by construction for SE (Step 8/9 use calls character-identical to `main.rs:1086-1087`
  and `:1880-1881`) and enforced for hardtrim by the enum (§4.2). This is the assumption whose
  violation is invisible to every *rejection* test, which is why V3 asserts filenames.
- **A10 (fixed).** `std::path::absolute` has POSIX semantics on the CI targets (ubuntu, macos).
  Windows semantics differ (UNC, drive-relative paths); the crate has no Windows CI target and
  none is claimed.

---

## 10. Validation

Every rejection case is paired with an acceptance case **on the same dispatch path and output
format**, and acceptance cases assert *content or filename*, never mere existence — v1's
weakest point, and the exact failure mode SESSION_HANDOFF §5 records four instances of.

**V1 — unit, `src/io.rs`.**
1. distinct outputs → `Ok`.
2. byte-identical outputs → `Err`, message names both.
3. case-only variants → `Err` (unreachable from an integration test: two case-variant paths
   cannot coexist on APFS — this is why the helper is in the library).
4. `./x_trimmed.fq.gz` vs `x_trimmed.fq.gz` → `Err` (defect 1).
5. `<cwd>/x_trimmed.fq.gz` vs `x_trimmed.fq.gz` → `Err` (defect 1).
6. planned path equal to an input path → `Err`, and the message is the **alias** wording, not
   the duplicate wording (defect 2).
7. `a/../x_trimmed.fq.gz` vs `x_trimmed.fq.gz` → `Ok`, pinning A8's residual.
8. `hint: Some(…)` appears in the duplicate message; `None` does not add it.
9. empty `planned` → `Ok`.

**V2 — integration, rejections.** Non-zero exit, stderr contains
`Output path collision (case-insensitive, for APFS/NTFS safety)` (or the alias wording for
2.9–2.10), **and the output directory is empty** — stronger than "the primary is absent", and it
pins the absence of the two misleading reports that made #383 hard to spot.
1. SE trim FASTQ, `sample.fastq.gz` + `sample.fq.gz`.
2. SE trim FASTQ, `./sample.fastq.gz` + `sample.fq.gz` (defect 1).
3. **SE trim FASTQ, `dirA/same.fastq.gz` + `dirB/same.fastq.gz` with `-o <third>`.** The
   highest-value single test in the plan: it is the only case where writing `None` instead of
   `output_dir` in Step 8 silently under-rejects, and every other SE case collides either way.
4. SE trim uBAM, colliding stems.
5. `--hardtrim5 20` FASTQ, same basename in two dirs, `current_dir(tempdir)`.
6. `--hardtrim3 20` FASTQ, same shape.
7. `--hardtrim5 20 --output-format ubam`.
8. `--hardtrim3 20 --output-format ubam`.
9. SE trim, output aliases an input (`s.fastq.gz` + `s_trimmed.fq.gz`) — alias message (defect 2).
10. `--paired`, output aliases an input (the §2.2(f) four-file invocation) — proves the
    converted paired site gained the protection.
11. `.gz` + `.bgz` sharing a stem (#382 regression pin).

**V3 — integration, acceptance.** All assert content or filename.
1. SE trim, `dirA/same.fastq.gz` + `dirB/same.fastq.gz`, **no** `-o` → exit 0, two outputs,
   `dirA`'s containing only `DIRA_` reads and `dirB`'s only `DIRB_`. Existence alone would have
   passed while #383 was live.
2. SE trim, single input → exit 0, one output.
3. `--hardtrim5 20 -o <tmp>`, two distinct stems → exit 0, both `*.20bp_5prime.fq.gz` present.
4. `--hardtrim3 20 -o <tmp>`, two distinct stems → exit 0, both `*.20bp_3prime.fq.gz` present.
   `--hardtrim3` has **zero** output-producing coverage in the repository today
   (`tests/integration_adapter2.rs:377-382` only asserts a `-a2` rejection message), and Step 7
   is a fresh copy of Step 6 into it.
5. `--hardtrim3 20 --output-format ubam`, two distinct stems → exit 0, both `*.20bp_3prime.bam`.
6. SE trim uBAM, two distinct stems → exit 0, both `<stem>_trimmed.bam`. Catches an
   over-rejecting candidate list, e.g. hashing the input path instead of the output.
7. `--paired`, two normal pairs with distinct stems → exit 0, four `_val_` outputs. Regression
   guard on the converted paired site.

**V4 — unit, `src/cli.rs`.** Duplicate SE input rejected with the dedicated message, not the
APFS/NTFS wording; a distinct-input SE command line still validates; the existing duplicate-pair
and R1≠R2 tests still pass. v1 had no test at all for its only intentional behaviour change.

**V5 — unit, A2 as a table.** For path pairs whose primaries differ, assert that
`report_name`, `json_report_name`, the demux base name, and the *default-configuration*
`<stem>_fastqc.zip` all differ too, across `--basename` / `--dont_gzip` / `-o` on and off.
A comment names the `--fastqc_args --outdir` exception so a future reader does not "fix" the
table by deleting a row. This replaces v1's V6, whose premise
(`report_name(a) == report_name(b)` with `a != b`) forced its own conclusion and could
therefore only ever pass.

**V6 — CI validation job**, after `ci.yml:640`, mirroring the existing guards' idiom exactly
(`set +e`, `rc=${PIPESTATUS[0]}`, `set -e`, `test $rc -ne 0`,
`grep -q "Output path collision"`, `test -z "$(ls -A <out>)"`):
1. SE trim collision — `illumina_10K.fastq.gz` copied as `sample.fastq.gz` and `sample.fq.gz`.
2. `--hardtrim5 30` collision — same basename in two subdirectories, with `-o` (not `cd`),
   since `-o` does not rescue this collision (§3.4) and the `ls -A` residue check needs a
   directory to inspect.
3. SE output-aliases-input — a copy named `sample.fastq.gz` plus one named
   `sample_trimmed.fq.gz`.

**V7 — gates and manual re-runs.** `cargo fmt --all -- --check`;
`cargo clippy --all-targets --release -- -D warnings`; `cargo test` from the crate root. Then
rebuild `target/release/trim_galore` (clippy does not refresh it) and re-run all six §2.2
reproductions — (a)–(f) must each exit non-zero with the right message — **plus** the positive
control from V3.1, which is the behaviour most at risk from an over-broad key and which no CI
step covers. Run (c) from a scratch directory, not the repo root.

---

## 11. Questions and ambiguities

**Open 1 — A6's message wording.** Rejection is settled; the exact sentence is not. Proposal:
mirror `cli.rs:535`'s shape — "Input file '<path>' was given more than once." plus one clause on
why that is refused.

**Open 2 — hint mechanism.** `hint: Option<&str>` is the least ceremonious way to fix the
hardtrim advice. Alternative: a small enum of remediation classes, if a reviewer prefers the
message text to live entirely in `io.rs`.

**Open 3 — pre-existing, not fixed.** `--hardtrim5 N --hardtrim3 M` is accepted by
`Cli::validate` but the `hardtrim5` branch returns early, so `--hardtrim3` is silently ignored.
Verified. Its own issue; the fix is a clap `conflicts_with`. Worth filing before this lands,
since a reviewer of these blocks will ask.

**Open 4 — a docs sentence.** `modes/hardtrim.md` advertises multi-file use without noting that
output goes to the CWD, which is exactly why the collision class is wide there. Docs-only,
separate commit, but it should not be lost.

**Out of scope**

1. The CWD-vs-`input.parent()` naming asymmetry (§2.3). The pre-flight is the right response;
   changing the naming moves output locations for every existing hardtrim user.
2. Refusing to clobber files already on disk that are *not* named as inputs — a `--force` /
   no-clobber policy. Distinct from defect 2, which is now in scope: there the user has named
   the file as an input, so writing it is unambiguously wrong and no policy is needed.
3. Symlink and `..` aliasing (A8).
4. Hoisting all collision checks above `ensure_output_dir`; it would change the tested paired
   contract (A7).
5. `--fastqc_args "-o DIR"` collisions (A2's named exception).
6. #384 (uppercase `.FASTQ.GZ`) — it changes what `norm_path` and `strip_fastq_extensions`
   mean, so it should follow this, not precede it.
7. A single mode → candidate-list `match` in `main()` instead of per-branch checks. Both
   reviewers raised it; it is the structural answer to "each new mode grows its own loop and
   three of eight forget the check", but it duplicates the dispatch logic it sits above. The
   helper's call-site enumeration (§4.1) is this plan's cheaper mitigation. Recorded as the
   direction, not adopted.
8. A writer-level registry of already-opened output paths (B's alternative 2). It would catch
   all of these plus the FastQC and symlink cases at one choke point, but fires after the first
   write, so it cannot deliver the "nothing written" contract the CI guards assert.
   Defence-in-depth later; explicitly not a substitute.

---

## 12. Revision history

**v1 → v2**, after dual independent plan review (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`); every
finding below was re-verified by the orchestrator against the binary before adoption.

Adopted from **A**: the absolutised collision key (A's Critical — v1's raw-string key let
`./x` and `x` through, so #383 reproduced through v1's own guard); the helper doc-comment
enumerating call sites; the sanity/format-before-collision ordering note.

Adopted from **B**: the output-aliases-input check (B's Critical — and `--paired` has the hole
today, which is why the four hand-rolled copies convert in this commit rather than later); A6
moving to `Cli::validate`; the hardtrim-specific remediation hint and the `docs/` check; the
`PathBuf` import correction (**A asserted the signatures compile; B was right that they do
not**); the `HardtrimEnd` enum in place of a duplicated `&str`; `--hardtrim3` having no
output-producing coverage at all.

Adopted from **both**: A2 narrowed with `--fastqc_args -o` as a named exception; V1's V6
replaced because it could not fail; acceptance cases added for the three unpaired arms;
V2 upgraded to assert an empty output directory; `current_dir` tests pinned to absolutised
fixture paths; the copy count corrected from five to four.

**Rejected:** B's claim that a bare `trim_galore *.gz` reaches defect 2. Not reproducible —
both zsh and bash collate `s_trimmed.fq.gz` before `s.fastq.gz`, which is the harmless order.
§2.2 states the narrower reachability that was measured.

---

## 13. Implementation notes

Implemented on branch `fix/383-output-collision-preflight` off `dev` @ `72624c7`.
All 13 steps done. Gates: `cargo fmt --check` clean, `cargo clippy --all-targets --release
-- -D warnings` clean, **520 tests pass** (398 lib + 122 integration, 20 of them new).
All eight §2.2 reproductions now exit non-zero; the §10 V3.1 positive control still
succeeds with zero crosstalk and a prior output left byte-intact.

Diff: 6 files changed, +503/−99, plus the new `tests/integration_output_collision.rs`.

### Deviations from the plan

**D1 — Step 10's predicate was wrong (behavioural, caught by existing tests).**
The plan specified `if !self.paired` for the duplicate-input check. `--clock` and
`--implicon` are paired modes that **never set `self.paired`**; they call
`validate_paired_input` separately at `cli.rs:953-956`, i.e. *after* the new check.
The guard therefore pre-empted their more specific "Read 1 and Read 2 appear to be
the same file" message, breaking `test_validate_clock_r1_equal_r2_within_pair_rejected`,
`test_validate_clock_duplicate_pair_rejected` and
`test_validate_implicon_duplicate_pair_rejected`. Predicate is now
`!self.paired && !self.clock && self.implicon.is_none()`. Neither reviewer caught this,
and neither did the plan — the three existing tests did.

**D2 — `guarded_inputs()` helper added.** Not named in the plan. §3.2 required the input
set to be `cli.input` plus `--passthrough`; rather than repeat that at eleven call sites,
one `main.rs` helper builds it. Behaviour as planned.

**D3 — two defect-1 integration tests were passing for the wrong reason.** As first
written they passed `-o`, under which the output path is built from `output_dir` + stem,
so the input's `./` or absolute spelling never reaches the output key and the plain
duplicate check caught them. Found by negative control 1 (below): both passed with the
absolutisation reverted. They now run without `-o`, so the output inherits the input's
spelling, and assert the directory afterwards holds exactly its two inputs.

**D4 — `tempdir()` in the new test file is canonicalised.** `std::env::temp_dir()` yields
`/tmp/…` while the child process's `getcwd` yields `/private/tmp/…` on macOS. The key is
lexical (A8), so it cannot see through that symlink and the absolute-vs-relative test
passed vacuously. This is A8 showing up in the harness rather than a defect in the fix,
and it is now a documented comment on the helper.

**D5 — CI step names truncate at `#383` in YAML** (an unquoted ` #` opens a comment), so
the Actions UI shows "… (issue". The existing `(issue #216)` step has the identical
quirk; left consistent rather than fixing only the new ones. Cosmetic, and a candidate
for a separate one-line sweep.

**D6 — `demux::demux_base_name` extracted (new function, beyond the plan).** Closing the
coverage audit's Gap 2 required asserting that the demux stem is a function of the primary
output path. That derivation was inline in `demultiplex`, so the assertion would have been a
copy of the formula — a test incapable of detecting drift in the code it checks. Extracted as
`pub fn demux_base_name(&Path) -> String` and called from both `demultiplex` and the test.
Pure refactor; verified behaviour-preserving by re-running `--demux` on
`test_files/demux_test.fastq.gz` (stems still `demux_test_trimmed_<sample>.fq.gz`, and the
`Processed sequences from file >…<` line still carries the full filename, which matters because
`--demux` is in the Perl byte-identity matrix). The separate `trimmed_name` binding that feeds
that message was restored after the first extraction attempt dropped it — caught at compile time.

### Post-audit gap closures

The coverage audit (`COVERAGE.md`) returned **INCOMPLETE — 2 items**, both plan shortfalls with
no behavioural effect, both now closed. See the addendum in `COVERAGE.md` for detail.

1. **Step 1's doc-comment sub-clause** — `norm_path` now names `collision_key`, and the stale
   "`main::run` paired-end" locator is gone. Written without rustdoc link brackets to avoid a
   new `private_intra_doc_links` warning.
2. **§10 V5's A2 table** — replaced with `distinct_primary_outputs_imply_distinct_secondary_outputs`,
   which asserts the plan's stated direction over the full `--basename`/`--dont_gzip`/`-o`
   matrix and covers the demux stem via D6; the original coarser-than test is retained as its
   complement. FastQC's zip name is documented as not-separately-assertable rather than faked.

### Post-code-review round (4 review passes, 3 coverage audits)

Coverage: two independent audits returned **COMPLETE** (25 items and 65 items). Code review:
0 Critical, 1 unanimous High. Applied, with the `..` decision taken by the user:

- **D7 — `collision_key` folds `..` lexically** (was A8's documented residual). Absolutise, then
  fold `Component::ParentDir` against the component stack, keeping POSIX `/..` == `/` and
  preserving an unfoldable leading `..`. Reason for the reversal: A8 was scoped out on the
  grounds that `..` is rare, but the reachable shape is the *same* mixed absolute/relative
  argument list §2.2(e) used to justify fixing `./x` — so #383 reproduced verbatim through its
  own fix. **A8 is now narrowed to symlinks only.**
- **D8 — one definition of file identity.** `collision_key` is `pub`; `--paired`'s R1≠R2 check,
  its duplicate-pair check, the duplicate-input check and the `--passthrough` alias check all
  use it instead of raw `==` / `norm_path`. Closes a silent-wrong-output bug: `--paired
  ./a_R1.fq a_R1.fq` exited 0 and wrote two byte-identical files labelled a validated pair.
- **D9 — secondary outputs are checked against inputs.** The unanimous finding: A2 justifies
  checking primaries for output-vs-**output** collisions and says nothing about
  output-vs-**input**, so a report or a `--demux` per-barcode file could overwrite a named input
  while every primary stayed distinct. Reports and demux outputs now join the candidate list;
  `guarded_inputs()` gains the `--demux` barcode file. New `demux::demux_output_paths` shares
  `demux_output_dir` / `output_filename` with the writer, and
  `demux_writes_exactly_the_planned_paths` pins the two together.
- **D10 — the CWD hint covers all four CWD-naming modes and replaces the generic advice.**
  `--clock`/`--implicon` name output into the CWD exactly as hardtrim does, so the generic
  "use different source directories or `--output_dir`" is false for them too. `HARDTRIM_HINT`
  became `CWD_OUTPUT_HINT`, `run_specialty_paired` takes a `hint` (passed for clock/implicon,
  `None` for `--clump_only --paired`, whose namer uses `input.parent()`), and the helper now
  substitutes the hint for the generic clause instead of appending to it.
- **D11 — the V5 fixture can now fail.** Two reviewers found
  `distinct_primary_outputs_imply_distinct_secondary_outputs` unfalsifiable: its inputs had
  pairwise-distinct basenames, and every secondary namer embeds `file_name()`, so the assertion
  held regardless of whether a namer honoured the directory. Confirmed by making `report_name`
  directory-blind and watching the test stay green. Added `d/same.fastq.gz` + `e/same.fastq.gz`;
  the same control now fails. The earlier negative control (constant return) was too crude to
  expose this.
- **D12 — mechanical.** `demux_base_name` lifted out of `demultiplex`'s doc-comment (D6 had
  re-parented it — the same defect Step 2 repaired in `main.rs`); stale key descriptions deleted
  at the paired and `run_specialty_paired` sites; `assert_rejected_cleanly`'s
  `unwrap_or_default()` replaced with a panic so a `read_dir` failure cannot pass vacuously;
  CHANGELOG corrected where it overclaimed defect 2's coverage, plus new entries for the CWD
  hint, the single identity definition, and A1's previously unstated case-fold trade-off.

Gates after the round: **530 tests** pass, `fmt --check` clean, `clippy -D warnings` clean,
`cargo doc` adds no warning. All eight original reproductions plus the two new ones (`..`,
self-mate pair) reject; both acceptance controls still pass with zero crosstalk and the prior
output byte-intact. All three demux routes verified closed individually.

**Still open, unchanged:** §11 Opens 3 and 4; symlink aliasing (A8, now the sole residual);
`--fastqc_args -o` (A2's named exception).

### Iteration log

**#1 — duplicate-input predicate.** Three `--clock`/`--implicon` validation tests failed
after Step 10. Excluded both modes from the new check (D1). All 398 lib tests green.

**#2 — defect-1 test rigour.** Negative control 1 (revert `collision_key` to
`norm_path` only) showed the unit test failing but both integration tests still passing.
Removed `-o` from those two cases and switched to an exact-directory-contents assertion;
they then failed under the control and passed once it was reverted (D3).

**#3 — symlink-induced vacuous pass.** After #2 the absolute-vs-relative test failed with
the fix in place, because the child's `getcwd` resolved `/tmp` → `/private/tmp` while the
test's absolute path did not. Canonicalised `tempdir()` (D4). 20/20 green.

### Negative controls run (per SESSION_HANDOFF §5)

| Control | Expectation | Result |
|---|---|---|
| `collision_key` → `norm_path` only | defect-1 tests fail | 1 unit + (after D3) 2 integration failed; reverted → all pass |
| input-key set emptied | defect-2 tests fail | 1 unit + 2 integration failed; reverted → all pass |
| three CI step bodies run locally on `illumina_10K.fastq.gz` | all three pass | all three PASS |
| `ci.yml` parsed | valid YAML | 40 steps in `validation` |

### Left undone, deliberately

§11 Opens 3 (`--hardtrim5` + `--hardtrim3` silently ignoring the 3′ trim) and 4 (the
`modes/hardtrim.md` sentence about CWD output) are unaddressed, as scoped. Both want
their own issue.
