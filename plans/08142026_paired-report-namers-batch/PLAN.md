# PLAN — #423 → #424 → #429, one PR

**Issues:** [#423](https://github.com/FelixKrueger/TrimGalore/issues/423) (namers),
[#424](https://github.com/FelixKrueger/TrimGalore/issues/424) (pre-flight for two arms),
[#429](https://github.com/FelixKrueger/TrimGalore/issues/429) (opaque write-side refusal)
**Branch base:** `dev` @ `d2ea293`
**Test baseline:** 707 macOS / 706 Linux (one linux-only test in `cli.rs`, two macOS-only in
`format.rs`; do not gate on a single number)
**Version:** v3 — dual plan review folded in (`PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md`); O1 decided

---

## 1. Goal

Three items, one PR, three commits:

1. **#423** — dedicated namers for the paired trimming-report paths, so the writer and the
   pre-flight candidate list cannot disagree. Pure refactor: **no path changes.**
2. **#424** — the output-collision pre-flight for the two single-file interleaved paired arms,
   the only write-capable dispatch paths without one.
3. **#429** — make the uBAM writer's refusal name the offending read name.

**On ordering:** v1 said "the order is not a preference". Both reviewers corrected that. #429 is
independent (it touches only `bam.rs`). #424 depends on #423 for *tidiness*, not correctness — it
could hand-roll the strings. So 423 → 424 is a quality dependency and 429 is free-floating. The
order stands because it keeps `main.rs` shorter and gives both arms one place to look, not
because anything breaks otherwise.

---

## 2. Context

### 2.1 #429 has narrowed since it was filed

Its body pairs the opaque failure with a partial output file. #428 landed as `d2ea293`, so the
residue half is gone — both reviewers confirmed empirically that a refused run now leaves only
the input. What remains is the message:

```
Error: failed to write BAM record
Caused by: invalid input parameter
```

### 2.2 The report-path construction sites

Corrected after review. `pair_output_dir` has **eight** callers outside `io.rs` — six in
`main.rs` plus `clump_only.rs:537` and `:1089` — and `io.rs` uses it internally seven times.
Four `main.rs` sites are in scope:

| Site | Kind | Shape |
|---|---|---|
| `main.rs:922-923` | **candidate** — FASTQ trim paired | per-mate `report_name`/`json_report_name`, `OutputSource::Input` |
| `main.rs:1742-1748` | **writer** — `run_paired` | per-mate `PairedReportFile` |
| `main.rs:2040-2041` | **candidate** — uBAM trim paired | same loop shape as `:922` |
| `main.rs:2367-2373` | **writer** — `run_ubam_output_paired_two_files` | same as `:1742` |

Out of scope, and the descriptions v1 had **swapped** (both reviewers caught this
independently):

- `main.rs:619` passes `&[r1, r2]` to `clump_report_candidates` → **two** clumping reports, one
  per mate, matching `clump_only_paired`.
- `main.rs:746` passes `std::slice::from_ref(&chunk[0])` → **one** report keyed on the pair's
  first input, matching `clump_only_paired_to_bam_one_pair`.
- `clump_only.rs:537`, `:1089` — the clump writers, same different-scheme bucket.

They are *not* the same shape as each other, which matters if anyone later tries to fold them
into one helper.

No seventh in-scope site: `write_paired_reports` constructs no paths (it takes
`PairedReportFile` descriptors), and SE reports call `report_name(input, output_dir)` directly.

Plus the two hand-rolled blocks, produced by no namer and used nowhere else:
`main.rs:1959-1966` (FASTQ out from interleaved uBAM) and `:2489-2496` (uBAM out).

### 2.3 The naming divergence that makes this dangerous

**`io.rs::report_name` keys on `file_name()`; both interleaved arms key on `file_stem()`.**

```
report_name(Path::new("sample.fq.gz"))  →  sample.fq.gz_trimming_report.txt
the arms' stem for "ip.bam"             →  ip_R1_trimming_report.txt
```

Divergent by design. So the natural implementation — put `interleaved_report_names` in `io.rs`
beside `report_name` and compose it from `report_name` — silently renames every interleaved
report to `ip.bam_R1_trimming_report.txt`, breaking MultiQC/nf-core glob expectations.

**And nothing in the repo would catch it.** `grep -rn "_R1_trimming_report" tests/ docs/
.github/` → zero hits. `integration_ubam.rs` pins the interleaved *primaries* against committed
references and says nothing about reports. Worse, after step 5 the candidate list and the writer
both derive from the same namer, so a joint rename keeps `created ⊆ planned` true and every
promoted tripwire case green. This is the single riskiest edit in the batch; §7 V1b exists for
it alone.

### 2.4 Candidate order, and what actually pins it

Order at both candidate sites is, per mate, txt then json — verified by both reviewers running
`TRIM_GALORE_DUMP_PLANNED_OUTPUTS=1`:

```
src/main.rs:928   a/r1_val_1.fq            a/r1.fq_trimming_report.txt
src/main.rs:928   a/r2_val_2.fq            a/r1.fq_trimming_report.json
src/main.rs:928   a/r1_unpaired_1.fq       a/r2.fq_trimming_report.txt
src/main.rs:928   a/r2_unpaired_2.fq       a/r2.fq_trimming_report.json
```

`integration_output_collision.rs:1040-1043` asserts a refusal names
`reads.fq_trimming_report.txt`, so **txt-before-json within a mate is pinned** — swap them and
the collision fires on `.json` and that test goes red. **Mate order is pinned nowhere**, and
neither is reports-after-primaries. #414's tripwire compares sets (`planned_keys:
BTreeSet<String>`), so it sees no reorder at all.

So the shape that would actually break a test is grouping by *kind* (`[all txt, all json]`),
not returning four paths per call. The per-mate signature is still right, but for a better
reason than v1 gave (§4 A2): both **writer** sites build a `PairedReportFile` per mate, so
per-mate is the only shape one function can serve all four sites with.

### 2.5 #424's premise, verified

Neither interleaved arm reads `cli.basename` — `run_paired_ubam_single_file` (1807-1988) and
`run_ubam_output_paired_single_file` (2411-2529) contain no occurrence. Both derive the stem with
`file_stem()`, which strips one `.`-suffix, while every appended suffix begins with `_`, so
`stem + suffix` can never equal `file_name`. Reviewer B confirmed this survives the `.gz` case:
`one.bam.gz` → `one.bam_val_1.fq.gz`.

**Consequence, and it changes §7:** a collision on these arms is **unreachable by construction**.
Output-vs-output is impossible (`_val_1`/`_val_2`, `_R1`/`_R2` all distinct under case-folding);
output-vs-input is impossible by the argument above; `guarded_inputs` is a singleton because
`--passthrough` needs two inputs and `--demux` is SE-only. The pre-flight is pure
defence-in-depth against the day something honours `--basename` there — which is the point, but
it means **no test can provoke the refusal**, and v1's V3 asked for exactly that.

### 2.6 Write sets, verified empirically by both reviewers

| Arm | Files written |
|---|---|
| `--paired one.bam` | `one_val_{1,2}.fq`, `one_R{1,2}_trimming_report.{txt,json}` |
| `+ --retain_unpaired` | above `+ one_unpaired_{1,2}.fq` |
| `+ --cores 4` | identical set |
| `+ --no_report_file` | primaries only |
| `--paired --output-format ubam one.bam` | `one_val.bam`, `one_R{1,2}_trimming_report.{txt,json}` |
| `+ --fastqc` | `+ one_val_fastqc.{html,zip}` — covered by the tripwire's FastQC allowance |
| **`--paired one.bam.gz`** | **`one.bam_val_{1,2}.fq.gz`, `one.bam_R{1,2}_trimming_report.*`** |

The last row is load-bearing: `gzip` is `!dont_gzip && is_gzipped(&cli.input[0])` and
`is_gzipped` is extension-only, so a `.gz`-*named* BAM produces `.fq.gz` outputs. **A candidate
list that hardcodes `.fq` passes all three promoted cases and mis-plans every `.gz` run** — it
would guard paths nobody writes while the real outputs go unguarded. Hence V5's extra case.

`--fastqc` is refused on the FASTQ interleaved arm (#421) and accepted on the uBAM one.
`--clumpify` is refused earlier (`--cores >= 2`), and at `--cores 2` the write set is unchanged.

### 2.7 What #424 will deliberately break

Three tripwire cases carry `.unplanned()` and assert that *nothing is dumped*:
`orphan_paired_fastq_from_interleaved_bam` (`:1015`),
`orphan_paired_fastq_retain_unpaired` (`:1029`),
`orphan_paired_ubam_from_interleaved_bam` (`:1040`).

**Five tests will go red, not three** — see §5 step 6. Both reviewers flagged that the extra two
would otherwise be misread as the promote signal misfiring.

---

## 3. Behaviour

**#423.** No output path changes anywhere: same paths, same order, same `OutputSource`
provenance, same `PairedReportFile` shape.

**#424.** Both arms gain a pre-flight before any I/O. Provenance is **`OutputSource::Input(input)`
for every candidate, primaries included** — see §4 O1. A refusal is unreachable today (§2.5).

**#429.** A read name the encoder rejects produces an error naming the offending name. The
accepted set does **not** change: #415 deliberately admits `@`, non-ASCII, VT and >254-byte names
read-side because none defeats the FASTQ id's splitters, and both reviewers confirmed all three
still run clean to `*_trimmed.fq` on the FASTQ path. This is a message fix at the write site.

### Edge cases

| Case | Behaviour |
|---|---|
| `--no_report_file` on either interleaved arm | report candidates omitted, primaries still planned |
| `--retain_unpaired` on the FASTQ interleaved arm | two more candidates; absent on the uBAM arm, which rejects the flag |
| A `.gz`-named BAM input | `.fq.gz` primaries and the same report names; candidates must take `ext` from the arm's own binding |
| A collision on either interleaved arm | unreachable by construction (§2.5); the machinery is there for the day `--basename` is honoured |
| A read name of exactly 254 bytes (uBAM out) | accepted; 255 refused |
| `@` in a read name | refused on uBAM out, accepted on FASTQ out — unchanged |

---

## 4. Signatures, and the one open scope decision

```rust
/// Text + JSON trimming-report paths for one mate of a two-file pair.
///
/// Both land in the pair's output directory — R1's parent unless `--output_dir`
/// overrides — while the filename stays keyed on the mate's own input.
pub fn paired_report_names(
    input: &Path, input_r1: &Path, output_dir: Option<&Path>,
) -> (PathBuf, PathBuf);
```

**O1 — decided (Felix, 2026-08-14): the struct.** Report-only namers would leave #424
hand-rolling the *primaries* (`format!("{}_val_1{}", stem, ext)` at `:1826-1827`,
`_unpaired_{1,2}` at `:1835-1836`, `_val.bam` at `:2430`) — so the PR whose first commit deletes
two copies of the report strings would add two copies of the primary strings. Both reviewers
converged on this shape independently.

The two arms write different primaries, so the reports are shared and the primaries are not:

```rust
/// The two mates' report paths for a single interleaved pair. `_R1`/`_R2`
/// infix because both derive from one input.
pub struct InterleavedReports {
    pub r1_txt: PathBuf, pub r1_json: PathBuf,
    pub r2_txt: PathBuf, pub r2_json: PathBuf,
}

/// Every path a FASTQ-output single-file interleaved paired run writes.
pub struct InterleavedFastqOutputs {
    pub val_1: PathBuf,
    pub val_2: PathBuf,
    /// `Some` under `--retain_unpaired`.
    pub unpaired: Option<(PathBuf, PathBuf)>,
    /// `None` under `--no_report_file`.
    pub reports: Option<InterleavedReports>,
}

/// Every path a uBAM-output single-file interleaved paired run writes.
pub struct InterleavedBamOutputs {
    pub val: PathBuf,
    pub reports: Option<InterleavedReports>,
}
```

Two properties this buys that report-only namers do not:

1. **The conditionals move into the type.** `--retain_unpaired` and `--no_report_file` are
   `Option` rather than an `if` repeated at each site, so a candidate list cannot disagree with a
   writer about whether a file exists.
2. **A forgotten field becomes a build failure.** Each type gets

   ```rust
   pub fn planned(&self, src: &OutputSource) -> Vec<PlannedOutput>
   ```

   whose body opens with an **exhaustive destructure** — `let Self { val_1, val_2, unpaired,
   reports } = self;` — so adding a field later cannot compile until it is either enumerated or
   explicitly discarded. That is the mechanism that earns #423's phrase; a `Vec` assembled by
   hand from `self.foo` accessors does not.

Also kept, for the two-file paired arms (four sites, §2.2):

```rust
/// Text + JSON trimming-report paths for one mate of a two-file pair.
///
/// Both land in the pair's output directory — R1's parent unless `--output_dir`
/// overrides — while the filename stays keyed on the mate's own input.
pub fn paired_report_names(
    input: &Path, input_r1: &Path, output_dir: Option<&Path>,
) -> (PathBuf, PathBuf);
```

Per-mate, because both **writer** sites build a `PairedReportFile` per mate — the only shape one
function serves all four sites with (§2.4).

---

## 5. Implementation outline

**Step 0 — capture the baselines.** Before touching anything: the created-file set for all seven
rows of §2.6, and the ordered candidate dumps for the two candidate sites. These are V1's
before-halves and cannot be reconstructed afterwards.

**Step 1 — `io.rs`: `paired_report_names`,** with the "both mates share R1's directory; each
filename stays keyed on its own input" invariant above the function.

**Step 2 — thread it through the four sites**, writers first (`:1742-1748`, `:2367-2373`) then
candidates (`:922-923`, `:2040-2041`), so the order-sensitive edits come last. **Remove the
orphaned `pair_dir` bindings** at `:883`, `:1740`, `:2037`, `:2365` — `unused_variables` is a
hard error under `-D warnings`, and at `:883` the binding sits 39 lines above its only uses.

**Step 3 — `io.rs`: the three interleaved types and their two constructors.** State on
`InterleavedReports` that the paths key on **`file_stem()`**, not `file_name()`, and why (§2.3).
Do **not** compose them from `report_name`. Each output type gets `planned(&self, src)` opening
with an exhaustive `let Self { .. } = self;` destructure (§4). While here, route both arms' `dir`
through `naming::pair_output_dir` — they each inline a copy of its body at `:1823-1825` and
`:2427-2429`.

**Step 4 — replace the two hand-rolled blocks** at `:1959-1966` and `:2489-2496`, and the
hand-rolled primaries at `:1826-1827`, `:1835-1836` and `:2430`, so each arm reads its paths from
the struct it was handed.

**Step 5 — #424's pre-flight.** Two structural points the plan v1 missed:

- **Hoist.** The report paths are currently computed *after* the whole trimming run
  (`:1951-1968`, `:2481-2498`). A pre-flight before any I/O needs them earlier, so compute them
  above the reader/writer setup and thread them down into the `if !cli.no_report_file` block.
- **Site.** Put the pre-flight at the **dispatch site** (`main.rs` ~`:853`, ~`:2004`) *before*
  `setup_trimming`, not at the top of the arm function. `setup_trimming` is adapter
  auto-detection over up to 1 M reads; a refusal after that has made the user wait for nothing.
  `cli`, `output_dir` and `gzip` are all in scope there. Build the struct at the dispatch site,
  pre-flight from `planned()`, and pass it into the arm — so the arm never re-derives a path and
  the hoist in the first bullet falls out for free.
- **Provenance: `OutputSource::Input(input)` for every candidate.** Not `Pair` — see O1 note
  below.
- **The candidate list is `planned()`**, never a second derivation of stem/ext/dir.

**Step 6 — the tripwire, which is five edits and one decision:**

1. Drop `.unplanned()` at `:1021`, `:1035`, `:1046`.
2. **Delete `fn unplanned` (`:267-270`)** — it is the only constructor of `Arm::Unplanned`, so
   leaving it turns `dead_code` into a CI failure under `-D warnings`. **Keep** the variant and
   `compare`'s early-return branch: they are the trap for a future no-pre-flight arm. Add one
   self-test driving `compare(…, Arm::Unplanned)` with non-empty `planned`, expecting the
   "promote the case" error — otherwise that branch is reachable by no test.
3. `the_site_inventory_is_complete` (`:1149-1164`): **11 → 13**. Its own message says *"do not
   just update the count"*, so items 4 and 5 come first.
4. `every_site_is_reached_by_some_case` (`:1253-1334`): add `["--paired", "ip.bam"]` and
   `["--paired", "--output-format", "ubam", "ip.bam"]`, and move both the array size and
   `assert_eq!(sites.len(), …)` to 13. `ip.bam` is already staged for every shape.
5. Four stale doc spots: module Limits (`:63-65`), the `Arm::Unplanned` variant doc
   (`:214-221`), the section header at `:1012` ("the two arms with no pre-flight"), and T3's doc
   (`:1116-1118`).

**Step 7 — #429.** The name **is** reachable at the wrap site, so no predicate duplication:
`BamWriter::write_record` holds it, and the failing call is
`write_alignment_record(...).context("failed to write BAM record")` (`bam.rs:641-643`).

- **Capture `record.id`, not `name`** — `name` is moved into `set_name` on the line above. This
  matches the precedent five lines up: `format!("while writing BAM record {:?}", record.id)`.
- **No ordinal.** `BamWriter` has no counter and the cited precedent
  (`reject_whitespace_in_qname`) has none either — it names the escaped name, the reason and a
  remedy. #429 asks for the name. Dropped from the plan.
- **Two distinct underlying errors, which decides the shape.** `@` gives
  `invalid input parameter`; a 255-byte name gives `out of range integral type conversion
  attempted`. So a pre-check would have to reproduce both the charset *and* the length rule —
  O2 closes in favour of wrapping.
- **Two-tier context, to avoid misdiagnosing I/O faults.** `write_alignment_record` also fails on
  disk-full and EPIPE, and unconditionally appending "read names must be…" would relabel a full
  disk as a read-name problem. Reviewer A proposed gating on `ErrorKind::InvalidInput`; that
  alone would leave the 255-byte case opaque, since a wrapped `TryFromIntError` does not land in
  `InvalidInput`. So: **always** attach the neutral half (`while writing BAM record {id:?}`) —
  true for every failure and already enough to fix #429's complaint — and attach the SAM-QNAME
  explanation only for the kinds that actually carry it. **Determine both kinds by probe as the
  first implementation action** (O3).

**Step 8 — docs + `CHANGELOG.md`.** #423/#424 are internal; #429 changes a user-visible message.
One changelog sentence.

---

## 6. Efficiency

Two `parent()` lookups per pair instead of one, twice per run — immeasurable. Do not take the
"pass the already-computed `pair_dir`" shortcut: taking `input_r1` is what makes the derivation
unrepresentable-if-wrong. #424 adds one hash of ≤8 paths.

The one user-visible latency point is step 5's siting: pre-flighting after `setup_trimming` would
make a refused run pay for a full adapter-detection scan first. With `InterleavedOutputs` built at
the dispatch site, #424 is *smaller* than v1's version, not larger — no second derivation.

---

## 7. Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| **V1a** | #423 changes no **candidate** path | diff the step-0 ordered dumps against post-refactor. **Compare field 3 only** (`… 2>&1 \| grep TRIM_GALORE_PLANNED \| cut -f3`) — field 2 is a `#[track_caller]` site whose line number the refactor shifts, so a byte diff is all noise. Env var is `TRIM_GALORE_DUMP_PLANNED_OUTPUTS` | identical, order preserved |
| **V1b** | **#423 changes no *interleaved* report path** — the gap both reviewers called Critical | (a) `io.rs` unit tests pinning both constructors' literal output: the `file_stem` convention, `--output_dir` override, a bare filename, a two-dot name (`ip.bam.gz`), and both `Option` states; (b) diff the step-0 **created-file sets** for all seven §2.6 rows | identical filenames. V1a is structurally blind here — those arms dump nothing today — and after step 5 a typo is self-consistent across writer and candidate list |
| **V1c** | `planned()` enumerates every field | a unit test asserting `planned()`'s length and contents for each `Option` combination, plus the exhaustive destructure making a new field a **build** error | a field added later cannot silently escape the candidate list |
| **V2** | Collision messages unchanged | `cargo test --release --test integration_output_collision` | 57/57. Pins txt-before-json within a mate, and nothing else (§2.4) — V1a's ordered diff is the real order guard |
| **V3** | #424's candidate lists have the right **content** | `DUMP_PLANNED_ENV` on both arms; assert all expected paths appear, including the `.gz` row | v1's V3 ("provoke a refusal") was unconstructible — §2.5 proves the refusal unreachable. This replaces it |
| **V3b** | The one-input refusal message reads correctly | `io.rs` unit test driving `preflight_output_collisions` with a synthetic one-input candidate list | names one source, no `X + X`. Reviewer B measured the real thing on a sibling arm |
| **V4** | #424 does not over-reject | the adversarial names both reviewers verified: `x_val_1.fq`, `y_R1_trimming_report.txt` | run succeeds, input survives, as today |
| **V5** | The promoted cases go red then green | run after step 5, before step 6 | **five** red: the three promoted plus `the_site_inventory_is_complete` and `every_site_is_reached_by_some_case`. Plus a new `--paired ip.bam.gz` case, so the candidate list cannot hardcode `.fq` |
| **V6** | #429 names the offending read name | uBAM-output runs with `@bad@name` and a 255-byte name, generated in-test | the message **names the read name**. Not "no `invalid input parameter`" — anyhow prints the whole cause chain, so a wrap keeps it as the cause; v1's expectation contradicted its own chosen approach |
| **V7** | #429 over-rejects nothing | the same names on the **FASTQ** path | succeed, as #415 intends. Both reviewers confirmed this already holds |
| **V8** | Suite green both platforms | `cargo test --release` | 707 macOS / 706 Linux plus new cases |
| **V9** | CI by name | `Validate vs Perl TrimGalore` + `Validate uBAM input` | green |

---

## 8. Assumptions

1. **Fixed:** the clump sites are a different scheme and stay out (§2.2, descriptions corrected).
2. **Fixed:** neither interleaved arm reads `cli.basename` (§2.5, verified by both reviewers).
3. **Fixed:** candidate order is pinned only as txt-before-json-within-a-mate; the tripwire
   compares sets (§2.4).
4. **Fixed:** #429 is a message fix; the accepted read-name set does not change.
5. **Fixed:** `interleaved_report_names` keys on `file_stem()`, not `file_name()` — the difference
   between a refactor and a silent rename (§2.3).
6. **Resolved, not pending:** the one-input refusal message reads correctly **under `Input`
   provenance**. Measured. It would read `one.bam + one.bam` under `Pair`.
7. **Fixed:** a collision is unreachable on both interleaved arms today (§2.5), so #424 is
   defence-in-depth and no test can provoke its refusal.
8. **Fixed:** `Arm::Planned` needs no machinery the promoted cases lack — planned non-empty,
   single site, all candidates inside the case root. The blockers are the file's own counters.
9. **Assumed:** #423 is a pure refactor. V1a **and V1b** are what make that checkable; V1a alone
   is not.

---

## 9. Notes from review — not filed

Per the filing bar, none has a user-visible symptom worth an issue; all are inside this PR.

- **Both arms inline a copy of `pair_output_dir`'s body** (`:1823-1825`, `:2427-2429`). Step 3
  routes them through it.
- **`write_paired_reports`' doc names two callers; there are four** (`main.rs:2533-2536`). One
  line while step 2 is touching all four.
- **Tripwire section headers cite site line numbers** (`── site 1006 ──`,
  `── sites 2047 / 2109 ──`); both shift under this batch.
- **`--gzip` (hidden, `cli.rs:424`) is not consulted** by the gzip resolution, which is
  extension-only. Pre-existing, plausibly deliberate Perl-era compat, out of scope — noted so
  §2.6's `.gz` row is not mistaken for it.
- **The `hint` argument for the two new pre-flight calls is unobservable** — only the
  duplicate-path messages use it and those are unreachable with one input. `Some(PAIRED_REPORT_HINT)`
  matches the siblings; its text would read oddly if it ever fired.
- **"Drift becomes unrepresentable" is an overstatement** for report-only namers: both sides
  still call the namer independently, so passing the wrong `input` remains possible. Only O1's
  struct earns the phrase.

---

## 10. Questions

- **O1 — answered (Felix, 2026-08-14): the struct.** §4 carries the shape. The diff is wider than
  #423's body describes, and the two behaviour-free properties it buys — conditionals in the type,
  a forgotten field as a build error — are what make the issue's "unrepresentable" literal.
- **O2 — closed in favour of wrapping.** Two distinct noodles failure modes mean a pre-check
  would duplicate a charset *and* a length rule, and the name is reachable at the wrap site.
- **O3 — Open, resolve by probe as the first implementation action:** which `io::ErrorKind` does
  each of the two failures carry? The answer decides which errors get the SAM-QNAME explanation
  appended. The neutral half (`while writing BAM record {id:?}`) is attached unconditionally
  either way, so #429 is satisfied regardless — this only affects how much the message explains.
- **O4 — Open, assumption taken: one PR, three commits**, so the pure refactor stays bisectable.
  Reviewer B noted #429 could ship independently in an hour; Felix chose one PR.

---

## 11. Self-review

**What v1 got wrong, and what caught it.** Four things, none of which I would have found by
re-reading: step 5 prescribed `OutputSource::Pair` on a one-input arm, which renders
`one.bam + one.bam` (B produced the string rather than inferring it); V1's whole safety case was
blind to the two sites steps 3-4 rewrite, which are also pinned by nothing else in the repo;
V3 asked for an input that §2.4 of the same plan proves cannot exist; and step 6 was one edit
where it is five plus a `-D warnings` failure.

**The thing most likely to go wrong now.** Composing the interleaved namer from `report_name`.
It is the natural move, it silently renames a user-visible file, and after step 5 both sides
agree with each other about the wrong name. V1b is the only guard.

**Second most likely.** A candidate list that hardcodes `.fq`. All three promoted cases pass; every
`.gz`-named-input run is mis-planned. Hence §2.6's last row and V5's extra case.

**Remaining risks.** (a) O3 is a probe, not a design risk. (b) `cargo llvm-cov` coverage of
`io.rs` will rise while `main.rs`'s falls as logic moves; not a gate, but the Coverage job's delta
will look odd. (c) The struct makes the diff wider than either issue body describes, so the
coverage audit should read §4 before judging scope.

---

## 12. Implementation notes

### #423 → `6569652` (previous session)

Steps 1–4 as written. 707 → 714 tests.

### #424 (this session)

Steps 5–6 as written, with three deviations.

**D1 — V5 showed four red, not five.** `the_site_inventory_is_complete` went red (it reads
`main.rs` and counts `preflight_output_collisions(` — 11 → 13). `every_site_is_reached_by_some_case`
did **not**, and cannot: it counts *distinct sites dumped by its own hand-listed shapes*, and
none of the 11 shapes was `--paired ip.bam`. Giving an unexercised arm a pre-flight is invisible
to it. Both plan reviewers predicted five; the correct number is four. Its real role is the one
the tripwire's own message describes — it goes red on a **count bumped without a shape added**,
which is the failure mode step 6 item 3 was ordered to avoid. The three promoted cases failed with
the intended text: *"this arm had no pre-flight, but one dumped 5 candidates from src/main.rs:2023
— promote the case to Arm::Planned"*.

**D2 — the three promoted cases were renamed.** `orphan_paired_*` → `interleaved_paired_*`;
"orphan" named the no-pre-flight state that no longer exists. Only plan artefacts referenced the
old names.

**D3 — §9's `write_paired_reports` doc fix was still outstanding** (both its own doc and
`PairedReportFile`'s said "two callers"; there are four). Done here, since #424 adds the second
single-file caller.

**Verification.** V3 dumped both arms in every `Option` state plus `-o out`; the `.gz` row
confirmed live as `ip.bam_val_1.fq.gz`, so the new tripwire case is a real guard and not a
vacuous one. V4: BAMs named `x_val_1.fq` and `y_R1_trimming_report.txt` both ran to exit 0 with
inputs byte-intact. V1c: `planned()` asserted over all four FASTQ `Option` combinations and both
BAM ones. V3b: `interleaved_refusal_names_one_input` shows the refusal reads *"named from ip.bam"*
with no `ip.bam + ip.bam`. Suite 714 → **719** macOS (+3 `io.rs`, +2 tripwire), fmt and clippy
clean.

**Harness note.** The first V3 dump returned six empty sections: `for args in "--paired ip.bam"`
does not word-split in **zsh**, so the whole string arrived as one argument and every run failed
silently. The `feedback_harness_false_negatives` memory names this exact shape.

### #429 (this session)

**O3 answered, and it overturns the plan's premise.** Probed with a throwaway
`tests/tmp_o3_probe.rs` (deleted after reading):

| Name | Outcome | Cause text | `io::ErrorKind` |
|---|---|---|---|
| `@bad@name` | refused | `invalid input parameter` | **`InvalidInput`** |
| 254 bytes | **accepted** | — | — |
| 255 bytes | refused | `out of range integral type conversion attempted` | **`InvalidInput`** |
| `réad` | refused | `invalid input parameter` | `InvalidInput` |
| `rea\x0bd` (VT) | refused | `invalid input parameter` | `InvalidInput` |

Step 7's third bullet assumed *"a wrapped `TryFromIntError` does not land in `InvalidInput`"*, and
built the two-tier design partly to work around that. It is false — noodles wraps it as
`InvalidInput` too, so **one** kind gate covers both read-name failures. Reviewer A's original
`InvalidInput` suggestion was right, and the plan's objection to it rested on an unverified
assumption. The two-tier shape still ships, because it is what keeps a full disk or a closed pipe
from being relabelled a read-name problem — but it is now one `if`, not a kind list.

**D4 — the shown name is capped at 60 chars.** A FASTQ id has no length limit, so the neutral half
echoed an unbounded untrusted string; the 255-byte case printed 256 characters on one line. `bam.rs`
already caps shown values at 60 with `…` in two places (the aux-tag-framing refusal and the
dropped-description notice), so this follows that. The byte count is appended when the cap bites,
since for the length rule the length *is* the diagnosis.

**Verification.** V6 asserts the exact shown form, truncation included, on both failures; it was
proved able to fail by planting the old `"failed to write BAM record"` string and watching it go
red. V7: `@`, 255 bytes, non-ASCII and VT all still run to `x_trimmed.fq` on the FASTQ path, and a
254-byte name is accepted on the uBAM path, so the refusal is the length rule and not an
off-by-one. Suite 719 → **722** macOS.

**Harness note.** `record.id` is used rather than `name`, which is moved into `set_name` on the
line above — as the plan said, and matching the `validate_and_normalize_seq_for_write` precedent
five lines up.

---

## 13. Revision history

**v3 (2026-08-14)** — O1 answered: the struct. §4 replaced with `InterleavedReports` +
`InterleavedFastqOutputs` + `InterleavedBamOutputs`, the conditionals moved into `Option` fields,
and `planned()` specified to open with an exhaustive destructure so a later field cannot escape
the candidate list without a build error. Steps 3-5 committed to that shape; step 4 now also
replaces the hand-rolled primaries. New V1c.

**v2 (2026-08-14)** — after dual plan review. Both reviewers returned "not safe to implement as
written", with four Criticals each and nine findings in common. Substantive changes: provenance
corrected to `OutputSource::Input` (v1 actively prescribed a message reading `X + X`); V1 split
into V1a/V1b after both found the interleaved half uncovered by anything; V3 replaced, having
been unconstructible; step 6 expanded from one edit to five plus the `fn unplanned` deletion that
`-D warnings` forces; step 5 gained the report-path hoist and moved the pre-flight to the dispatch
site so a refusal precedes adapter detection; step 7 lost the nonexistent record ordinal, gained
`record.id` (since `name` is moved) and a two-tier context after I found that A's suggested
`InvalidInput` gate would leave the 255-byte case opaque; §2.6 gained the `.gz` row; §2.2's clump
descriptions un-swapped and the site count corrected to eight; the "hard order" claim softened to
a quality dependency; and O1 opened as a real scope decision.

**v1 (2026-08-14)** — initial plan.
