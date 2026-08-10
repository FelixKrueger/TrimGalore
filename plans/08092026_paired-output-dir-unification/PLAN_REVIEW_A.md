# Plan Review A — #398 one output directory per pair

**Reviewer:** A (independent; no shared state with Reviewer B)
**Target:** `plans/08092026_paired-output-dir-unification/PLAN.md` r1
**Tree:** `dev` @ `38deb1ee`, working tree clean of tracked changes
**Verdict:** **REVISE** (three Critical findings; the plan's shape is sound)

Verification legend: **[CONFIRMED]** — I read or ran the specific thing named, output
quoted or summarised. **[DERIVED]** — followed from code I read plus the plan's own
steps, but not executed (a patched build does not exist). **[SUSPECTED]** — reasoned
only.

A note on the binary I used as an oracle: `target/release/trim_galore` was stamped
`3f1b0fe` when I started, which is **not** an ancestor of HEAD — it is the pre-squash
commit of the same #408 change, and `git diff 3f1b0fe 38deb1e -- src/main.rs src/cli.rs`
is confined to the `--rename`/uBAM guard, touching no naming or pre-flight code
**[CONFIRMED]**. I then rebuilt via `cargo test --release` and re-ran the load-bearing
repro against a binary stamped `38deb1e`; results were identical. Validation 9's
"unpatched build" control must be rebuilt from HEAD rather than reusing whatever is
in `target/`.

---

## 1. Logic review

### 1.1 The plan's core claims all hold

I verified every line-numbered assertion in the plan. All of them are correct:

| Plan claim | Verdict |
|---|---|
| Five primary paired namers share one directory expression, `io.rs:296-298 / 337-339 / 368-370 / 473-475 / 522-524` | **[CONFIRMED]** identical `output_dir.map(\|d\| d.to_path_buf()).unwrap_or_else(\|\| input_r1.parent().unwrap_or(Path::new(".")).to_path_buf())` in all five; `_input_r2` unused in the two BAM ones |
| `report_name` `io.rs:551`, `json_report_name` `io.rs:567`, `clumping_report_name` `io.rs:535` derive filename and directory independently | **[CONFIRMED]** |
| `passthrough_output_name` `io.rs:404-409` anchors on the carrier's own parent; `io.rs:385` already calls it one of "the three pair outputs" | **[CONFIRMED]** |
| Writers: `main.rs:1729-1738`, `main.rs:2352-2358`, `clump_only.rs:535-536`, `clump_only.rs:1084` | **[CONFIRMED]** all four exist as described |
| Candidates: `main.rs:906-912`, `main.rs:2026-2027`, `main.rs:893-903`, `clump_report_candidates` `main.rs:121-139` | **[CONFIRMED]** |
| `PAIRED_REPORT_HINT` at `main.rs:57-62` | **[CONFIRMED]**, and it is shared by three call sites (`main.rs:918`, `main.rs:2034`, `main.rs:608`) |
| The inverted guard at `tests/integration_output_collision.rs:1029-1041`, comment as quoted | **[CONFIRMED]** verbatim |
| A4 — `grep -c report src/specialty.rs` → 0 | **[CONFIRMED]** |
| The repro: `--paired A/reads.fq B/reads.fq` exits 0 with reports beside their own mates | **[CONFIRMED]**, exactly the layout in the plan |
| A3 — the report paths really are in the candidate list, so the new collision is caught | **[CONFIRMED]** by the `-o` proxy: refused, `out/` empty, message `the outputs named from A/reads.fq and B/reads.fq would be written to the same file, out/reads.fq_trimming_report.txt` — #397 provenance names both inputs as the plan promises |
| Multi-pair anchors per pair | **[CONFIRMED]** `main.rs:864` `chunks(2)` + `chunk[0]` as R1 into namers that use `input_r1.parent()`; a live 2-pair run put pair 2's primaries in pair 2's R1 directory |
| `--passthrough`'s carrier is the third directory | **[CONFIRMED]** live: `Passthrough: C/i1.fq → C/i1_passthrough.fq` while primaries went to `A/` |

Arm coverage is, as far as I can enumerate it, **complete**. I walked every paired
dispatch branch: trim FASTQ→FASTQ (`main.rs:861`), trim BAM(1 file)→FASTQ
(`run_paired_ubam_single_file`, `main.rs:1795` — reports built from one stem +
one `dir`, `main.rs:1946-1955`, already unified), trim FASTQ→uBAM (`main.rs:2004`),
trim BAM(1)→uBAM (`run_ubam_output_paired_single_file`, `main.rs:2396`), clump
FASTQ→FASTQ (`main.rs:602`), clump FASTQ→BAM Shape A (`main.rs:717`), clump
BAM(1)→BAM Shape B (`main.rs:681`), `--clock`/`--implicon`, `--hardtrim5/3`.
`--demux` is `--paired`-incompatible (`src/cli.rs:1865` comment; barcode outputs are
handled on the SE path in `planned_secondary_outputs`, `main.rs:103-112`), so it is
correctly absent. `--retain_unpaired` is already R1-anchored (`io.rs:368-370`) and is
already a candidate (`main.rs:878-888`).

### 1.2 CRITICAL — a **second** inverted guard exists and the plan does not mention it

`tests/integration_output_collision.rs:1206-1240`, `clump_report_candidates_do_not_over_reject`:

```rust
/// Over-rejection guard: same filename in two dirs WITHOUT -o is legal — both
/// primaries land in R1's directory, each report beside its own mate. Pins that
/// the candidates honour output_dir = None; a builder that resolved reports
/// into one directory would over-reject this.
```

That is the same sentence, almost word for word, as the guard the plan built its whole
"read this before implementing" section around — but on the `--clump_only --paired`
arm. It asserts `ok` for `--clump_only --paired A/reads.fq B/reads.fq` and then
`assert_dir_holds_only(&dir.join("B"), &["reads.fq", "reads.fq_clumping_report.txt"])`.

**[CONFIRMED]** it passes today (`cargo test --release --test integration_output_collision -- --exact … → 3 passed`) and **[CONFIRMED]** that invocation exits 0 today with the split layout. Under the plan it must refuse, so this test goes red.

Concrete failure this causes: the plan's own §Context argues the single most likely way
to get this change wrong is an implementer "fixing" a red over-rejection guard by
narrowing the candidate list. The plan inoculates against that for the trim twin and
leaves the clump twin as an unannounced red test — the exact trap, one arm over. It
also breaks validation 10, which claims to cover "the inverted guard" in the singular.

**Fix:** name both guards in step 9 and validation 10; the clump twin's blocks 1-2
(there are none — it is a single-block test) mean the whole test body is rewritten to
assert the refusal, mirroring `clump_paired_rejects_shared_report_name` (`:1050`) but
without `-o`.

### 1.3 CRITICAL — the change **removes** a refusal, and a test asserts that refusal

`tests/integration_output_collision.rs:1247-1275`,
`clump_paired_rejects_shared_mate_report_without_output_dir`:

```
--clump_only --paired p/a.fq d/x.fq q/b.fq d/x.fq
```

**[CONFIRMED]** today this is refused:
`two outputs named from d/x.fq would be written to the same file, d/x.fq_clumping_report.txt`
(the 2c branch, `io.rs:160-170`), and the test asserts `!ok` plus three empty-of-output
directory listings.

**[DERIVED]** after the change, pair 1 anchors on `p/` and pair 2 on `q/`, so the four
report paths become `p/a.fq_…`, `p/x.fq_…`, `q/b.fq_…`, `q/x.fq_…` and the four
primaries `p/a_clumped_1.fq`, `p/x_clumped_2.fq`, `q/b_clumped_1.fq`, `q/x_clumped_2.fq`
— **eight distinct paths, no collision**. The run succeeds and the test fails on
`assert!(!ok)`.

The trim-mode twin of the same shape (`--paired p/a.fq d/x.fq q/b.fq d/x.fq`) is
**[CONFIRMED]** refused today with the identical message and has **no test**, so it will
silently flip from refuse to succeed.

Two things follow, and neither is in the plan:

1. The plan's Behavior section is **one-directional** — it describes exactly one new
   refusal and says nothing about refusals disappearing. §Behavior 7 and the edge-case
   table therefore under-describe the change. The maintainer approved "the new
   collision is accepted as a loud refusal"; they were not shown that an existing
   deliberate refusal (added by #391, guarded by a test with a written rationale)
   goes away.
2. The new behaviour is, on the facts, *fine* — the shared mate produces two distinct
   reordered outputs in two distinct directories and nothing is overwritten. So this is
   a "record the decision and update the test", not a defect. But it must be recorded:
   a red test asserting a refusal that the plan never predicted is precisely the
   situation where an implementer reaches for the wrong fix.

### 1.4 CRITICAL — a new output-vs-**input** hazard, which is #409's data-loss shape

The plan discusses only report-vs-report collisions. Moving a report into R1's parent
also makes it possible for a report path to land **on one of the inputs**.

**[CONFIRMED]** on `38deb1ee`:

```
$ trim_galore --paired A/reads.fq_trimming_report.txt B/reads.fq     # exit 0
A/  reads.fq_trimming_report_val_1.fq  reads_val_2.fq
    reads.fq_trimming_report.txt                 <- R1's input, intact (@RPT_ reads)
    reads.fq_trimming_report.txt_trimming_report.{txt,json}
B/  reads.fq  reads.fq_trimming_report.{txt,json}
```

R2's report is `B/reads.fq_trimming_report.txt` today — harmless. After the change it
becomes `A/reads.fq_trimming_report.txt`, **which is R1's input path**. If the candidate
builder at `main.rs:906-912` moves with the writer, the pre-flight's output-vs-input
branch (`io.rs:147-156`) refuses — correct. If it does not, the writer destroys a named
input after the run has read it. That is #409 verbatim
(`se_trim_ubam_rejects_output_that_aliases_a_report_input`, `tests/…:1411`), on a new
arm, created by this change.

No validation in the plan targets this direction. Validation 8 ("compare every file
created against that arm's candidate list") would only catch it if the fixture happened
to contain a report-named input, and it does not say to use one. Given the family
history this needs its own test per paired report arm, in both the refusing and the
`--no_report_file` acceptance direction (the shape `tests/…:1411` and `:1440` already
model).

### 1.5 IMPORTANT — `OutputSource` attribution: #397 explicitly left this to #398

`main.rs:872-875`, written by #397:

```rust
// Uniform Pair for paired primaries (#397 decision C2): _val_1 takes its
// stem from R1 but its directory from R1 too, and _val_2 mixes both — naming
// one mate would encode a claim about which supplies the directory, which is
// exactly what #398 may change.
```

The previous change left a note saying #398 has to decide this, and the plan does not
mention it. After the change, both report candidates (`main.rs:908-910`) and the
passthrough candidate (`main.rs:901`) take their **filename** from one input and their
**directory** from R1 — the precise condition `OutputSource::Pair` was introduced for
(`io.rs:94-96`).

My recommendation is to **keep `OutputSource::Input` and say why in the plan**, because
the message quality argues for it and switching would make it worse:

- Keeping `Input`: two distinct sources, one path → the 2a branch (`io.rs:171-180`),
  which prints *"the outputs named from A/reads.fq and B/reads.fq would be written to
  the same file, …"* **[CONFIRMED]** live via the `-o` proxy. Followable: rename one input.
- Switching to `Pair(A,B)`: both candidates get the same source identity
  (`OutputSource::identity`, `io.rs:103-111`) → the **2c** branch, which prints *"List
  each input once — a file that appears in more than one pair is reported once per
  pair"* **[CONFIRMED]** live (that is what the shared-mate clump case prints). That
  advice is wrong here: there are two inputs and renaming one does fix it.

So the finding is not "the plan chose wrong" — it is that the plan is silent on a
decision a previous change explicitly deferred to it, and the naive "improvement"
degrades a user-facing message that shipped last week.

### 1.6 IMPORTANT — the newly-refused population is larger than the plan states

Behavior §7, the edge-case table, the changelog instruction (step 11) and the hint
reword (step 8) all describe the new collision as *"two mates share a filename"* —
i.e. within one pair. The actual rule is broader: **any two report-keyed inputs whose
filenames fold-equal and whose pairs anchor on the same directory.**

**[CONFIRMED]** the cross-pair case is real and newly refused:

```
$ trim_galore --paired A/s.fq B/t.fq A/u.fq B/s.fq      # exit 0 today
A/ s_val_1.fq s_val_2.fq t_val_2.fq u_val_1.fq s.fq_trimming_report.{txt,json} u.fq_…
B/ s.fq t.fq  s.fq_trimming_report.{txt,json}  t.fq_trimming_report.{txt,json}
```

Pair 1's R1 and pair 2's R2 share the filename `s.fq`; both pairs anchor on `A/`.
After the change both reports resolve to `A/s.fq_trimming_report.txt`. The primaries
are `A/s_val_1.fq` and `A/s_val_2.fq` — distinct — so **the reports are the only
collision** and this is a genuinely new refusal in a shape the plan does not describe.

I checked the obvious variant where the two *R2s* share a filename
(`--paired A/x.fq B/y.fq A/z.fq C/y.fq`) and it is **[CONFIRMED]** already refused
today, on the primary `A/y_val_2.fq` — so that one is not new. The R1-vs-R2 crossing
above is the one that matters.

Consequences for the plan: the hint reword must not say "mates", the changelog must
not either, and the edge-case table should carry this row.

Related, on the hint: `PAIRED_REPORT_HINT` is attached to the whole
`preflight_output_collisions` call, so it is printed for **primary** collisions on those
arms too — **[CONFIRMED]**, I saw it appended to an `A/y_val_2.fq` refusal. So the
reword has to read correctly for a primary collision as well, where "a shared input
directory" is equally not the cause. Something on the shape of *"every output of a pair
is written to R1's directory (or `--output_dir`) and named from the input filename
alone, so two inputs sharing a filename collide there"* covers both; the plan's
instruction ("reword so the stated cause matches the new rule") does not flag that the
constant is shared with primary-collision messages.

### 1.7 IMPORTANT — A6's stated mechanism is wrong, and the wrong reading is a trap

A6: *"`Path::parent()` on a bare filename yields `Some("")`, and the existing
`unwrap_or(Path::new("."))` already handles it."*

The first half is right and the second is not: because `parent()` returns `Some("")`,
`unwrap_or(".")` **does not fire** — the empty path is what carries the behaviour, and
`Path::new("").join("x")` is `"x"`, not `"./x"`. The `unwrap_or(".")` is effectively
unreachable for any real input path.

The conclusion (extraction is safe) still holds, but the muddled mechanism invites a
concrete regression: if `pair_output_dir` normalises the empty case to `"."` — which
reads like a tidy-up and which A6's wording actively suggests is already happening —
then **every** paired path printed to stderr gains a `./` prefix, across all five
primary namers plus the reports and the carrier. Today **[CONFIRMED]**:

```
  Output R1: r1_val_1.fq
Report: r1.fq_trimming_report.txt
```

Collision detection would survive it (`lexical_normalise` absolutises and drops
`Component::CurDir`, `io.rs:46-63`), so nothing would fail loudly — it is a silent
cosmetic regression on every paired run. Restate A6 to say what actually happens, and
add the assertion to validation 5.

### 1.8 IMPORTANT — step 6's rationale for `--passthrough` is factually wrong

Step 6: *"Prefer adding an explicit `output_dir` argument at the call site over changing
the namer, so the single-pair SE-shaped uses stay put."*

Two errors. `passthrough_output_name` **already** takes `output_dir` as its second
parameter (`io.rs:387-391`) — there is no argument to add; the minimal change is to pass
`Some(&pair_dir)`. And there are no SE-shaped uses to protect: **[CONFIRMED]** by grep,
the only non-test call sites are `main.rs:895` (candidate) and `main.rs:1499` (writer),
both on the paired path, and `--passthrough` requires `--paired` with exactly two inputs.

Worse, the minimal change leaves debris the plan does not account for:

- The namer's `output_dir == None` branch becomes **unreachable from production**, while
  its docstring (`io.rs:375-386`) still describes deriving from the carrier's own parent.
- Four of its six unit tests pass `None` and assert the carrier lands beside the carrier
  input — `io.rs:853-858`, `:860-865`, `:874-879`, `:888-900`. They keep passing while
  pinning a path production can no longer produce, i.e. they become **vacuous**. This
  repo has a written norm against exactly that:
  `tests/integration_output_collision.rs:803-804` — *"if the two ever drift the guard
  would silently protect paths the run never touches — invisible to every rejection
  test."*

See §5 for the alternative I would recommend instead.

### 1.9 Blocks 1-2 of the trim guard are genuinely unaffected

Both use inputs in one directory passed as bare relative names with `run_in`'s cwd set
there (`tests/…:1007-1014`, `:1018-1027`) **[CONFIRMED]** by reading the helper at
`:54-64`. R1's parent is `""`, which is already every input's parent, so no path moves.
The plan's instruction (rewrite block 3, keep 1-2) is correct — subject to §1.7: if
`pair_output_dir` returns `"."`, blocks 1-2 still pass (they assert `is_file()` on
absolute joins), so they will *not* catch that regression.

### 1.10 Smaller logic points

- **The clump-paired-BAM row of the writer/candidate table is a no-op, not a change.**
  `clump_only.rs:1082-1084` keys the single per-pair report on `inputs[0]` — R1 — so its
  directory is *already* R1's parent, and `main.rs:740-744` already passes
  `from_ref(&chunk[0])` **[CONFIRMED]**. Listing it as an arm to change is harmless but
  costs the implementer a puzzled detour; say "already unified, keep it that way".
- **Validation 8 says "all four arms"; the table lists five.** With the row above being
  a no-op the true count of changing arms is four, so the number happens to work out —
  but the plan should not leave the reader to reconcile it.
- **`clump_report_candidates` has five call sites, not two.** `main.rs:616` (paired
  FASTQ), `:653` (SE FASTQ), `:698` (paired Shape B, one input), `:740` (paired BAM
  Shape A), `:798` (SE BAM) **[CONFIRMED]**. Step 7 says "the paired callers pass the
  pair directory and the SE caller passes `output_dir`" — three of the five are neither
  cleanly "the SE caller" nor need a pair directory. Adding a required parameter makes
  this a compile error rather than a silent miss, so the risk is low, but the sentence
  should enumerate.
- **The candidate site for clump-paired-FASTQ is a closure, not the helper.** The table
  cites `main.rs:134` (inside `clump_report_candidates`' body); the site that has to
  change is the closure at `main.rs:609-622`, invoked per pair by `run_specialty_paired`
  (`main.rs:2637-2653`). Cite that.
- **`--fastqc` needs nothing, and the plan should say so.** **[CONFIRMED]** live:
  `--paired --fastqc A/reads.fq B/reads.fq` puts `reads_val_1_fastqc.{html,zip}` in `A/`
  beside the trimmed output — `fastqc::run` (`src/fastqc.rs:37-59`) is handed the
  **output** path, so it follows the primaries, and the carrier's FastQC will follow the
  carrier automatically once the carrier moves. Add it to the out-of-scope list so the
  next reader does not have to re-derive it. (Separately, FastQC outputs are in no
  candidate list at all — pre-existing, `#414`'s business, not this plan's.)
- **`--basename` cannot mask the new collision, and cannot create one.** Reports ignore
  `--basename` (`io.rs:551-564`), and `--basename` is rejected for multi-pair, so no
  cross-pair primary interaction is reachable. The plan's edge-case row is correct.

---

## 2. Assumptions

| | Verdict |
|---|---|
| **A1** five primaries share the expression → pure extraction | **[CONFIRMED]** by reading all five. Sound. |
| **A2** report namers separate filename from directory | **[CONFIRMED]** `io.rs:535-580`. Sound. |
| **A3** report paths already in the candidate list for all four paired report arms | **[CONFIRMED]** empirically for the trim FASTQ arm via the `-o` proxy (refusal + provenance + empty `out/`); **[CONFIRMED]** by reading for the other three (`main.rs:2023-2029`, `main.rs:616-620`, `main.rs:740-744`). The accepted-refusal decision rests on solid ground. |
| **A4** `--clock`/`--implicon` write no reports | **[CONFIRMED]** `grep -c report src/specialty.rs` → 0. |
| **A5** SE paths never call the paired namers | **[CONFIRMED]** by grep: the paired namers' only non-test callers are the paired dispatch branches. |
| **A6** bare-filename `parent()` | mechanism **wrong**, conclusion right — see §1.7. |

### Unstated assumptions the plan relies on

1. **That the change only ever adds refusals.** False — §1.3.
2. **That report-vs-report is the only new collision class.** False — output-vs-input,
   §1.4.
3. **That `OutputSource` attribution stays correct.** Unexamined, and #397 left a note
   asking — §1.5.
4. **That `pair_output_dir` will reproduce the empty-path case exactly.** Load-bearing
   for stderr output on every paired run — §1.7.
5. **That the docs paths in step 10 exist.** They do not as written: the tree has
   `Docs/src/content/docs/guide/outputs.md` (paired table at lines 19-25),
   `…/guide/paired-end.md`, `…/modes/clump-only.md` **[CONFIRMED]** by `find`. Also
   **`Docs/src/content/docs/modes/passthrough.md` is missing from the plan's docs list**
   even though decision 2 moves the carrier — its output table at `:25` and the
   `--basename` sentence at `:28` describe the three outputs with no directory
   statement. (Its worked example at `:85-99` passes `--output_dir /tmp`, so the example
   itself is safe **[CONFIRMED]**.)
6. **That "File naming matches v0.6.x exactly, so existing pipelines continue to work
   without changes"** (`Docs/…/guide/outputs.md:6`) survives. It is a filename claim, not
   a directory claim, so it is not falsified — but it sits four lines above the table the
   plan is editing and will read as a compatibility promise about the thing that just
   changed. Worth a qualifying clause.

### A fact that strengthens decision 1, which the plan does not have

**Perl 0.6.11 wrote every output into one directory — the CWD — never beside the input.**
**[CONFIRMED]** from the tagged source: `git show 0.6.11:trim_galore` line 608,
`my $output_filename = (split (/\//,$filename))[-1];` strips the input's directory, and
line 620 opens `$output_dir.$report`, where `$output_dir` is `''` unless `-o` was given
(lines 3258-3272). Primaries go through the same `$output_dir.$out1` at lines 152-175.

So Perl had **no** asymmetry to inherit: it was already one directory per run. "Beside
the input" is a v2.x invention, and #398 moves back toward Perl's shape while choosing a
different anchor. Three implications:

- Decision 1 is *more* defensible than the plan argues — it restores a property the
  original had.
- The changelog should not claim or imply v0.6.x directory parity, because R1's parent
  is not what Perl did either.
- It is worth the maintainer knowing that the CWD anchor — which `--hardtrim5/3`,
  `--clock` and `--implicon` already use (`Docs/…/modes/hardtrim.md:42`,
  `modes/clock.md:49`, `modes/implicon.md:32`: *"never beside the input"*) — is the
  Perl-faithful one. **I am not reopening decision 1**; R1's parent is a reasonable
  choice and it is the maintainer's call. But the plan's open question ("should
  `pair_output_dir` also back the `--hardtrim5/3` CWD rule?") is answered on richer facts
  than the plan has: after this change the codebase has three directory rules (CWD for
  specialty, R1's parent for paired, the input's parent for SE) where Perl had one.
  That is worth one honest sentence in the plan rather than a "recommend no".

---

## 3. Efficiency

Nothing to report. One `PathBuf` per pair at path-derivation time, and if callers compute
`pair_dir` once and pass `Some(&pair_dir)` the change is net *fewer* allocations than
today's per-namer rebuild. No new I/O, no new syscalls, nothing inside a per-record loop.
The plan's "Nil" is right.

`pair_output_dir` returning `PathBuf` by value is the right shape — a `Cow<Path>` to dodge
one allocation per pair would be noise at this call frequency.

---

## 4. Validation sufficiency

Validations 1-11 are well-constructed for the failure mode the plan is worried about
(writer/candidate drift on the arms it lists), and validation 3's whole-directory-listing
assertion and validation 9's inverted control are both the right instincts. The gaps are
all in classes the plan did not know about.

**Answer to validation 7, done rather than restated.** The Perl matrix is **unaffected**,
for two independent reasons **[CONFIRMED]** by reading `.github/workflows/ci.yml`:

- Every output-producing invocation passes `-o`. The only `trim_galore` invocations
  without `-o` are `:466` (odd-count rejection), `:481` (single-FASTQ-under-`--paired`
  rejection), `:490` (R1==R2 rejection), `:504` (`--basename` multi-pair rejection) — all
  four assert non-zero exit and write nothing — plus `:551`/`:562` (`--clock`, a
  CWD-output specialty mode with no reports) and `:687` (SE alias test).
- Independently, every fixture lives in `test_files/`, one shared directory, so R1's
  parent *is* R2's parent for every paired invocation in the file.

The Perl-parity steps themselves (`:383-398`, `:411-417`, `:764-789`) all pass `-o` to
both binaries. Nothing to change, and no fixture needs editing.

### Gaps

| # | Gap | Why it matters |
|---|---|---|
| **G1** | No validation for output-vs-**input** (a report landing on the other mate's input) | §1.4 — reachable **[CONFIRMED]**, and drift there is silent data loss, not a refusal. Needs a test per paired report arm plus the `--no_report_file` acceptance sibling. |
| **G2** | No validation for the newly-refused **cross-pair** shape | §1.6 — `--paired A/s.fq B/t.fq A/u.fq B/s.fq`, where reports are the only collision. |
| **G3** | No validation for the **un-refused** shared-mate shape | §1.3 — an existing test asserts the opposite; both the clump and the (untested) trim variant need a decision and coverage. |
| **G4** | Validation 10 names one inverted guard; there are two | §1.2. |
| **G5** | Nothing asserts stderr paths keep their exact spelling | §1.7 — a `"."`-returning helper is a silent `./`-prefix regression on every paired run. Cheapest home: extend the existing bare-name assertions. |
| **G6** | Nothing exercises the reworded hint against a **primary** collision | §1.6 — the constant is shared, so a reword tuned to reports can end up nonsense on a `_val_2` collision. |
| **G7** | Validation 8's fixtures are unspecified | It checks "every file created == the candidate list", which is the right property, but with same-directory fixtures it is satisfied vacuously — the paths do not move. Each arm's fixture must have R2 in a *different* directory from R1, and one variant must include a report-named input (G1). |

One trap the plan gets right and should keep: validation 11's instruction to count the
delta against the 578 baseline rather than eyeballing "ok". Note the corollary for
targeted runs — `cargo test --test integration_output_collision -- --exact <names>` and
check the passed count; I got `3 passed; 45 filtered out`, which is how you know the
filter bit.

---

## 5. Alternatives for the implementation shape

The three behaviour decisions are fixed; these are about how the code is arranged.

**5.1 Make the writer and its candidate builder call one shared derivation, per arm
(recommended).** The plan's central stated risk is that eight independent call sites (four
writers, four candidate builders) must each be given the same directory. Its mitigation is
discipline — "step 7 is deliberately worded to land with each writer". But the #383 /
#388 / #391 / #409 family is *exactly* what discipline-based coupling produces. A helper
per arm that both sides call removes the class instead of guarding it:

```rust
// io.rs
pub fn paired_report_paths(r1: &Path, r2: &Path, output_dir: Option<&Path>) -> [PathBuf; 4]
```

or, better shaped to the existing types, a function returning the two `PairedReportFile`s
that `main.rs:1729-1738` and `main.rs:2351-2360` both build by hand today, with the
candidate builders consuming the same function's paths. The two trim arms are already
near-duplicates of each other (they were pulled into the shared `write_paired_reports`
for precisely this reason, per the comment at `main.rs:1716-1717`) — the *path*
derivation is the half that was left un-shared. This is a bigger diff than the plan's,
and it is the diff that makes instance #6 structurally impossible on these arms rather
than merely tested for.

**5.2 If the plan's shape is kept, put the invariant where it can be read.** The plan
ends with "no signature change to the three shared report functions", which means nothing
at `report_name`'s definition (`io.rs:550-564`) tells the next person that paired callers
must supply the pair directory. Per the repo's own comment policy — *"put an invariant on
the thing it constrains, not in the function that relies on it"* — a one-line doc note on
the three report namers is worth more than the paragraph in `pair_output_dir`'s docstring,
because the failure mode is a *future* paired caller passing `output_dir` straight through.

**5.3 Fold `input_r1` into `passthrough_output_name` rather than passing `pair_dir` at the
call site (recommended over step 6).** Every other paired namer takes `input_r1` and
derives the directory internally; making the carrier's namer match
(`passthrough_output_name(input_r1, input_passthrough, output_dir, basename, gzip)`)
(a) puts the rule in one place instead of two call sites, (b) keeps the docstring honest,
and (c) forces the four `None`-passing unit tests to be updated rather than quietly
becoming vacuous (§1.8). It is a slightly larger diff and a strictly safer one.

**5.4 Consider `pair_output_dir` for the SE namers too, or say why not.** The SE namers'
`match output_dir { Some(d) => d.join(f), None => input.parent()…join(f) }` computes the
same path as the paired `map/unwrap_or_else` form, so one helper could serve all twelve
namers with no behaviour change. The plan's choice to leave SE alone is the lower-risk
one and I would keep it — but the plan should say that the *expression* is shared with SE
and that only the *anchor argument* differs, otherwise the next person reads
`pair_output_dir` as encoding something paired-specific that it does not.

---

## 6. Action items

### Critical — the plan should not go to implementation without these

1. **§1.2** — Add `tests/integration_output_collision.rs:1206-1240`
   (`clump_report_candidates_do_not_over_reject`) to step 9 and validation 10 as the
   second inverted guard, with the same "rewrite, do not narrow" warning. **[CONFIRMED]**
   it passes today and its invocation exits 0 with the split layout.
2. **§1.3** — State in Behavior that the change also **removes** refusals, and decide the
   fate of `clump_paired_rejects_shared_mate_report_without_output_dir`
   (`tests/…:1247-1275`), which asserts a refusal that disappears. Note the untested
   trim-mode twin flips the same way. **[CONFIRMED]** pre-change refusal;
   **[DERIVED]** post-change success.
3. **§1.4** — Add the output-vs-**input** class to Behavior and to Validation: a mate
   named like the other mate's prospective report. **[CONFIRMED]** exits 0 today; after
   the change it must refuse, and a candidate that fails to move destroys an input
   (#409's shape). Cover each paired report arm plus a `--no_report_file` sibling.

### Important

4. **§1.5** — Record the `OutputSource` decision explicitly (recommend: keep `Input`,
   because `Pair` drives the message into the 2c branch whose "list each input once"
   advice is wrong for two-input collisions). `main.rs:872-875` names #398 as the change
   that has to decide this.
5. **§1.6** — Restate the newly-refused population as "any two report-keyed inputs whose
   filenames fold-equal and whose pairs anchor on the same directory", and fix the hint
   reword, the changelog sentence and the edge-case table accordingly. The cross-pair
   shape is **[CONFIRMED]** newly refused. Note that `PAIRED_REPORT_HINT` is also printed
   for primary collisions, so the reword must fit both.
6. **§1.7** — Rewrite A6 to say what actually happens (`parent()` → `Some("")`, so
   `unwrap_or(".")` never fires), and require `pair_output_dir` to return the empty path
   rather than `"."`. Add the stderr-spelling assertion to validation 5.
7. **§1.8 / §5.3** — Drop step 6's incorrect rationale (`passthrough_output_name` already
   has an `output_dir` parameter; there are no SE-shaped uses — **[CONFIRMED]** by grep)
   and prefer folding `input_r1` into the namer, which also keeps `io.rs:853-900`'s four
   `None`-passing unit tests meaningful instead of vacuous.
8. **§2.5** — Fix the docs paths (`Docs/src/content/docs/…`) and add
   `modes/passthrough.md` to step 10, since decision 2 moves the carrier. Consider a
   qualifying clause on `guide/outputs.md:6`'s "matches v0.6.x exactly".
9. **§5.1** — Consider the shared writer+candidate derivation. The plan's own diagnosis of
   the risk argues for removing the coupling rather than documenting it.

### Optional

10. **§1.10** — Mark the clump-paired-BAM row as already unified (a no-op) rather than an
    arm to change; reconcile validation 8's "four arms" with the five-row table;
    enumerate `clump_report_candidates`' five call sites; cite `main.rs:609-622` (the
    closure) rather than `main.rs:134` for the clump-paired-FASTQ candidates.
11. **§1.10** — Add `--fastqc` to the out-of-scope list with the reason (**[CONFIRMED]**
    anchored on the output file, so it follows the primaries and the carrier for free).
12. **§2 (Perl)** — Add the v0.6.11 fact: Perl stripped the input's directory and wrote
    everything to `$output_dir`, defaulting to CWD **[CONFIRMED]** at
    `git show 0.6.11:trim_galore` lines 608/620/3258-3272. It strengthens decision 1,
    changes what the changelog may claim, and gives the "should `--hardtrim` share this
    helper?" open question a real answer instead of a preference.
13. **§4** — Record validation 7's answer in the plan (CI is unaffected: every
    output-producing invocation passes `-o`, and every fixture shares `test_files/`), so
    nobody re-runs it.
14. Note that `target/release/trim_galore` was stale when this review began; validation
    9's control build must come from HEAD.

---

## 7. Verdict

**REVISE.**

The plan is unusually good on the thing it set out to be good at: the extraction is
genuinely behaviour-preserving (I checked all five namers), the writer/candidate table is
accurate as far as it goes, A3 — the assumption the whole accepted-refusal decision rests
on — holds empirically, and the instinct to treat the inverted guard as a decision record
rather than a broken test is right.

What it misses is that the change's behavioural footprint is wider than "one new refusal
in one place". There is a **second** inverted guard on the clump arm; an existing
deliberate refusal **disappears**; and moving reports into R1's parent opens an
output-vs-**input** path that is the exact shape of the data-loss bug fixed three commits
ago. Each of the three is a red or newly-green test that the plan does not predict, in a
change whose own §Context argues that unpredicted red tests are how this bug family
recurs. They are all fixable on paper — none of them threatens the three decisions — but
the plan should carry them before an implementer meets them.
