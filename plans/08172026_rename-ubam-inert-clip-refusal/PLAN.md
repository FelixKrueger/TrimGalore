# Plan — narrow the #408 `--rename` + uBAM refusal so an inert Read 2 clip flag no longer triggers it

**Status:** v2, revised after dual plan review. Awaiting an implementation trigger. Not implemented.
**Issue:** [#450](https://github.com/FelixKrueger/TrimGalore/issues/450).
**Base:** `origin/dev` @ `4f7b907` (83 ahead of `master`).
**Origin:** M-3, opened for decision in `plans/08172026_clip-r2-se-warning/PLAN.md` §13. Felix chose
**narrow the guard** on 2026-08-17.

## Revision history

**v2 (2026-08-17)** — after `PLAN_REVIEW_reviewer-A_g1.md` and `PLAN_REVIEW_reviewer-B_g1.md`. Both
reviewers verified all eight load-bearing claims as TRUE and found no defect in the proposed behaviour.
Three things changed materially:

1. **V4's mutation pointed the wrong way** (A's C1). It weakened the guard and expected tests that assert
   *success* to fail — they would have stayed green, and the honest reading of that green run is "the
   narrowing is not load-bearing", the opposite of the truth. V4 and V5 were also the same mutation.
   Both rewritten. B had explicitly blessed V4-V6 as well-chosen; A carried the mutation through to the
   specific assertion and found the sign error. **This is the finding the two reviewers disagreed on.**
2. **The docs check was a false-negative grep** (A's C2 = B's C1). The pattern was the *error-message*
   string; the docs paraphrase. Three surfaces state the gate and the narrowing falsifies all three.
   Step 13 is now three named edits.
3. **Scope widened to the full fold-in** (A's I1 = B's I1, and Felix's call). Read ownership sits at six
   sites, not three. The decisive argument is A's risk ordering: as scoped in v1 the plan rewrote the one
   display site that **no test pins** and left alone the two that are pinned in both directions.

Also from review: the truth table gained two rows (I4 both), V7 gained a valid-cell set and a
non-degeneracy assertion (A's I5), `any_effective`'s docstring was narrowed (B's I2), `user_clips()` is
doc-scoped (B's I5), and eight citations were corrected. A's X1 — a genuine **under**-refusal on the
`--rrbs` path — is filed as [#459](https://github.com/FelixKrueger/TrimGalore/issues/459) and stays out
of scope.

---

## Goal

A single-end run of

```
trim_galore --rename --output-format ubam --clip_R2 5 se.fastq
```

currently warns that `--clip_R2` is not used in this mode and is being ignored, then **refuses the run
because `--clip_R2` is set**. It announces that it is ignoring the flag, then refuses on account of it.
(Both messages are single unwrapped lines; `cli.rs:1168-1170` emits one line.)

After this change the run **proceeds**, keeping the warning. Every other refusal the #408 guard performs
is unchanged, and each is pinned by a test a mutation can turn red.

Secondary goal, and the reason this is not a one-line edit: **read ownership of a clip flag is encoded at
six places**, and this defect is what that duplication produced. One `ClipFlag` enum becomes the single
source, and the predicate `paired || !read_2` — currently spelled four times — becomes one method.

---

## Context

> **Line-number provenance.** Every number below is against **`origin/dev` @ `4f7b907`**, read via
> `git show origin/dev:<file>` into a temp copy. **Do not read the working tree.** A parallel session
> has committed repeatedly on top of this base and renamed its branch; at the time of writing `HEAD` was
> `992235e` on `ci/clumpify-memory-guard-439`. Files that differ from `origin/dev`: `src/main.rs`
> (`+2` before line 220, so the guard reads `:398` there and `:396` here), `src/parallel.rs` (`+13`, so
> the R2 trim call reads `:496` there and **`:495`** here), `src/clump.rs`, `src/clump_only.rs`,
> `src/fastqc.rs`, `tests/integration_preflight_tripwire.rs`, `examples/mk_clump_fixture.rs`,
> `CHANGELOG.md`, `.github/workflows/ci.yml`, `.github/clumpify-memory-baseline.tsv`.
> `src/library.rs`, `src/cli.rs`, `src/trimmer.rs`, `src/report.rs`, `src/fastq.rs` and `tests/` agree in
> both. **`CHANGELOG.md` is in that list and step 12 edits it, so expect a conflict once `dev` advances.**

### The guard as shipped — `src/main.rs:389-401`

```rust
// #408 — format-gated, so it cannot live in `Cli::validate()`; sited after the
// two structural pair checks above so those report their own defect first. The
// clip-flag list must match every `--rename`-driven `append_to_id` site
// (`trimmer.rs` clip_5/clip_3, `specialty.rs` hardtrim) — without one of them
// set, `--rename` appends nothing and there is nothing to lose.
// The clip list reads through `effective_clips()` so a `--library` preset
// (#440) counts here exactly as the four flags it stands in for would.
if cli.rename
    && matches!(cli.output_format, trim_galore::cli::OutputFormat::UBam)
    && (cli.effective_clips().any_set() || cli.hardtrim5.is_some() || cli.hardtrim3.is_some())
    && input_formats.iter().any(|f| !matches!(f, InputFormat::UnalignedBam))
```

The guard sits at 4-space indent directly in `main()` — the preceding `if cli.paired` block opens at
`:372` and closes at `:387` — so a single-end run reaches it. `ResolvedClips::any_set()`
(`src/library.rs:141`) ORs all four values including the Read 2 pair, so an inert Read 2 flag satisfies
the clip term. `cli.validate()` runs first at `:242`, which is why the warning prints before the bail.

**The guard's own comment is the specification its condition violates.** It states the criterion as "the
clip-flag list must match every `--rename`-driven `append_to_id` site … without one of them set,
`--rename` appends nothing and there is nothing to lose." On a single-end run with only Read 2 flags set,
no `append_to_id` site is reachable.

### Why the Read 2 values are categorically unreachable on a single-end run

`src/trimmer.rs:246-255` selects the clip value by `is_r2`, and `append_to_id(":clip5:…")` /
`(":clip3:…")` fire only inside the two `if let Some(n)` arms that follow (`:264`, `:273`).

There are **nine** production `trim_read` call sites. Three pass `is_r2 = true`:

| Site | Enclosing function |
|---|---|
| `trimmer.rs:480` | `run_paired_end` (`:381`) |
| `trimmer.rs:731` | `run_paired_end_to_bam` (`:704`) |
| `parallel.rs:495` | `classify_paired` (`:479`) |

The other six pass `false`: `trimmer.rs:319`, `:479`, `:645`, `:730`, `parallel.rs:494`, `:1111`. Every
entry into the three paired functions is gated on `cli.paired` (`main.rs:862`, `:893`, `:2046`, `:2072`;
`run_paired_end_parallel` is reached only from `main.rs:1606` and `:1885`, both inside those). So
`!cli.paired` ⟹ no `is_r2 = true` call ⟹ `config.clip_r2` and `config.three_prime_clip_r2` are never
read for trimming ⟹ nothing is appended. Every other mention of those fields in `src/` is config
construction or report plumbing.

Both reviewers independently confirmed the enumeration of `true` sites is complete, including that
`specialty.rs`'s four `--rename` appends (`:62`, `:105`, `:164`, `:217`) are hardtrim — covered by the
guard's own hardtrim terms — that `specialty.rs:310-311`/`:409-410` are `--clock`/`--implicon` UMI
appends whose modes are refused for uBAM output at `cli.rs:630-640`, and that `demux.rs` contains no
`append_to_id`.

This is a **categorical** absence, and that is the boundary of this change. The guard is deliberately
conservative about the *statistical* case: `clip_5prime` (`src/fastq.rs:110-119`) returns `None` when
`seq.len() <= n`, so even a live Read 1 flag appends nothing for short reads — and the error text says
the whole run is refused rather than decided per record precisely because of that. This change removes
only the case where the site cannot fire for **any** read.

### The repo has already accepted this principle one step back

`tests/integration_ubam_out.rs:1078`, `rename_without_clip_flag_is_accepted_into_ubam`:

> `--rename` alone appends nothing (`append_to_id` is reached only under a clip flag), so refusing it
> would break wrappers that pass the flag unconditionally.

A single-end run whose only clip flag is a Read 2 one is in exactly that class.

### No existing test asserts the behaviour being changed

`--clip_R2` / `--three_prime_clip_R2` appear nowhere in `tests/integration_ubam_out.rs`. The four
refusal tests use `--clip_R1` (`:1006`), `--library emseq` (`:1044`), `--hardtrim5` (`:1172`),
`--hardtrim3` (`:1204`); the two negative tests (`:1232`, `:1263`) use `--clip_R1`. All six behave
identically after the change.

### The six sites encoding read ownership

| # | Site | Value source | Read-2 rule today | Pinned by a test? |
|---|---|---|---|---|
| 1 | `ResolvedClips::flag_summary(paired)` — `library.rs:122-138` | resolved | `if read_2 && !paired { continue }` (`:130`) | **Yes, 3-4 ways**: `library.rs:322-325`, `:331-334`, `integration_library.rs:272`, `report.rs:1399` |
| 2 | The inert-flag warning ladder — `cli.rs:1138-1172` | **raw `self.clip_r*`** (Assumption 3) | R1 entries `None` unless `all_four`; R2 always live | **Order: nowhere.** Presence yes, order no |
| 3 | `ResolvedClips::any_set()` — `library.rs:141`, sole caller `main.rs:398` | resolved | **none — this is the defect** | No |
| 4 | `ClipOverride.flag: &'static str` + `.read_2: bool` — `library.rs:87-94`, populated by four `pick(…)` calls with hard-coded literals at `library.rs:186-199` | — | four hand-written bools | `every_override_knows_which_read_it_belongs_to`, `library.rs:296-317` |
| 5 | `main.rs:1128` — `.filter(\|o\| cli.paired \|\| !o.read_2)` | resolved | open-coded predicate | **Yes, both directions**: `integration_library.rs:339`, `:371` |
| 6 | `report.rs:287` — `.filter(\|o\| config.paired \|\| !o.read_2)` | resolved | the same predicate again | **Yes**, as above |

The predicate `paired || !read_2` is spelled **four** times (sites 1, 5, 6, and the new one). The v1 plan
covered sites 1-3 only, which meant rewriting the single site nothing pins while leaving the two
best-protected ones alone.

---

## Behavior

### B1 — the narrowed guard

1. Prerequisites unchanged: `--rename`, `--output-format ubam`, at least one input detected non-uBAM.
2. The clip term becomes read-aware: a Read 2 clip value counts only on a paired run. The
   `hardtrim5`/`hardtrim3` terms are untouched — hardtrim genuinely appends.
3. Refusal text unchanged. No new message anywhere.

Truth table for FASTQ input + `--rename --output-format ubam`. **All twelve rows become tests** (V3):

| # | Invocation | Today | After | Test |
|---|---|---|---|---|
| 1 | `--clip_R2 5` (single-end) | refused | **runs** + warning | step 9 (new) |
| 2 | `--three_prime_clip_R2 7` (single-end) | refused | **runs** + warning | step 9 (new) |
| 3 | `--clip_R1 3` (single-end) | refused | refused | existing `:1006` |
| 4 | `--three_prime_clip_R1 3` (single-end) | refused | refused | **step 10a (new — untested today)** |
| 5 | `--clip_R2 5 --paired r1 r2` | refused | refused | step 10 (new) |
| 6 | `--three_prime_clip_R2 7 --paired r1 r2` | refused | refused | step 10 (new) |
| 7 | `--clip_R1 3 --clip_R2 5` (single-end) | refused | **refused** | **step 10b (new — see below)** |
| 8 | `--library <all five>` (single-end) | refused | refused | existing `:1044` + step 11 |
| 9 | `--hardtrim5 30` | refused | refused | existing `:1172` |
| 10 | `--hardtrim3 30` | refused | refused | existing `:1204` |
| 11 | `--rename` with no clip flag | runs | runs | existing `:1078` |
| 12 | uBAM input + any clip flag | runs | runs | existing `:1109` |

**Row 7 is the one whose failure mode is the original #408 data-integrity defect** rather than an
over-refusal. A natural mis-reading of the narrowing — "skip the refusal when an inert Read 2 flag is
present" instead of "do not *count* an inert Read 2 flag" — lets it run, and `:clip5:` is then appended
for R1 into a BAM read name. Neither of step 9's tests has an R1 flag and `:1006` has no R2 flag, so
nothing else catches it. Step 11's preset loop covers the predicate incidentally (a preset populates all
four values), but incidentally is not the same as explicitly.

**Row 4 was missing from v1 entirely** and is untested before or after: `git grep -F --
"three_prime_clip_R1" -- tests` never hits `integration_ubam_out.rs`. v1's step 9 cited the #447 lesson
that one representative per branch misses the second entry, and then omitted the second entry of the R1
pair.

### B2 — one shared read-ownership source, all six sites

`ClipFlag` in `src/library.rs` becomes the single place that says which read a clip flag belongs to, and
`ClipFlag::applies(paired)` the single spelling of the predicate. `ClipOverride.flag` becomes a
`ClipFlag` and its `read_2` field is dropped, which makes the four hand-written bools at
`library.rs:186-199` unrepresentable.

**Output-neutral, verified:** the JSON override writer (`report.rs:896-900`) emits only `flag`,
`preset_value`, `user_value` — **`read_2` is not serialised**. `report.rs:898` becomes
`json_escape_string(o.flag.name())` and `report.rs:292` becomes `o.flag.name()`, so `name()` must return
exactly today's literals, leading `--` included, or the JSON and the text report change bytes.

Sites 1, 2, 4, 5 and 6 are **pure refactors**; V7 measures that. Site 3 is the only behaviour change.

Each display site keeps **its own emission order** — see Signature. `ALL` is not a display order.

### B3 — edge cases

| Case | Handling |
|---|---|
| `--library <preset>` on a single-end run | All five presets set `clip_r1 ≠ 0` (`library.rs:69-75`: EmSeq 10, Accel 10, Zymo 10, ScBs 6, Pbat 8), so the Read 1 term keeps the refusal for **every** preset. Step 11 tests all five |
| `--rrbs` directional auto-clip (`clip_r2 = 2`) | Cannot reach the guard: gated on `cli.paired` and set inside `setup_trimming` (`main.rs:1141`, function at `:1051`, earliest call `:880`), after the guard at `:396`. **But it therefore also escapes the guard, which is a separate under-refusal — tracked as [#459](https://github.com/FelixKrueger/TrimGalore/issues/459), deliberately not fixed here.** No change in this commit |
| `--paired` with one interleaved uBAM | `paired == true` so Read 2 values count, and the input-format term excludes the shape anyway. Unchanged twice over |
| `--clip_R2 0` on a single-end run | Runs, as any inert Read 2 value now does |
| `--clip_R1 0` (any run) | **Still refused, and still wrong** — `clip_5prime` guards on `n > 0`. Same defect class, different question. Out of scope, see Q-1 |
| `--clock` / `--implicon` + clip flags | Refused for uBAM output at `cli.rs:630-640`, so they never reach the guard |

---

## Signature

```rust
// src/library.rs

/// One of the four end-clipping flags.
///
/// An enum rather than a string table: the read a flag belongs to is one match
/// arm rather than a fact each consumer re-derives, and a fifth flag would not
/// compile until every consumer handled it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClipFlag { ClipR1, ClipR2, ThreePrimeClipR1, ThreePrimeClipR2 }

impl ClipFlag {
    /// The closed set. Iterated by `any_effective` and by the drift tests; no
    /// display site iterates it, so its order carries no meaning.
    pub const ALL: [ClipFlag; 4] = [ /* R1, 3'R1, R2, 3'R2 */ ];

    /// Flag name as the user would type it. These four strings reach the JSON
    /// report, so they are output: `--clip_R1`, `--clip_R2`,
    /// `--three_prime_clip_R1`, `--three_prime_clip_R2`.
    pub fn name(self) -> &'static str;

    /// A Read 2 flag.
    pub fn read_2(self) -> bool;

    /// This flag names a read the run processes.
    pub fn applies(self, paired: bool) -> bool { paired || !self.read_2() }
}

/// One clip flag given on the command line over a preset that also sets it.
pub struct ClipOverride {
    pub flag: ClipFlag,          // was `&'static str` + a separate `read_2: bool`
    pub preset_value: usize,
    pub user_value: usize,
}

impl ResolvedClips {
    /// The value in force for `flag`, preset folded in.
    pub fn value(&self, flag: ClipFlag) -> Option<usize>;

    /// True when a clip value in force belongs to a read this run processes.
    /// Says nothing about whether the *mode* honours clip flags — the caller's
    /// hardtrim terms cover that.
    pub fn any_effective(&self, paired: bool) -> bool;
}

impl UserClips {
    /// The value the user typed for `flag`, before any preset.
    pub fn value(&self, flag: ClipFlag) -> Option<usize>;
}
```

```rust
// src/cli.rs
impl Cli {
    /// The values as typed. For the inert-flag warning only, which must not
    /// report flags a preset supplied; every consumer of the values *in force*
    /// goes through `effective_clips()`.
    pub(crate) fn user_clips(&self) -> crate::library::UserClips;
}
```

`any_effective`'s doc is deliberately narrow. The wider claim — "can reach a clip site" — is **false on a
shipped path**: under `--hardtrim5 30 --clip_R1 3` nothing reaches a clip site (`specialty.rs` reads no
clip flag) yet the method returns `true`; the refusal is carried by the guard's separate hardtrim term.
Overstating it would reproduce, one level down, the exact defect this plan exists to fix.

**`ALL` is not a display order.** The two display sites emit the four flags differently, and both orders
are user-visible:

| Site | Order today |
|---|---|
| Warning ladder, `cli.rs:1154-1166` | `clip_R1, three_prime_clip_R1, clip_R2, three_prime_clip_R2` |
| `flag_summary`, `library.rs:124-129` | `clip_R1, clip_R2, three_prime_clip_R1, three_prime_clip_R2` |

No single shared order preserves both, so **each display site keeps its own explicit list of variants**
and shares only `name()`, `read_2()`, `applies()` and `value()`. Read ownership is the fact worth
unifying; emission order is presentation, and unifying it would be an unforced user-visible change
inside what is meant to be a pure refactor.

---

## Implementation outline

1. **`src/library.rs` — add `ClipFlag`** with `ALL`, `name()`, `read_2()`, `applies(paired)`. `read_2()`
   is `matches!(self, ClipR2 | ThreePrimeClipR2)` — no string inspection (Assumption 4). Widen the
   module doc comment, which currently scopes the module to `--library` presets; `ClipFlag` serves the
   ladder and the #408 guard, neither of which involves a preset.
2. **`src/library.rs` — `ResolvedClips::value` and `UserClips::value`.** Plain matches on the variant.
3. **`src/library.rs` — replace `any_set()` with `any_effective(paired)`,** iterating `ALL` with
   `applies()`. Delete `any_set`: one caller, no test references. It is a `pub` item in a `pub` module of
   a crates.io-published crate, so removal is technically a public-API change — deliberate, because
   leaving a correctly-named method with the defective behaviour is the worse hazard, and the lib target
   exists for tests rather than downstreams.
4. **`src/library.rs` — `ClipOverride.flag: ClipFlag`, drop `read_2`.** Change `pick`'s first parameter
   to `ClipFlag` and drop its fourth; the four call sites at `:186-199` pass `ClipFlag::ClipR1` etc.
5. **`src/library.rs` — `flag_summary` uses `applies()` and `value()`,** keeping its own order
   (`clip_R1, clip_R2, three_prime_clip_R1, three_prime_clip_R2`). Output must not change; it is pinned
   3-4 ways, so a mistake here fails CI immediately.
6. **`src/cli.rs` — extract `Cli::user_clips()`** with the doc comment above, `pub(crate)`. Note there is
   **no** duplicated `UserClips` literal today (`cli.rs:1228` is the only one); the extraction prevents
   the duplicate step 7 would otherwise create.
7. **`src/cli.rs:1138` — rewrite the warning ladder** over its own order, reading
   `self.user_clips().value(f)`. The condition becomes `value.is_some() && (all_four || f.read_2())`,
   equivalent arm by arm to the shipped nested construction: R1 entries warn only under a specialty mode,
   R2 entries warn whenever `inert` is `Some`. Both reviewers verified the equivalence; V7 confirms it.
8. **`src/main.rs:1128` and `src/report.rs:287` — `.filter(|o| o.flag.applies(cli.paired))`.**
   `main.rs:1131` and `report.rs:292` become `o.flag.name()`.
9. **`src/report.rs:898` — `json_escape_string(o.flag.name())`.** Keep the escape call for uniformity
   even though the four names are known-safe.
10. **`src/main.rs:398` — `cli.effective_clips().any_effective(cli.paired)`.** The comment block at
    `:389-395` is already seven lines against CLAUDE.md's two-line maximum; **collapse it to one line
    stating the criterion** rather than adding an eighth. The mechanism belongs in the commit message.
11. **`src/library.rs` tests — the drift guards.** (a) For every variant,
    `name().ends_with("_R2") == read_2()` — flag name as an independent oracle, legitimate in a test and
    nowhere else. **Replace** `every_override_knows_which_read_it_belongs_to` (`:296-317`) with it rather
    than leaving two near-identical drift tests; the old one cannot survive, since `ClipOverride` no
    longer has a `read_2` field. (b) `ALL` contains all four variants exactly once, via an exhaustive
    match. Note the honest limit: an exhaustive match forces a *visit* to this test when a fifth variant
    appears; nothing forces the variant into `ALL`, whose length is the literal `4`. (c) `name()` returns
    the four exact strings including the leading `--`, since they reach the JSON report and step 11a's
    `ends_with` oracle passes with or without the dashes.
12. **`tests/integration_ubam_out.rs` — rows 1, 2 (new, must run).** `--clip_R2` and
    `--three_prime_clip_R2` on single-end FASTQ, each asserting exit 0, the warning present, the BAM
    written, and no `:clip5:`/`:clip3:` in any read name. **Two different output directories or input
    stems**, since both runs would otherwise write `sp_trimmed.bam`.
13. **`tests/integration_ubam_out.rs` — rows 4, 5, 6, 7 (new, must refuse).** Row 4 is
    `--three_prime_clip_R1` single-end; rows 5 and 6 are the paired mirrors; row 7 is
    `--clip_R1 3 --clip_R2 5` single-end.
14. **`tests/integration_ubam_out.rs` — row 8**: a loop over `emseq`, `accel`, `zymo`, `scbs`, `pbat`.
    `:1044` stays.
15. **`tests/integration_inert_clip_flags.rs` — pin the ladder's emission order** (this file is the
    dedicated home for the ladder and must otherwise stay green **unmodified** by step 7; it also
    contains `preset_alone_on_single_end_does_not_warn` at `:200`, which already test-enforces
    Assumption 3). In `each_inert_flag_gets_its_own_line` (`:219-240`), assert
    `err.find("--clip_R2 5") < err.find("--three_prime_clip_R2 7")`. Add one **all-four-flags** case
    under `--hardtrim5` asserting the four positions ascend — today no test sets all four at once, so
    the ladder's full 4-element order is entirely unexercised.
16. **`docs/src/content/docs/guide/outputs.md:100` and `:145`** — one clause each: the condition is not
    "a clipping flag is set" but a clipping flag that can apply to this run's shape.
17. **`src/cli.rs:312-313`** — the `--help` text omits the clip-flag condition entirely, so it already
    overstates the refusal (`--rename` alone is accepted, `integration_ubam_out.rs:1078`). Felix's call:
    fix it here. There is **no `--help` snapshot test anywhere in `tests/`**, so this text drifts
    silently and permanently — worth one sentence in the commit message.
18. **`CHANGELOG.md`** — one sentence under Unreleased as its own bullet citing #450. Every neighbouring
    bullet is a paragraph and #440's runs 23 lines, so one sentence reads as a deliberate contrast; it is
    correct per the global rule. Expect a conflict, since the parallel session has also edited this file.

---

## Efficiency

Nothing is on a per-record path. `any_effective` runs once per invocation over a 4-element `const` array;
`flag_summary`, the ladder and the two override filters likewise, over at most four items. `ClipFlag` is
`Copy` and one byte, and replacing `ClipOverride`'s `&'static str` + `bool` with it makes the struct
smaller. No allocation is added.

---

## Integration

- **Reads:** `cli.paired`, `cli.rename`, `cli.output_format`, the four clip flags, `cli.library`,
  `input_formats` — all already read at this point in `main()`.
- **Writes:** nothing new. No output file, report field or JSON key changes. The JSON `overrides[].flag`
  strings must stay byte-identical (step 9).
- **Order:** the guard stays after the two structural pair checks; the ladder stays at `cli.rs:1138`,
  still the last fallible position.
- **Downstream:** the JSON report records what was *requested* (#445's invariant) and is untouched.
- **Perl byte-identity:** `Cli::validate()` runs on every invocation, so `Validate vs Perl TrimGalore`
  is load-bearing for step 7.

---

## Assumptions

1. **Fixed rule.** `!cli.paired` ⟹ no `trim_read(…, is_r2 = true)` call. Verified across all nine
   production call sites by me and independently by both reviewers. A future single-end path passing
   `true` would falsify the narrowing; there is no compile-time guard, which is why this is written here
   and named in the commit message.
2. **Fixed rule.** All five presets set a non-zero `clip_r1` (`library.rs:69-75`), so no preset is Read
   2-only. A future Read 2-only preset would correctly begin to proceed on a single-end run.
3. **Fixed rule, and a trap.** The ladder must keep reading **raw** `self.clip_r*`, never
   `effective_clips()`: a preset populates all four `Option`s, so the resolved view would warn about
   flags the user never typed on every single-end `--library` run. Already test-enforced by
   `preset_alone_on_single_end_does_not_warn` (`integration_inert_clip_flags.rs:200`) — stronger than v1
   claimed. The shared table must not quietly re-unify the *value source* while unifying read ownership.
4. **Fixed rule.** Production code must not derive read ownership from a flag *string*. Step 11a uses the
   name as an oracle, which is legitimate there and nowhere else.
5. **Configurable.** The refusal text is unchanged. Making it mention that inert flags do not trigger it
   is a wording decision, not part of the narrowing.
6. **Environment.** Implement in a **separate worktree explicitly based on `origin/dev`** — both
   `isolation: "worktree"` and a bare `git worktree add` base off the default branch (`master`), 83
   commits stale. Verify with `git worktree list` immediately after creating it. The hazard is no longer
   "uncommitted edits in the shared tree" (v1's diagnosis): the parallel session has **committed**, so
   the working tree is clean and reading a working file shows a *newer commit* — which looks like clean
   ground and is not. Disk at time of writing: 103 GiB free, 89% used; a worktree's `target/` costs ~2 GB
   and V7 needs two builds.

---

## Validation

| # | What to verify | How | Expected |
|---|---|---|---|
| V1 | The defect reproduces at base | Build `4f7b907`; `--rename --output-format ubam --clip_R2 5` on a one-record single-end FASTQ | Warning **and** non-zero exit citing `--rename is refused` |
| V2 | The narrowed case runs, identical to omitting the flag | Same invocation after the change, in **two separate output dirs** so the two `sp_trimmed.bam` files cannot collide; compare with `assert_ubam_eq` (`integration_ubam_out.rs:108`, ignores `@PG`) | Exit 0, warning present, records identical, no `:clip5:`/`:clip3:` in any read name |
| V3 | Every row of B1's twelve-row table | One test per row (steps 12-14 plus the five existing) | Each matches the "After" column |
| V4 | **The narrowing is load-bearing** | Mutation: make `any_effective` count Read 2 values regardless of `paired` — i.e. reinstate `any_set`'s behaviour | Rows 1 and 2 (step 12) go **red** |
| V5 | **The paired direction is load-bearing** | Mutation: drop `paired \|\|` from `applies`, so Read 2 values never count | Rows 5 and 6 (step 13) go **red** |
| V6 | Read ownership is load-bearing | Mutation: flip `read_2()` for `ClipR2` only | Step 11a goes red |
| V6b | The R1 arm is load-bearing | Mutation: flip `read_2()` for `ThreePrimeClipR1` only | Row 4 (step 13) goes red |
| V7 | Sites 1, 2, 4, 5, 6 are pure refactors, **emission order and JSON included** | Differential: base-commit binary vs new binary. Capture stderr, `*_trimming_report.txt` and `*_trimming_report.json`; temp files + `cmp`, **not** `diff <(…) <(…)` (blocked in the sandbox). See the cell rules below | Every pair identical except rows 1 and 2 |
| V8 | Nothing else moved | `cargo test --no-fail-fast` (without it the run stops after the lib target and hides ~283 integration results), `cargo fmt --all -- --check`, `cargo clippy --all-targets --release -- -D warnings` | ≥ 782 + new tests passing, 0 failures, clean |
| V9 | Perl parity | The `Validate vs Perl TrimGalore` CI job | Green |

**V4 and V5 point in opposite directions, and that is the whole point.** v1 had V4 weaken the guard and
expect tests that assert *success* to fail. A mutation that loosens a condition can only make more runs
proceed, so those tests would have stayed green — and a green V4 reads as "the narrowing is not
load-bearing", the opposite of the truth. Any future edit to this table must check that each mutation's
direction and its expected failure point the same way.

**V7 cell rules** (A's I5 — a matrix that silently degenerates is worse than no matrix):

- **State the valid-cell set explicitly, never as a product.** `--library` is refused with
  `--hardtrim5`, `--hardtrim3`, `--clock` and `--implicon` (`cli.rs:1115-1135`), so "5 modes × 5 presets"
  contains 20 impossible cells. `--clock`/`--implicon` additionally need the clock fixtures and two files
  *without* `--paired`.
- **Two binaries emitting the same bail compare equal**, so an impossible cell reports PASS while testing
  nothing. Each cell must first assert it produced either a `WARNING: … is not used in this mode` line or
  a `Clipping in force:` line, and only then compare.
- **Include cells with two or more clip flags set**, in both an inert set and a live one: a one-flag cell
  cannot reveal a reordering. Include one all-four-flags cell.
- **Site 1 needs only the preset cells** — `flag_summary` is reached only under `if let Some(preset)`
  (`main.rs:1122`, `report.rs:275-276`). Non-preset cells exercise site 2 alone.
- Neither report nor stderr carries a timestamp, elapsed time or git hash (`config.version` is
  `CARGO_PKG_VERSION`), so byte-comparison is feasible.

**Docs discovery** — record the pattern, since v1's miss came from grepping the error string:
`git grep -n -i -E 'rename.*(refus|reject)|(refus|reject).*rename|clipping flag is set' -- docs README.md src`.

Harness rules for every step: `while IFS= read -r` (zsh does not word-split), `grep -F --` for patterns
starting with `--`, one summary line per item so a truncated output cannot pass as a completed
verification, and check exit status alongside any measured value.

---

## Questions or ambiguities

All resolved. Kept with their answers so a reviewer can see what was decided.

- **Q-1 — `--clip_R1 0` also over-refuses. Leave unfiled, out of scope** (Felix, 2026-08-17).
  `clip_5prime` guards on `n > 0`, so a zero value can never append, yet it triggers the refusal.
  **Not to be quietly folded in** — the diff is one predicate, which is exactly why it needs saying.
- **Q-2 — CHANGELOG: separate bullet citing #450** (Felix), not folded into #440's entry. Step 18.
- **Q-3 — issue filed as #450** (Felix), so the fix PR cites it.
- **Q-4 — refactor scope: the full six-site fold-in** (Felix, after both reviewers rejected the
  three-site middle). The rejected alternatives were the minimal five-line predicate with the
  consolidation filed separately, and keeping three sites with the rest scoped out.
- **Q-5 — X1, the `--rrbs` under-refusal: filed as #459, out of scope** (Felix). B3's row says so.
- **Q-6 — the `--help` overstatement at `cli.rs:312-313`: fix in this commit** (Felix). Step 17.
- **Resolved during planning.** Whether the narrowing weakens a data-safety guard: no. The hazard the
  error text describes is categorically absent in the narrowed case, and the guard's own comment states
  that criterion.

---

## Self-Review

**Efficiency.** Nothing added to a per-record path; `ClipOverride` gets smaller.

**Logic.** Traced the guard's four terms and confirmed only the clip term changes. Confirmed the ladder's
reformulation arm by arm, then declined to trust the reading alone — V7 measures it, and both reviewers
reproduced it independently. Confirmed the guard is not inside an `if cli.paired` branch, which is what
makes the single-end case reachable at all.

**Edge cases.** Adjusted four times. `--clip_R1 0` found while checking `clip_5prime` and held out as
Q-1. The interleaved-uBAM shape found to be excluded twice over. Row 7 (mixed ownership) and row 4
(`--three_prime_clip_R1`) added from review — row 4 is the sharper miss, because v1 cited the very lesson
it then broke. And the `--rrbs` row's verdict changed from "no change" to "under-refuses, tracked as
#459": true about the mechanism, and v1 was silent about the consequence.

**Integration.** Two defects were caught in this plan rather than in the code, and they are worth
separating because they are different failure modes:

1. *Unifying the loop instead of the fact* (caught in v1's own self-review). Making all sites iterate one
   shared `ALL` would have silently reordered user-visible output, since the two display sites emit
   different orders. **Unifying a duplicated fact is safe; unifying the loop that consumes it is not.**
2. *A mutation whose expectation pointed the same way as the mutation* (caught by Reviewer A). V4 could
   not have failed. The general check: if a mutation loosens a condition, it cannot turn a
   should-succeed test red.

**Remaining risks.**
1. **The ordering fix must land at both display sites.** Site 1 is pinned 3-4 ways so a mistake fails CI;
   site 2's order is pinned nowhere today, which is why step 15 adds committed assertions rather than
   relying on V7. That converts the top risk from "measured once before merge" to "held by CI".
2. **Assumption 1 is a whole-codebase claim** verified by enumerating call sites at one commit. No
   compile-time guard exists.
3. **`name()`'s four strings are now output.** They reach the JSON report, so a typo is a silent report
   change rather than a compile error. Step 11c pins them; the `ends_with("_R2")` oracle alone would not.
4. **This is the third consecutive change to this cluster** (#445, #447, this one), and the scope grew
   under review rather than shrank. That is the shape of a design wanting one consolidation — which is
   what Q-4 now chooses to do, so a fourth visit should be read as a signal that something structural
   remains, not as more of the same.

---

## §12 — Implementation notes (2026-08-17)

Implemented in worktree `wt-450` on branch `fix/450-inert-clip-refusal`, based explicitly on
`origin/dev` @ `4f7b907` and verified with `git worktree list` before any edit. **Not committed** —
awaiting review. 8 files, **+438/−88**, **788 tests passing** (781 base + 7 new), `fmt` clean,
`clippy --all-targets --release -D warnings` clean.

All 18 outline steps landed. Every step's file and line matched the plan, which is what the corrected
provenance note bought — no line-number surprises during implementation.

### Verification results

| # | Result |
|---|---|
| V1 | **Reproduced at base by execution.** `--clip_R2 5 --rename --output-format ubam` on a one-record SE FASTQ: warning present, refusal present, **exit 1**, no BAM written |
| V2 | **Confirmed, with a positive control.** After the change: exit 0, warning kept, refusal absent, `sp_trimmed.bam` written. Its BAM records are **byte-identical** to a `--rename`-only run, read name `withspace`, no `:clip5:`/`:clip3:` anywhere. A third run that genuinely clips (`--clip_R1 5`) reports DIFFER against the same comparator, so the IDENTICAL verdict is not vacuous |
| V3 | All twelve truth-table rows are tests; the five new test functions were each confirmed to have run **by name**, not inferred from a total |
| V4 | **Tighten** (`any_effective` counts Read 2 always ≡ `any_set`) → `rename_into_ubam_accepted_for_inert_read_2_clip_flag` **red**, 2 failures total |
| V5 | **Loosen** (`applies` drops `paired \|\|`) → `rename_into_ubam_still_refused_where_a_clip_can_apply` **red**, 6 failures total |
| V6 | `read_2()` false for `ClipR2` → the oracle **red**, 11 failures total |
| V6b | `read_2()` true for `ThreePrimeClipR1` → the oracle **and** row 4 **red**, 7 failures total |
| V7 | **27 cells, 102 artifacts compared, 1 differing cell, 0 vacuous.** The one difference is `se_rename_ubam_r2_only` — the #450 shape. Its three nearest neighbours (`--clip_R1` SE, `--clip_R2` paired, `--library emseq` SE) are exit 1 on both sides with identical stderr |
| V8 | `cargo test --no-fail-fast` 788/0, `fmt --check` clean, `clippy` 0 warnings — all re-run on the final tree after rustfmt |
| V9 | Not run locally (CI job). `Cli::validate()` was edited, so this is the one to watch on the PR |

**V4 and V5 fail different tests, which is the empirical proof they are no longer the same mutation.**
That was v1's defect: as written, V4 would have printed a green run and the honest reading of it was
"the narrowing is not load-bearing".

**V5's blast radius is the clearest evidence the fold-in worked.** Breaking `applies` alone turned red
six tests across `flag_summary`, the paired report, the override announcements and the #408 guard —
five sites failing from one edit is what "one shared predicate" means. Under v1's three-site scope the
same mutation could only have reached one of them.

### Deviations from the plan

| # | Deviation | Why |
|---|---|---|
| D1 | **Step 10's comment is two lines, not one.** | One line could not carry both the criterion and the siting, and dropping the criterion would re-lose the fact #450 turned on. Two is CLAUDE.md's maximum and this is its stated reason. The `append_to_id` site list the old comment carried moves to the commit message |
| D2 | **V7 needed two normalisations the plan did not anticipate.** Reviewer A checked that no report carries a timestamp and concluded byte-comparison was feasible; two things still differ by construction. (a) stderr carries the `--version` provenance line, and the two binaries were built minutes apart. (b) both reports record `command_line`, which begins with **argv[0]** — the binary's own absolute path, i.e. which worktree it was built in. Relative input paths, which the plan did specify, do not cover argv[0]. Both are erased by an explicit two-rule normaliser that names what it erases; nothing else is touched | Without it every one of the 27 cells "differed" and the differential proved nothing |
| D3 | **V7 gained a `CONTROLS` set.** The plan's non-degeneracy rule flagged 7 cells as vacuous. Three were refusal cells whose real marker (`--rename is refused`) was missing from the marker list; four are cells where every flag given is live, so no site under test is reached — regression ballast, not evidence. Naming them beats letting them fail the rule silently | The rule was right and its marker list was incomplete; a "7 of 27 prove nothing" report understated the differential |
| D4 | **`pick`'s fourth parameter is gone rather than converted.** The plan said change the first parameter to `ClipFlag` and drop the fourth; worth recording that this makes #447's "set `read_2` at the call sites, not from the flag string" decision obsolete rather than merely satisfied — the variant now carries the fact, so passing a wrong bool is unrepresentable | Strictly stronger than the invariant it replaces |
| D5 | Added `any_effective_discounts_read_2_on_single_end` to `library.rs`, not in the outline | The behaviour change deserved a unit test at the predicate, not only integration coverage. It is red under V4, V5, V6 and V6b — the most sensitive of the new tests |

### Iteration log

- **#1** `cargo build` filtered through `grep` returned empty output and a blank exit status. Re-ran writing to a log and checking `$?` directly: 0 errors, 0 warnings. A filtered build log is not a build result.
- **#2** `cargo fmt --check` failed on `pick`'s reformatted closure. Ran `cargo fmt --all`; re-check clean. rustfmt also rewrapped three of my new asserts.
- **#3** First V2 attempt was denied by the permission gate for `rm -rf`. Re-ran with a fresh directory name instead — `unlink` and unique paths, per CLAUDE.md.
- **#4** V7 round 1: all 27 cells differed, on the build-timestamp provenance line. Inspected one JSON diff before normalising anything and found the second, unrelated cause (`command_line` carrying argv[0]). Normalising blind would have hidden whichever cause I had not yet found.
- **#5** V7 round 2: 1 differing cell, 7 vacuous. Refined the marker list and added `CONTROLS` (D3). Round 3: 1 differing cell, 0 vacuous.
- **#6** Mutation matrix run as a script that snapshots, mutates, tests and restores per mutation, then asserted the pristine predicates are back and the snapshot file is gone — a subagent in an earlier session left a mutation applied while reporting a clean tree.
