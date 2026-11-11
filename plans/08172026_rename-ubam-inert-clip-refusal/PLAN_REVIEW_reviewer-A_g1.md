# Plan review — Reviewer A (generation 1)

**Plan:** `plans/08172026_rename-ubam-inert-clip-refusal/PLAN.md`
**Reviewed against:** `origin/dev = 4f7b907` — the plan's stated base, unchanged throughout the review.
Every source claim was checked via `git show origin/dev:<file>` into `$TMPDIR/rev/`, never the working
file.
**Constraints honoured:** read-only; no `cargo build` / `test` / `clippy`; no branch, checkout or stash
operations. Three items could only be settled by compiling or running the binary and are listed as open
verification items at the end.

> **The shared tree moved twice while this review was in progress** — `HEAD` went
> `1b2579e` → `4e47b53` → `992235e`, and the branch was renamed from `feat/clumpify-guards-439` to
> `ci/clumpify-memory-guard-439`. None of the verdicts below are affected, because all of them were
> taken from `origin/dev` temp copies rather than the working file. See finding **I6**, which the moves
> strengthen rather than change.

---

## 1. Verification log — the plan's eight load-bearing claims

One line per item, with the file:line I read.

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1 | `!cli.paired` ⟹ no `trim_read(…, is_r2 = true)` ⟹ Read 2 clip values never read ⟹ no `append_to_id` from a Read 2 flag | **TRUE** | All 9 production `trim_read` sites enumerated; the 3 passing `true` are `trimmer.rs:480` (`run_paired_end`), `trimmer.rs:731` (`run_paired_end_to_bam`), `parallel.rs:495` (`classify_paired`) — every one in a paired function, and all three paired entry points are called only under `if cli.paired` (`main.rs:862`, `main.rs:893`, `main.rs:2046`/`2072`). `config.clip_r2` / `three_prime_clip_r2` are read for trimming at exactly `trimmer.rs:247` and `:252`, inside `trim_read`; every other mention in `src/` is report plumbing. `specialty.rs`'s four `--rename` appends (`:62`, `:105`, `:164`, `:217`) are hardtrim, covered by the guard's own hardtrim terms; `specialty.rs:310-311`/`:409-410` are `--clock`/`--implicon` UMI appends, unconditional but rejected for uBAM output at `cli.rs:630-640`; `demux.rs` contains no `append_to_id`. The enumeration of `is_r2 = true` sites is complete. |
| 2 | The #408 guard is reachable on a single-end run — not inside `if cli.paired` | **TRUE** | Guard at `main.rs:389-415` sits at 4-space indent directly in `main()` (`main.rs:217`); the preceding `if cli.paired` block opens at `:372` and closes at `:387`. `cli.validate()` runs first at `:242`, which is why the warning prints before the bail. |
| 3 | The two display sites emit the four flags in different orders | **TRUE** | `flag_summary` array `library.rs:124-129`: `clip_R1, clip_R2, three_prime_clip_R1, three_prime_clip_R2`. Ladder array `cli.rs:1154-1166`: `clip_R1, three_prime_clip_R1, clip_R2, three_prime_clip_R2`. No shared order preserves both. |
| 4 | All five presets set a non-zero `clip_r1`, so every preset keeps the refusal on single-end | **TRUE** | `library.rs:69-75`: EmSeq 10, Accel 10, Zymo 10, ScBs 6, Pbat 8. |
| 5 | No existing test asserts the behaviour being changed | **TRUE** | `--clip_R2` / `--three_prime_clip_R2` appear nowhere in `tests/integration_ubam_out.rs`. The four refusal tests use `--clip_R1` (`:1006`), `--library emseq` (`:1044`), `--hardtrim5` (`:1172`), `--hardtrim3` (`:1204`); the two negative ones (`:1232`, `:1263`) use `--clip_R1`. All six still behave the same after the change. |
| 6 | `value.is_some() && (all_four \|\| flag.read_2())` is arm-by-arm equivalent to the shipped nested construction | **TRUE** | `cli.rs:1154-1172`. R1 arms: shipped `if all_four { self.clip_r1 } else { None }` warns iff `all_four && is_some()`; proposed `is_some() && (all_four \|\| false)` — identical. R2 arms: shipped `self.clip_r2` warns iff `is_some()`; proposed `is_some() && (all_four \|\| true)` — identical. The printed `{n}` is the same raw value in both. |
| 7 | `any_set()` has one caller and no test references it | **TRUE** | `git grep any_set` across `origin/dev` returns exactly `library.rs:141` (definition) and `main.rs:398` (sole caller). No test mentions it. `pub mod library` in `lib.rs:17` on a `pub` struct means no dead-code warning would fire if it were left behind — step 3's rationale is sound. |
| 8 | `clip_5prime` guards on `n > 0` and `seq.len() > n` | **TRUE** | `fastq.rs:110-119` (`if n > 0 && self.seq.len() > n`), `clip_3prime` the same at `:123-133`. The categorical/statistical distinction the plan draws on top of this is sound. |

Two further checks I ran that support the plan: I confirmed the RRBS auto-clip cannot reach the guard
(`main.rs:1141` is inside `setup_trimming`, which starts at `:1051` and is called at `:880`, `:988`,
`:1038`, `:2059`, `:2117`, `:2168` — all after the guard at `:396`), and I confirmed neither report nor
stderr carries a timestamp, elapsed time or git hash (`config.version` is `CARGO_PKG_VERSION`,
`report.rs:205`/`:801`), so V7's byte-comparison is feasible in principle.

**Every load-bearing claim holds.** The findings below are about scope, validation design, and three
things the plan states about its own environment and about the docs that are not true.

---

## 2. Critical

### C1 — V4's mutation cannot produce V4's expected result, and V4 and V5 are the same mutation

`any_effective(paired)` is `… && (paired || !flag.read_2())`. V4 says "drop `paired ||`", which makes
the predicate `value.is_some() && !flag.read_2()`: a Read 2 value then **never** counts.

- On single-end `--clip_R2 5`, that mutation makes the guard **not** fire — which is exactly what
  step 9's two new tests assert. They stay **green**. V4's stated expectation ("V2's two tests go
  red") is wrong.
- V5's mutation, "make `any_effective` skip Read 2 values unconditionally", *is* that same mutation,
  with the right expectation attached (the paired rows go red).

So the plan's two core pieces of evidence collapse into one, and the direction it most needs to prove
— that the narrowing is load-bearing rather than vacuous — has no mutation behind it. The mutation
that turns step 9's tests red is the opposite one: **make `any_effective` count Read 2 values
regardless of `paired`**, i.e. reinstate `any_set`'s behaviour. Rewrite as:

| # | Mutation | Expected |
|---|---|---|
| V4 | `any_effective` ignores `paired` and counts Read 2 values always (≡ `any_set`) | step 9's two single-end tests go red |
| V5 | drop `paired \|\|` so Read 2 values never count | step 10's two paired tests go red |

This matters beyond bookkeeping: run as written, V4 reports "tests stayed green", and the honest reading
of that is "the narrowing is not load-bearing" — the opposite of the truth.

### C2 — Step 13's docs check is a false negative. The refusal is documented in three places and all three become inaccurate

`grep -rn "rename is refused" docs/` greps for the **error string**, which appears only in
`src/main.rs:404`. The *behaviour* is documented, in the user's own words, at three places I found by
grepping for `--rename` instead:

- `docs/src/content/docs/guide/outputs.md:100` — "`--rename` — only when at least one input is FASTQ,
  and **only when a clipping flag is set**; uBAM input is accepted"
- `docs/src/content/docs/guide/outputs.md:145` — "`--output-format ubam` rejects it when **at least one
  input is FASTQ** and **a clipping flag is set**: …"
- `src/cli.rs:312-313` — the `--help` text: "Refused with `--output-format ubam` when any input is
  FASTQ" — which is *already* over-broad today (`--rename` alone is accepted; see
  `rename_without_clip_flag_is_accepted_into_ubam`, `integration_ubam_out.rs:1078`), and this change
  widens the gap.

Step 13 as written ("this may be a no-op; if the refusal is undocumented, note it as a gap and do
**not** open a new docs section") instructs the implementer to conclude the refusal is undocumented and
ship two published sentences stating the superseded rule. The docs site auto-deploys from `dev`.

There is no `--help` snapshot test anywhere in `tests/`, so `cli.rs:312` will drift silently and
permanently.

Replace step 13 with the three exact locations and the corrected condition — "a clipping flag that can
apply to this run", or whatever wording Felix prefers. No new section is needed; all three are
one-clause edits.

---

## 3. Important

### I1 — Read ownership lives at five production sites, not three; the plan rewrites the least-tested one and leaves the two best-tested alone

Missing from the plan's table: **`ClipOverride`** (`library.rs:87-94`), whose two fields are
`flag: &'static str` and `read_2: bool` — precisely the pair `ClipFlag` is being introduced to own. It
is populated by four `pick(…)` calls each passing a hardcoded string and a hardcoded bool
(`library.rs:187-200`), and `read_2` is consumed as `filter(|o| paired || !o.read_2)` at **two**
production sites: `main.rs:1128` (the banner) and `report.rs:287` (the text report). `o.flag` is also
written straight into the JSON report at `report.rs:898`.

After the plan lands as scoped:

- `library.rs` holds **two** independent read-ownership tables 60 lines apart, one of which (`pick`'s
  call sites) is still a hand-written bool literal per flag — which is where the drift risk actually
  lives, since a `ClipFlag` argument would make it unrepresentable.
- The predicate `paired || !read_2` is spelled **four** times: `library.rs:130` (negated),
  `main.rs:1128`, `report.rs:287`, and the new `any_effective`.
- Step 8a's drift test duplicates the existing `every_override_knows_which_read_it_belongs_to`
  (`library.rs:296-317`), which already uses `flag.ends_with("_R2")` as the oracle over the same four
  flags. The repo ends with two near-identical drift tests instead of one.

Concretely: make `ClipOverride.flag` a `ClipFlag` and drop its `read_2` field; give `ClipFlag` one
predicate `applies(self, paired) -> bool` (= `paired || !self.read_2()`) that serves `any_effective`,
`flag_summary`'s skip, and both override filters; `report.rs:898` becomes
`json_escape_string(o.flag.name())`. That is the refactor the plan's secondary goal describes.

**The risk ordering is inverted by the current scope.** The sites the plan leaves alone are the
protected ones:

- `flag_summary`'s exact string *and* order is pinned three times — `library.rs:322-325` and `:331-334`,
  `tests/integration_library.rs:272`, `src/report.rs:1399`.
- The override filter is pinned in both directions — `single_end_announces_no_read_2_clip`
  (`integration_library.rs:339`) and `paired_still_announces_r2_overrides` (`:371`).
- The warning ladder's **order is pinned nowhere** (see I2).

So the plan rewrites the only unprotected display site and skips the two that would tell it immediately
if it got them wrong.

### I2 — The ladder's emission order has no test, and V7 is a one-shot manual diff that leaves none

`each_inert_flag_gets_its_own_line` (`tests/integration_inert_clip_flags.rs:219-240`) counts warning
lines and asserts both flag names appear. It never asserts that `--clip_R2` precedes
`--three_prime_clip_R2`. And no test anywhere sets all four clip flags at once:
`inert_clip_flag_warns_with_a_mode_accurate_reason` (`:60-154`) passes `R1_PAIR` on the four specialty
rows and `R2_PAIR` on the single-end row — two flags per row, never four — so the ladder's full
4-element order is **entirely unexercised** today.

The plan gives its longest section to preserving that order, then protects it with a pre-merge
differential that costs two release builds (~4 GB of `target/`, against the plan's own "89% used" note)
and leaves nothing behind for the fourth visit to this cluster that Remaining Risk 3 anticipates.

Cheap and permanent instead: add to `each_inert_flag_gets_its_own_line` an assertion that
`err.find("--clip_R2 5") < err.find("--three_prime_clip_R2 7")`, and add one all-four case under
`--hardtrim5` asserting the four positions ascend. That converts V7 from load-bearing evidence into an
optional belt-and-braces check.

### I3 — The plan never names `tests/integration_inert_clip_flags.rs`

V7 refers to "the #447 matrix" abstractly. The file exists, is 240 lines, and pins the ladder across
five modes plus four quiet/negative cases including
`preset_alone_on_single_end_does_not_warn` (`:200`) — the test that enforces Assumption 3's
raw-vs-resolved trap. The implementer should be told by name: this file must stay green **unmodified**
by step 6, and it is the natural home for I2's order assertions.

### I4 — V3 claims all 10 truth-table rows as tests; steps 9-11 deliver 8. `--three_prime_clip_R1` is tested nowhere

Coverage after the plan:

- row 3 `--clip_R1` single-end — existing `:1006`
- row 7 `--library` — existing `:1044` (emseq) + step 11 (all five)
- row 8 hardtrim5/3 — existing `:1172`, `:1204`
- row 9 `--rename` alone — existing `:1078`
- row 10 uBAM input — existing `:1109`
- rows 1, 2 — step 9; rows 5, 6 — step 10
- **row 4 `--three_prime_clip_R1` — nowhere, before or after.** `git grep -F -- "three_prime_clip_R1" -- tests` never hits `integration_ubam_out.rs`.

Step 9's own rationale cites "the #447 lesson that a table testing one representative per branch misses
the second entry" — and then omits the second entry of the R1 pair. One extra invocation closes it, and
it is the arm a mutation to `read_2()` for `ThreePrimeClipR1` would otherwise walk straight through.

### I5 — V7's matrix contains impossible cells, and two identical bail messages compare equal

`--library` is refused with `--hardtrim5`, `--hardtrim3`, `--clock` and `--implicon`
(`cli.rs:1115-1135`), so "5 modes × the five presets" contains 20 cells that cannot run;
`--clock`/`--implicon` additionally need the clock fixtures and two files *without* `--paired`
(the shape `integration_inert_clip_flags.rs:12-13` exists to cover).

Two binaries emitting the same "`--library` is not compatible" bail compare **equal**, so those cells
report PASS while testing nothing — a matrix that silently degenerates into comparing error strings.
This is the same class as the harness false-negatives already recorded in this repo: a parseable output
from a run that did nothing looks like data.

Fix: each cell must first assert it produced either a `WARNING: … is not used in this mode` line or a
`Clipping in force:` line, and only then compare the two outputs. State the valid-cell set explicitly
rather than as a product.

### I6 — Assumption 6 and the line-number provenance note are stale, and one is wrong in a misleading way

- **Nothing is uncommitted, and the branch Assumption 6 names no longer exists.** The parallel session
  committed three times during this review — `1b2579e`, `4e47b53`, `992235e` — and renamed its branch
  from `feat/clumpify-guards-439` to `ci/clumpify-memory-guard-439`. `git status --porcelain` shows no
  modified tracked files. So the hazard is not "uncommitted edits in the tree" as Assumption 6
  describes; it is that reading a working file now shows a *newer commit*, which looks like clean
  ground and is not. Assumption 6's remedy (a worktree explicitly based on `origin/dev`, verified with
  `git worktree list`) is still exactly right — only its diagnosis needs restating.
- Those commits touch `src/clump.rs`, `src/clump_only.rs`, `src/fastqc.rs`, `src/main.rs`,
  `src/parallel.rs`, `tests/integration_preflight_tripwire.rs`, `examples/mk_clump_fixture.rs`,
  **`CHANGELOG.md`**, `.github/workflows/ci.yml` and `.github/clumpify-memory-baseline.tsv` — not the
  three files Assumption 6 lists. `CHANGELOG.md` is the one that matters: step 12 edits it too, so a
  conflict is expected once `dev` advances. Worth one line in Assumption 6.
- **The note's claim that `src/parallel.rs` is "unmodified, so their numbers agree in both" is false.**
  `1b2579e` adds 13 lines to it, and the plan's own citations `parallel.rs:496` and `parallel.rs:1114`
  are the **working-tree** numbers. On `origin/dev` they are `495` and `1111`. An implementer following
  the plan's stated convention will read the wrong lines and may "correct" a true claim into a false
  one. (The `main.rs` `+2` offset the note does give is correct: `if cli.rename` is `:398` in the
  working file and `:396` on `origin/dev`, and that held across all three commits.)
- The plan says "all six production call sites". There are nine (three `true`, six `false`); the six
  named are the three `true` sites plus the three single-end `false` sites, omitting the R1 halves of
  the paired calls. Harmless to the argument — the three `is_r2 = true` sites are enumerated
  completely, which is the load-bearing part — but the count is wrong.

Assumption 6's worktree guidance itself is right and worth keeping.

---

## 4. Optional

- **O1 — three more citations off:** the preset table is `library.rs:69-75` (plan: 73-79); the clip
  selection is `trimmer.rs:246-255` (plan: 246-256); the ladder array is `cli.rs:1154-1166`
  (plan: 1154-1167).
- **O2 — `ALL`'s job is stated two ways.** Signature: "no display site iterates this" and "`ALL`
  therefore has exactly one job: letting step 8's drift tests iterate the closed set". Efficiency:
  "`ALL.iter()` over a `const` array compiles to the same straight-line code as the current
  hand-written tuple arrays" — which implies `any_effective` iterates it. Decide: if it does, `ALL` has
  a production consumer and its doc comment must say so; if it does not, `ALL` can be `#[cfg(test)]`.
  As written the implementer guesses.
- **O3 — step 8b claims a compile-time property it does not quite have.** An exhaustive `match` over an
  iterated `ClipFlag` does fail to compile when a fifth variant appears, which forces a visit to the
  test — but nothing forces the fifth variant *into* `ALL`, whose length is the literal `4`. Either
  spell the test out or weaken the claim to "adding a variant forces a visit to this test".
- **O4 — pin `name()`'s exact returns in the Signature,** `"--clip_R1"` with the leading `--`. Three
  user-visible strings compose it: `flag_summary`'s `format!("{flag} {v}")` (`library.rs:134`), the
  ladder's `"WARNING: {flag} {n} …"` (`cli.rs:1169`), and the JSON `overrides[].flag`
  (`report.rs:898`). Step 8a's `ends_with("_R2")` oracle passes with or without the dashes.
- **O5 — step 5's rationale is inaccurate.** There is no "duplicated `UserClips` literal" today;
  `cli.rs:1228` is the only one. The extraction prevents a duplicate that step 6 would otherwise
  create — say that.
- **O6 — the Goal renders the warning as two wrapped lines with a hanging indent.** `cli.rs:1168-1170`
  emits one unwrapped line. Cosmetic, but a spec that shows invented line breaks invites a test that
  asserts them.
- **O7 — V2 needs two output directories or two input stems.** Both runs write `sp_trimmed.bam` from
  `sp.fastq`; `assert_ubam_eq` (`integration_ubam_out.rs:108`) compares two produced files fine, but
  they must not collide first.
- **O8 — step 12's "one sentence" will read as a deliberate contrast.** Every neighbouring Unreleased
  bullet is a paragraph; the #440 entry runs 23 lines. One sentence is right per the global rule —
  just noting it is a visible break with the file's house style, not an omission.
- **O9 — `main.rs:389-395` is already a 7-line comment block** against CLAUDE.md's "one line, two
  maximum". Step 7 adds an eighth. Consider collapsing the block to one line stating the criterion and
  moving the rest to the commit message, which is where the plan already says the mechanism belongs.
- **O10 — `ClipFlag` in `library.rs` is an altitude mismatch** with that module's doc comment
  ("Library-prep presets (`--library <name>`)"): neither the ladder nor the #408 guard has anything to
  do with presets. `UserClips` already lives there, so it is consistent — just widen the module doc in
  the same commit.

---

## 5. Out-of-scope observation worth filing separately

### X1 — The same guard *under*-refuses on the `--rrbs` directional path, and B3 records the mechanism without noticing the asymmetry

`--paired --rrbs --rename --output-format ubam r1.fastq r2.fastq`, with no clip flag and no preset:

1. the guard's clip term is false and both hardtrim terms are false, so the run **proceeds**
   (`main.rs:396-401`);
2. `setup_trimming` then sets `clip_r2 = Some(2)` (`main.rs:1141-1147`) and `config.rename = cli.rename`
   (`main.rs:1232`);
3. the uBAM paired path uses that config (`setup_trimming` at `main.rs:2117` →
   `run_ubam_output_paired_two_files` at `:2311`) — `--rrbs` is not among §3.4a's uBAM exclusions;
4. `trim_read(r2, config, true)` reads `config.clip_r2` (`trimmer.rs:247`) and calls
   `append_to_id(":clip5:…")` (`trimmer.rs:264`);
5. `append_to_id` appends at the **end** of the ID when there is no tab (`fastq.rs:174-182`), so for a
   header carrying a description the annotation lands after the description and the BAM QNAME keeps
   only the first token.

That is the exact hazard the #408 error text describes, in the direction opposite to #450: an
under-refusal where #450 is an over-refusal. It also cuts against the plan's central argument — that the
guard's condition should match the set of reachable `append_to_id` sites — because the RRBS auto-clip is
a reachable site the condition does not cover. B3's row currently reads "It cannot reach the guard in
either direction. No change", which is true about the mechanism and silent about the consequence.

Not for this commit — Q-1 sets the precedent that a neighbouring defect of the same class stays out. But
B3's row should say "under-refuses; tracked separately" rather than "no change", and this clears Felix's
issue-filing bar (silent loss of a user-requested annotation, no error, exit 0). I could not run the
binary to observe step 5 end-to-end, so treat that step as derived from the code path.

---

## 6. Alternatives

### A1 — The behaviour fix is five lines and touches no display site

```rust
pub fn any_effective(&self, paired: bool) -> bool {
    self.clip_r1.is_some()
        || self.three_prime_clip_r1.is_some()
        || (paired && (self.clip_r2.is_some() || self.three_prime_clip_r2.is_some()))
}
```

That is the whole of B1, and it makes V7, the ordering section, and Remaining Risk 1 unnecessary.
Everything else in the plan serves B2 — and B2 as scoped resolves three of five sites (I1) while
carrying all of the display-order risk. There are two coherent positions and the plan sits between
them:

**(a)** ship the predicate now, file the consolidation separately with the full five-site scope; or
**(b)** do the enum properly — `ClipOverride.flag: ClipFlag`, one `applies(paired)` predicate serving
all four consumers, `main.rs:1128` and `report.rs:287` folded in, and the existing
`every_override_knows_which_read_it_belongs_to` **replaced** rather than duplicated.

(b) is the larger diff but lands on protected ground: every site it adds is already pinned in both
directions by a test, unlike the ladder. This is a scoping argument about the refactor, not a
re-opening of narrow-vs-suppress-vs-document, which Felix settled.

### A2 — If the enum stays as scoped, at least make the predicate shared, not just the ownership bit

The plan shares `read_2()` and leaves each consumer to re-spell `paired || !f.read_2()`. Adding
`applies(paired)` costs one method, changes no emission order, and removes the fourth spelling of the
predicate the plan exists to de-duplicate. This is orthogonal to the ordering concern and carries none
of its risk.

---

## 7. Open verification items (could not be settled without building or running)

1. **V8's baseline.** `origin/dev` carries 782 `#[test]` attributes (492 in `src/`, 290 in `tests/`),
   so "≥ 781" is plausible, but I did not run `cargo test --no-fail-fast` to confirm it as a pass
   count.
2. **X1's final step** — the annotation actually vanishing from the BAM QNAME — is derived from the code
   path, not observed.
3. **Step 9's two new tests exiting 0 end-to-end.** Nothing in `Cli::validate()` §3.4a or in `main()`
   blocks `--rename --output-format ubam --clip_R2 5 se.fastq` once the guard is narrowed, but I did
   not run it.

---

## 8. Summary

The plan's technical core is sound: all eight load-bearing claims verified, the narrowing is
well-motivated by the guard's own stated criterion, and the self-review's catch about unifying the
consumption loop rather than the fact is a genuinely good one. Three things need fixing before
implementation:

1. **C1** — V4's mutation is the wrong direction, so the plan's central claim has no evidence behind it.
2. **C2** — the docs check greps for the wrong string; the refusal is documented three times and all
   three become inaccurate.
3. **I1/I2** — the refactor covers three of five sites and rewrites the one display site nothing pins,
   while V7 substitutes a one-shot manual diff for a two-line test assertion.

I4, I5 and I6 are each a small edit to the plan. The Optional list is line numbers and wording.
