# Plan review — Reviewer B (generation 1)

**Plan:** `plans/08172026_rename-ubam-inert-clip-refusal/PLAN.md`
**Base verified against:** `origin/dev` @ `4f7b907`, read via `git show origin/dev:<file>` into pinned
temp copies. No file in the working tree was read for any claim, no build or test was run.
**Verdict:** the behaviour change is correct and its reasoning is sound. Eight of the eight
load-bearing claims hold in substance. One Critical defect (step 13 rests on a false-negative grep and
would ship two contradicting docs sentences), five Important findings, four Optional.

---

## 0. Environment correction — the plan's Assumption 6 is now stale

The parallel session has **committed**. `HEAD` is `1b2579e` (ahead of `origin/dev`) and the working
tree is clean — `git status --porcelain --untracked-files=no` returns nothing. `git diff origin/dev
--stat` shows the changed files as `src/clump.rs`, `src/clump_only.rs`, `src/fastqc.rs`, `src/main.rs`,
`src/parallel.rs`, `tests/integration_preflight_tripwire.rs`.

Two consequences:

- Assumption 6's description ("uncommitted edits to `src/clump.rs`, `src/main.rs` and
  `tests/integration_preflight_tripwire.rs`") no longer matches reality, and the set is wider than
  listed: **`src/parallel.rs` and `src/fastqc.rs` also differ from `origin/dev`**. `src/parallel.rs`
  matters here because the plan cites it for Assumption 1.
- Assumption 6's *conclusion* is unaffected and still correct: implement in a separate worktree
  explicitly based on `origin/dev`, and verify with `git worktree list`. A bare `git worktree add`
  bases off `master`.

Anyone re-checking this plan's `main.rs` / `parallel.rs` line numbers must use `git show
origin/dev:<file>`, not the working file.

---

## 1. The eight load-bearing claims

| # | Claim | Verdict |
|---|---|---|
| 1 | `!cli.paired` ⟹ no `trim_read(…, is_r2 = true)` ⟹ Read 2 clip values never read ⟹ no append site can fire | **TRUE in substance; the citation is wrong in three ways** — see 1.1 |
| 2 | The #408 guard is reachable single-end (not inside `if cli.paired`) | **TRUE** |
| 3 | The two display sites emit the four flags in different orders | **TRUE** |
| 4 | All five presets set a non-zero `clip_r1` | **TRUE** (citation off — see 1.4) |
| 5 | No existing test asserts the behaviour being changed | **TRUE, but the plan omits the one test file step 6 puts at risk** — see 1.5 |
| 6 | `value.is_some() && (all_four \|\| flag.read_2())` ≡ the shipped nested construction | **TRUE, arm by arm** |
| 7 | `any_set()` has one caller, no test references it | **TRUE in-tree**; it is also published public API — see O2 |
| 8 | `clip_5prime` guards on `n > 0` and `seq.len() > n` | **TRUE**, `src/fastq.rs:110-119`, exact |

### 1.1 Claim 1 — mechanism confirmed, enumeration miscounted

What I checked, and it all holds:

- Production `trim_read` call sites, complete: `is_r2 = true` at `src/trimmer.rs:480`,
  `src/trimmer.rs:731`, `src/parallel.rs:495`; `is_r2 = false` at `src/trimmer.rs:319`,
  `src/trimmer.rs:479`, `src/trimmer.rs:645`, `src/trimmer.rs:730`, `src/parallel.rs:1111`.
- The three `true` sites sit in `run_paired_end` (`trimmer.rs:381`), `run_paired_end_to_bam`
  (`trimmer.rs:704`) and `classify_paired` (`parallel.rs:479`). Every entry into those is gated on
  `cli.paired`: `main.rs:862` → `run_paired_ubam_single_file`, `main.rs:893` → `run_paired`,
  `main.rs:2046` → `run_ubam_output_paired_single_file`, `main.rs:2072` →
  `run_ubam_output_paired_two_files`; `run_paired_end_parallel` is reached only from `main.rs:1606`
  and `main.rs:1885`, both inside those paired functions.
- `config.clip_r2` / `config.three_prime_clip_r2` are read in production **only** at
  `src/trimmer.rs:247` and `:252`, inside `trim_read` under `is_r2`. Every other hit is config
  construction (`main.rs:1229/1231`) or report serialisation (`main.rs:1491/1493`, `:2279/2281`,
  `:2666/2668`, `report.rs`), none of which clips.
- The briefing asked me to consider `specialty.rs`, `demux.rs`, the uBAM output paths and the
  interleaved-uBAM paired path. `git grep 'clip_r1\|clip_r2\|trim_read\|TrimConfig' -- src/specialty.rs
  src/demux.rs` returns **nothing** — neither module touches a clip flag. The `:clip5:`/`:clip3:`
  append sites are exactly `trimmer.rs:264`/`:273` plus `specialty.rs:62/105/164/217` (hardtrim,
  covered by the guard's own `hardtrim5/hardtrim3` terms). `specialty.rs:310/311/409/410` append a
  clock/implicon UMI, not a clip, and both modes are refused with `--output-format ubam`
  (`cli.rs:632`, `cli.rs:638`), so they cannot reach the guard. The uBAM output paths route through
  the same `trimmer::run_*` functions. **The enumeration of `true` sites is complete.**

Three citation errors, none of which changes the conclusion but one of which sits in the plan's most
load-bearing sentence:

- `parallel.rs:496` should be **`:495`**.
- `parallel.rs:1114` should be **`:1111`**.
- Assumption 1 says "Verified against all six production call sites". There are **eight** — the plan
  omits `trimmer.rs:479` and `trimmer.rs:730`, the R1 halves of the two paired calls. The three `true`
  sites are what the claim needs, and those are all present, so the argument survives; but a sentence
  whose whole force is "I enumerated the closed set" should have the count right, because the next
  reader's only way to re-verify is to reproduce the enumeration.

### 1.4 Claim 4 — true, citation misses three of the five presets

`library.rs:69-75` is the value table: EmSeq `(10,10,10,10)`, Accel `(10,15,10,10)`, Zymo
`(10,10,10,10)`, ScBs `(6,6,6,6)`, Pbat `(8,8,8,8)`. All five have `clip_r1 ≠ 0`, so every preset keeps
the refusal on a single-end run. The plan cites `library.rs:73-79` in both B3 and Assumption 2; that
range starts at ScBs and runs into the `PresetClips` struct literal, so a reader following the
citation sees two of the five presets and none of the three that are `10`.

### 1.5 Claim 5 — true for #408, but the plan never names the file that pins step 6

The four #408 refusal tests are exactly where the plan says: `tests/integration_ubam_out.rs:1006`
(`--clip_R1 3`), `:1044` (`--library emseq`), `:1172` (`--hardtrim5`), `:1204` (`--hardtrim3`). No test
anywhere combines `--clip_R2`/`--three_prime_clip_R2` with `--rename`. Claim confirmed.

But the plan's Validation section never names **`tests/integration_inert_clip_flags.rs`** (240 lines,
5 tests), which is the dedicated test file for the warning ladder step 6 rewrites. V7 gestures at it as
"the #447 matrix"; the file is never identified. That matters for I3 below, because what the file does
and does not assert is the whole question for step 6.

Reassuringly, `preset_alone_on_single_end_does_not_warn` (`:201`) directly pins Assumption 3 — a
`--library emseq` single-end run must emit no inert warning, which is exactly what would break if the
ladder were switched to `effective_clips()`. The plan's Assumption 3 is therefore already
test-enforced, which is stronger than the plan claims for it.

### 1.6 Claim 6 — equivalent, checked entry by entry

Shipped (`cli.rs:1154-1172`) vs proposed, per entry:

| Entry | Shipped warns iff | Proposed `value.is_some() && (all_four \|\| read_2())` |
|---|---|---|
| `--clip_R1` | `all_four && clip_r1.is_some()` | `clip_r1.is_some() && (all_four \|\| false)` ✓ |
| `--three_prime_clip_R1` | `all_four && tpc_r1.is_some()` | `tpc_r1.is_some() && (all_four \|\| false)` ✓ |
| `--clip_R2` | `clip_r2.is_some()` | `clip_r2.is_some() && (all_four \|\| true)` ✓ |
| `--three_prime_clip_R2` | `tpc_r2.is_some()` | `tpc_r2.is_some() && (all_four \|\| true)` ✓ |

The printed line is `{flag} {n}`, identical under both. Cross-checked against the three tests that
constrain the ladder (`inert_clip_flag_warns_with_a_mode_accurate_reason`,
`read_1_clip_flags_do_not_warn_on_single_end`, `each_inert_flag_gets_its_own_line`) — all four arms
behave the same. The plan's instruction to verify by reading *and* by V7 is right, and the reading
holds.

---

## 2. Critical

### C1 — Step 13 rests on a false-negative grep. The refusal *is* documented, in wording the change makes wrong.

Step 13 says:

> `grep -rn "rename is refused" docs/` currently returns nothing, so this may be a no-op; if the
> refusal is undocumented, note it as a gap and do **not** open a new docs section in this commit.

The grep pattern is the **error-message** wording. The docs use "rejects it". Searching for the
concept instead of the string finds two live sentences, both of which state the clip-flag gate and
both of which the narrowing falsifies:

- `docs/src/content/docs/guide/outputs.md:100` —
  "`--rename` — only when at least one input is FASTQ, and **only when a clipping flag is set**; uBAM
  input is accepted"
- `docs/src/content/docs/guide/outputs.md:145` —
  "And `--output-format ubam` rejects it when **at least one input is FASTQ** and **a clipping flag is
  set**: …"

After the change the condition is not "a clipping flag is set" but "a clipping flag that can reach a
read on this run's shape is set". A single-end run with only `--clip_R2` satisfies the documented
condition and will no longer be refused. Both sentences become false.

There is a **third** surface the plan does not consider at all, because step 13 only looks in `docs/`
— the `--help` text at `src/cli.rs:312-313`:

> "Refused with --output-format ubam when any input is FASTQ, because a BAM read name cannot hold the
> annotation once the header carries a description"

This one omits the clip-flag condition entirely, so it is **already** an overstatement today —
`rename_without_clip_flag_is_accepted_into_ubam` (`integration_ubam_out.rs:1078`) proves `--rename`
alone is accepted. The change widens an existing gap rather than opening a new one, so I would not
insist it be fixed here; but the plan should decide about it explicitly rather than never see it. Note
that the error message at `main.rs:404` shares the same omission, and Assumption 5 already rules the
error text out of scope by Felix's decision — which means `outputs.md` is currently the **only** place
the clip-flag gate is written down anywhere. That makes step 13 load-bearing, not a probable no-op.

**Ask:** replace step 13 with a named edit to `outputs.md:100` and `:145` (one sentence each), and add
a one-line decision on `cli.rs:312-313`. Re-run the discovery with a concept-level pattern, e.g.
`git grep -n -i -E 'rename.*(refus|reject)|(refus|reject).*rename|clipping flag is set' -- docs
README.md`, and record the pattern in the plan so the next visit does not repeat the miss.

This is the failure mode `feedback_corrective_copy_defect` describes: the claim ("docs are silent") and
the artifact (two sentences that say otherwise) are two separate things, and only the claim was
checked.

---

## 3. Important

### I1 — Read ownership is encoded in **six** places, not three. The secondary goal is not met as stated.

The plan's table lists three sites and B2 promises `ClipFlag` "becomes the single place that says which
read a clip flag belongs to". Three more sites will still say it independently after the change:

| Site | Code | What it encodes |
|---|---|---|
| 4 | `ClipOverride.read_2` field, `library.rs:93`, populated by `pick(…, false/true/false/true)` at `library.rs:187`, `:188`, `:193`, `:199` | The same fact, as four hard-coded literals |
| 5 | `main.rs:1128` — `for o in clips.overrides.iter().filter(\|o\| cli.paired \|\| !o.read_2)` | The consumer predicate |
| 6 | `report.rs:287` — `.filter(\|o\| config.paired \|\| !o.read_2)` | The same consumer predicate, again |

Three things make this more than a completeness nit:

1. **`ClipOverride.read_2`'s doc comment is already the plan's new doc comment.** `library.rs:92-93`
   reads "A Read 2 flag, which a single-end run cannot honour"; the plan's `ClipFlag::read_2()` reads
   "A Read 2 flag reaches no clip site on a single-end run." After the change a reader finds two
   public items in the same module asserting the same fact in nearly the same words. That is the exact
   confusion the refactor exists to remove.
2. **`ClipOverride.flag: &'static str`** (`library.rs:89`, "Flag name as the user would type it, e.g.
   `--clip_R1`") duplicates the new `ClipFlag::name()` table, with the four literals repeated a third
   time at `library.rs:187-199`.
3. **`library.rs:296-317`, `every_override_knows_which_read_it_belongs_to`, already implements step
   8a's oracle** (`o.read_2 == o.flag.ends_with("_R2")`) for `ClipOverride`. Step 8a duplicates an
   existing test's idea for a new type while leaving the old one in place.

Sites 5 and 6 also show that the predicate the plan needs is already written twice in the codebase, in
precisely the form `paired || !read_2`. That argues the shared item should be a predicate, not just an
accessor:

```rust
/// A flag applies on a run of this shape.
pub fn applies(self, paired: bool) -> bool { paired || !self.read_2() }
```

with `any_effective` becoming `ALL.iter().any(|f| f.applies(paired) && self.value(f).is_some())`, and
sites 5 and 6 calling `f.applies(paired)` instead of open-coding it. As proposed, `read_2()` alone
leaves each of the four consumers to re-derive `paired || !read_2` for itself.

**The output-neutral consolidation is available and small.** Changing `ClipOverride` to
`{ flag: ClipFlag, preset_value, user_value }` needs `o.flag.name()` at `report.rs:898` and
`main.rs:1131`, and `o.flag.applies(paired)` at `main.rs:1128` / `report.rs:287`. I checked the JSON
writer (`report.rs:893-902`): it emits only `flag`, `preset_value`, `user_value` — **`read_2` is not
serialised**, so this changes no output byte, and `every_override_knows_which_read_it_belongs_to`
becomes unnecessary because step 8a covers the same table.

**Ask:** either fold sites 4-6 in (my recommendation — it is output-neutral and it is what B2 already
promises), or add a Q-4 that scopes them out explicitly, the way Q-1 scopes out `--clip_R1 0`. What
should not happen is silence: the plan is otherwise scrupulous about naming what it declines, and an
unstated omission here will read as an oversight to the next visitor, who is the person remaining-risk
#3 is written for.

### I2 — `any_effective`'s proposed docstring overstates its predicate, in the same shape as the defect being fixed

The Signature gives it as:

> `/// True when a clip value in force can reach a clip site.`

That is false on a shipped path. Under `--hardtrim5 30 --clip_R1 3`, no clip flag reaches a clip site —
`specialty.rs` reads no clip flag at all (verified: `git grep 'clip_r1\|clip_r2' -- src/specialty.rs`
is empty), and the warning ladder says so at `cli.rs:1139-1140`. Yet `any_effective(false)` returns
`true`. The refusal is still correct, but only because the guard's separate `hardtrim5.is_some()` term
carries it — the method's stated contract is violated and the violation is masked.

The plan's own strongest argument is that "the guard's own comment is the specification its condition
violates". Shipping a new method whose docstring overstates its condition reproduces that pattern one
level down, and the next person to reuse `any_effective` outside the guard will trust the docstring.

**Ask:** narrow the doc to what the method actually decides — read shape only, e.g. "True when a clip
value in force belongs to a read this run processes. Says nothing about whether the *mode* honours
clip flags; the caller's `hardtrim` terms cover that." One line, and it makes the method safe to reuse.

### I3 — The warning ladder's order is pinned by **no** committed test; `flag_summary`'s is. The plan's top risk is guarded only by a one-off.

The plan's remaining risk #1 is that the per-site ordering fix lands at one site and not the other, and
its mitigation is V7, a manual differential run once at implementation time. I checked what CI would
catch afterwards, and the two sites are asymmetric:

- **`flag_summary` order is pinned four ways.** `library.rs:320` asserts the exact string
  `"--clip_R1 10 --clip_R2 15 --three_prime_clip_R1 10 --three_prime_clip_R2 10"`; `library.rs:329`
  asserts the single-end form; `tests/integration_library.rs:196` asserts
  `"Clipping in force: --clip_R1 8 --three_prime_clip_R1 8"`; `integration_library.rs:275` asserts the
  full four-flag string through the text report. A reordering here fails CI immediately.
- **The warning ladder's order is pinned nowhere.** `each_inert_flag_gets_its_own_line`
  (`integration_inert_clip_flags.rs:220-240`) counts matching lines and then calls
  `err.contains("--clip_R2 5")` and `err.contains("--three_prime_clip_R2 7")` — order-agnostic.
  `inert_clip_flag_warns_with_a_mode_accurate_reason` (`:60-154`) likewise uses `contains` plus a line
  count. **A swap of the two warning lines passes the entire suite.**

So the site the plan is most worried about is the one with no regression test, and V7 evaporates the
moment it is run.

**Ask:** make the ordering a committed assertion rather than a one-off measurement. In
`each_inert_flag_gets_its_own_line`, collect the matching lines and assert the `--clip_R2` line precedes
the `--three_prime_clip_R2` one — or, better, add a two-flag cell to
`inert_clip_flag_warns_with_a_mode_accurate_reason`'s R1 set under a specialty mode so both orders are
pinned in one place. That converts the plan's top self-declared risk from "measured once" to "held by
CI", which is what the plan already achieved for `flag_summary` and should match.

### I4 — No proposed test sets two clip flags of *different* read ownership

B1's truth table has ten rows and V3 turns each into a test, but no row sets an R1 and an R2 flag at
once. The missing row is `--clip_R1 3 --clip_R2 5` single-end → **must still refuse**, and it is the
only row whose failure mode is the original #408 data-integrity defect (silent annotation loss) rather
than an over-refusal.

It matters because a natural mis-reading of the narrowing — "skip the refusal when an inert Read 2 flag
is present" rather than "do not *count* an inert Read 2 flag" — produces exactly that shape: it would
let `--clip_R1 3 --clip_R2 5` run, and `:clip5:` would then be appended for R1 into a BAM read name.
Neither of the plan's new R2-only tests has an R1 flag, and `:1006` has no R2 flag, so neither catches
it.

Not Critical, because step 11's five-preset loop covers the predicate today: a preset populates all
four values, so `--library emseq` single-end *is* an "R1 live + R2 inert" case and must refuse. But
that coverage is incidental — it depends on how presets resolve, not on the flag pair — and by the
plan's own #447 lesson an explicit row costs one test.

**Ask:** add an eleventh row and a test for it.

### I5 — `Cli::user_clips()` is a new public bypass of the `effective_clips()` contract

`effective_clips()`'s docstring (`cli.rs:1219-1224`) is emphatic: "Every consumer of the clip values
goes through here, so a preset cannot reach one code path and miss another." Assumption 3 correctly
warns the implementer not to route the ladder through it. But the fix is to add a **public** method on
`Cli` returning the raw pre-preset values, which is a door in the other direction: a future consumer
that reaches for `user_clips()` silently loses the preset, and that is a wrong-output bug rather than a
cosmetic one.

**Ask:** give `user_clips()` a doc comment that scopes it — something like "The values as typed. For
the inert-flag warning only, which must not report flags a preset supplied; every consumer of the
values *in force* goes through `effective_clips()`." Per CLAUDE.md the constraint belongs on the item
it constrains, which is `user_clips()`, not Assumption 3 of a plan file. Consider `pub(crate)` if
nothing outside the crate needs it.

---

## 4. Optional

- **O1 — citation slips.** `parallel.rs:496` → `:495`; `parallel.rs:1114` → `:1111`; "all six
  production call sites" → eight; `library.rs:73-79` → `:69-75`. Also worth naming
  `tests/integration_inert_clip_flags.rs` explicitly in V7 instead of "the #447 matrix".
- **O2 — `any_set` is published public API.** `src/lib.rs:17` is `pub mod library;` and `Cargo.toml`
  publishes `trim-galore` v2.3.0 to crates.io, so removing `ResolvedClips::any_set` removes a public
  item. Step 3's justification ("exactly one caller and no test references it") is an in-tree argument
  for an out-of-tree surface. I agree with removing it — the alternative leaves a correctly-named
  method with the defective behaviour, which is the worse hazard, and the lib target exists for tests
  rather than for downstreams. Just say so in one clause so the decision is visible.
- **O3 — `read_2()`'s doc comment states a consequence, not the fact.** "A Read 2 flag reaches no clip
  site on a single-end run" is what `any_effective` concludes; `read_2()` itself only knows which read
  the flag names. Splitting them keeps each comment true of its own item, and avoids the near-verbatim
  collision with `ClipOverride.read_2`'s existing doc (I1).
- **O4 — V8's test count is plausible.** 782 `#[test]` attributes across `src` and `tests` at
  `origin/dev`, consistent with "≥ 781". Not confirmable without running cargo.

---

## 5. Judged on its own terms

**What the plan gets right, specifically.** The behaviour change is correct and the argument for it is
the strongest part of the document: the categorical/statistical distinction is exactly the right
boundary and is backed by the real guard in `clip_5prime` (`n > 0 && seq.len() > n`); the appeal to
`rename_without_clip_flag_is_accepted_into_ubam` as established precedent is apt rather than
rhetorical; and the Self-Review's caught defect — that unifying the *iteration* would have silently
reordered user-visible output — is a real defect, correctly diagnosed, and the resulting rule
("unifying a duplicated fact is safe; unifying the loop that consumes it is not") is the reason the
`ClipFlag` design is the right shape rather than the wrong one. Q-1 is held out cleanly with the
temptation named. B3's dismissal of the `--rrbs` case is right for both reasons it gives: the
auto-clip is gated on `cli.paired` (`main.rs:1141`) and lives in `setup_trimming` (`main.rs:1051`),
whose earliest call is `main.rs:880` — after the guard at `:396`.

**Is `ClipFlag` the right shape?** Yes, and I looked for a simpler alternative. A `[(&str, bool);
4]` const table would re-introduce the string-carried ownership Assumption 4 rules out and would not
give exhaustiveness. The enum is right. The one design change I would make is I1's: carry `ClipFlag`
inside `ClipOverride` and add `applies(paired)` alongside `read_2()`, which is what makes the enum a
single source rather than a fourth one, and is output-neutral.

**Would Validation catch what matters?** Mostly. V4/V5/V6 are well-chosen and each names a mutation
that can genuinely go red — V5 in particular catches the dangerous direction (a paired Read 2 flag
ceasing to refuse would silently lose annotations). The gaps are I3 (the ladder's order is measured,
not tested) and I4 (no mixed-ownership row). V7 is the right idea and its two refinements — temp files
plus `cmp` rather than `diff <(…) <(…)`, and cells with two or more flags set — are both correct and
both learned. One scope narrowing for V7: `flag_summary` is reached only when a preset is present
(`main.rs:1122` and `report.rs:275-276` both guard on `if let Some(preset)`), so site 1's differential
needs only the preset cells; the non-preset cells exercise site 2 alone.

**Efficiency.** Nothing to add. Once-per-invocation over a four-element `const` array, `Copy`, no
allocation. The claim is correct and the section is appropriately short.

---

## 6. Open verification items — could not settle in this tree

The working tree holds another session's committed work and cargo holds a build lock, so I ran no
build, test, or binary. These remain unverified by me:

1. **V1/V2 by execution** — that the defect reproduces at `4f7b907` and that the narrowed run produces
   a BAM identical to the no-clip-flag run. `assert_ubam_eq` does exist
   (`tests/integration_ubam_out.rs:108`) and `sp_trimmed.bam` is confirmed as the correct SE uBAM
   output name (asserted at `:1090`), so V2 is implementable as written.
2. **V4/V5/V6 mutations actually going red** — reasoned through, not executed.
3. **V8's exact counts** and a clean `clippy -D warnings`.
4. **V7's differential** — the only check for I3, and the reason I3 asks for a committed assertion
   instead.

---

## 7. Action items

**Critical**
1. C1 — rewrite step 13 as named edits to `docs/src/content/docs/guide/outputs.md:100` and `:145`;
   decide explicitly about the `--help` text at `src/cli.rs:312-313`; record the concept-level grep
   pattern in the plan.

**Important**
2. I1 — extend the read-ownership table to `ClipOverride.read_2` (`library.rs:93`, `:187-199`),
   `main.rs:1128` and `report.rs:287`; either fold them in via `ClipOverride { flag: ClipFlag, … }`
   plus `ClipFlag::applies(paired)` (output-neutral — `read_2` is not serialised), or scope them out
   as a stated Q-4.
3. I2 — narrow `any_effective`'s doc comment to read shape only; it does not know whether the mode
   honours clip flags.
4. I3 — pin the warning ladder's emission order in a committed test; today nothing does, while
   `flag_summary`'s order is pinned four ways.
5. I4 — add the mixed-ownership row (`--clip_R1 3 --clip_R2 5` single-end → still refused) to B1 and
   V3.
6. I5 — doc-scope `Cli::user_clips()` on the method itself; consider `pub(crate)`.

**Optional**
7. O1 — fix the four citations; name `tests/integration_inert_clip_flags.rs` in V7.
8. O2 — note in one clause that `any_set` is published public API and that removal is deliberate.
9. O3 — split `read_2()`'s doc from `any_effective`'s consequence.
10. Update Assumption 6: the parallel session has committed, the working tree is clean, and
    `src/parallel.rs` / `src/fastqc.rs` now also differ from `origin/dev`. The worktree conclusion
    stands.
