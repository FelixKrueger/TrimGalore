# PLAN_REVIEW_B — `r1_unpaired` / `r2_unpaired` cores parity (#544)

**Reviewer:** B (independent; no contact with Reviewer A)
**Plan reviewed:** `/Users/fkrueger/Github/TrimGalore/plans/11112026_unpaired-counter-cores-parity/PLAN.md`
**Tree:** `dev` @ `86aa301`, tracked files clean. No source edited; measurements taken with the
pre-built `./target/release/trim_galore` and read-only greps.

---

## Verdict up front

**The fix itself is correct and the call-graph argument holds.** I re-traced every link and every
claim about `classify_paired` → `process_pairs` → the two call sites, and the equivalence the plan
rests on — `Some(writer_up_r*) && r*_ok` ≡ the serial predicate at `src/trimmer.rs:584-595` — is
true. Assumption 5, which the plan flags as unverified, is **true**. §2's measurement table
reproduces exactly. Steps 1–3 are mechanical and I would ship them as written.

**The test plan is where this needs work, and two of the findings are blocking.** Step 4's fixture
shape and Validation 7 are mutually unsatisfiable; and Step 5's `--retain_unpaired` half cannot be
hosted where the plan suggests, because the CLI forbids that combination. Separately, the plan
missed the *other* existing test that asserts these very fields — which happens to be the right home
for Steps 5–7 and already implements the Step 6 and Step 7 patterns the plan proposes inventing.

---

## 1. Logic review

### 1.1 Verified in source — the plan's core argument stands

Everything in this subsection I read directly. All of it checks out.

| Plan claim | Verified at | Result |
|---|---|---|
| Unconditional increments in `classify_paired`'s `TooShort` arm | `src/parallel.rs:571-576` | ✅ exact |
| `classify_paired` has exactly one caller | call at `src/parallel.rs:628`; `fn process_pairs` at `:602` | ✅ one caller |
| `PairOutcome::Unpaired` has exactly one consumer | constructed `:577`, matched `:642` | ✅ |
| `writer_up_r1` / `writer_up_r2` are parameters | `src/parallel.rs:613-614` | ✅ exact |
| `Unpaired` arm already writes under `Some(w) && r*_ok` | `src/parallel.rs:649-661` | ✅ |
| `process_pairs` has two call sites | `src/parallel.rs:386` (gzip), `:415` (plain) | ✅ exact |
| Both pass `Some` exactly when `retain_unpaired` | `:370-379` (gzip), `:430-438` (plain) | ✅ exact |
| `retain_unpaired = unpaired_r1_path.is_some()` | `src/parallel.rs:108` | ✅ exact |
| `PairValidationStats::merge` is a field-wise sum | `src/report.rs:148-149` | ✅ |
| Serial increments only when the writer is open | `src/trimmer.rs:584-595` | ✅ exact |
| `gen_multiome_rows` makes both mates short | `src/parallel.rs:2283-2309`, doc `:2275-2282` | ✅ exact |

I also closed a route the plan does not mention but which could have broken the argument:
`process_paired_batch` (`src/parallel.rs:338`) is the only function that calls `process_pairs`, and
it has exactly one caller of its own (`src/parallel.rs:142`, the worker loop). `--clumpify` routes
through `run_paired_end_parallel` like any `--cores >= 2` run (`src/main.rs:1735`), so there is no
second batch processor with its own rescue logic. **Every parallel paired record reaches the counter
through `classify_paired`.** The fix therefore covers the whole parallel surface, including
`--clumpify`. Confidence: high.

I closed one further route for completeness: `run_paired_end_to_bam` (`src/trimmer.rs:708`) has its
own `TooShort` arm (`:805-809`) which neither rescues nor counts, with the comment
*"`--retain_unpaired` is rejected upstream — no rescue"* — correct, per `src/cli.rs:652-657`, and it
passes an `unpaired_sentinel` rather than real cutoffs. Out of scope and consistent.

### 1.2 Assumption 5 is TRUE — the plan's open risk #3 is closed

The plan quotes a source comment and asks the implementer to confirm. Confirmed:

- `src/cli.rs:948-953` — `--passthrough is incompatible with --retain_unpaired (passthrough requires
  strict pair semantics in v1)`.
- Regression test: `src/cli.rs:2021` `test_passthrough_rejects_retain_unpaired`.

So §3's edge-case row and §11's remaining risk #3 can both be marked resolved. **But this same fact
creates finding 2.1 below** — the plan treats the exclusion as a reason not to think about the
combination, when it is also a constraint on where the Step 5 test can live.

### 1.3 The serial path rescues only in the `TooShort` arm — scope is right

The brief asked whether the serial path rescues on other filter results. It does not.
`src/trimmer.rs:564-607`: `TooManyN` carries the explicit comment *"N-filter: always discard entire
pair, no rescue"*, and `TooLong` has no rescue block either. `src/filters.rs:82-115` returns
`TooManyN` before the length check and `TooLong` after it, and only `TooShort` carries `r1_ok` /
`r2_ok`. The parallel path mirrors this exactly (`src/parallel.rs:556-586`). **The plan's scope —
the `TooShort` arm only — is correct.** Confidence: high.

### 1.4 Count-after-write ordering: matches serial, and unobservable either way

The brief asked whether placing the increment after `r1.write_to(*w)?` matters. The serial path does
increment after the write (`src/trimmer.rs:586-588`), so the plan's Step 2 matches it. Whether it
matters: no. On either path a write error propagates `Err` out of the trimming entry point before any
report is produced — the parallel worker's error reaches the main thread via the result channel and
aborts the run. There is no path on which a report is written after a failed rescue write, so
count-before vs count-after is not observable. The plan's reasoning ("keeps the count false if the
write fails") is sound but describes a state nothing can read. **No change needed.** Confidence: high.

### 1.5 §2's measurement table — rows 1 and 3 reproduced

Re-derived with `./target/release/trim_galore` on `test_files/BS-seq_10K_R{1,2}.fastq.gz`
(outputs under `/private/tmp/claude-501/-Users-fkrueger-Github-TrimGalore/aecc6671-429c-462c-ad9c-38f8e4d35fbe/scratchpad/m1/`):

| Run | `--cores 1` | `--cores 4` |
|---|---|---|
| `--paired` (no clip, no rescue) | text block **absent**; JSON `r1_unpaired: 0, r2_unpaired: 0` | `Unpaired read 2 kept: 4`; JSON `r2_unpaired: 4` |
| `--paired --retain_unpaired` | `Unpaired read 2 kept: 4`; `*_unpaired_2.fq.gz` holds **4 records**, `*_unpaired_1.fq.gz` holds **0** | identical |

Both exits were 0 and both greps non-empty (checked explicitly, per the empty-input hazard). So:
**4 is the true rescue count**, the counters do equal the records written when rescue is on, and the
`--cores >= 2` no-rescue run reports a correct number for a question nobody asked. §2's three
conclusions are all supported. Confidence: high (measured).

One detail the table omits, relevant to Validation 6: the text block lands **only in the R2 trimming
report** (`BS-seq_10K_R2.fastq.gz_trimming_report.txt`). The R1 text report has no pair-validation
section at all, while the JSON carries `pair_validation` in **both** files. See finding 3.3.

### 1.6 BLOCKING — Step 4's fixture shape makes Validation 7 unsatisfiable

This is the one real logic defect.

`src/filters.rs:100-107` derives the flags as `r1_ok: r1.len() >= unpaired.r1`,
`r2_ok: r2.len() >= unpaired.r2`. Step 4 says to reshape `gen_multiome_rows` so that *"some records
have **only R1** short, which produces `r1_ok == false, r2_ok == true`"*. That is correct — and it is
the **only** shape Step 4 asks for.

Consequence: on that fixture `r1_ok` is `false` for every pair that reaches the `Unpaired` arm, and
pairs that do not reach the arm never touch either counter. So `r1_unpaired == 0` on both paths,
before and after the fix, **whatever the R1 increment does**. Validation 7 ("Revert the R1 increment
alone → a test fails") cannot be satisfied. The same hole exists in the manual BS-seq runs: §2's own
table shows `r1_unpaired` at 0 in every row.

Reverting the R1 increment alone would therefore pass the whole suite, which is exactly the class of
silent gap Validation 9 was added to close — reproduced on the other mate.

The fix is a fixture with **three** shapes rather than one, all of which `gen_multiome_rows` can
produce deterministically by index:

| Shape | Effect with `length_cutoff = 20`, `UnpairedLengths { r1: 35, r2: 35 }` | Covers |
|---|---|---|
| R1 16 bp, R2 48 bp | `r1_ok = false, r2_ok = true` | `r2_unpaired` → Validation 8 |
| R1 48 bp, R2 16 bp | `r1_ok = true, r2_ok = false` | `r1_unpaired` → **Validation 7** |
| R1 16 bp, R2 16 bp | both `false` | the both-fail row §3 says the current fixture pins |

Keeping the third shape matters independently: §3's edge-case table lists both-fail as *"the case the
current fixture pins"*, and Step 4 as written would un-pin it, trading one blind spot for another.
Confidence: high.

### 1.7 BLOCKING — Step 5's `--retain_unpaired` run cannot live in the passthrough parity test

Step 5 says to *"extend or add alongside `test_passthrough_serial_parallel_parity`"* so parity is
asserted twice, *"once with `--retain_unpaired`"*. §10 then calls extend-vs-sibling *"implementer's
choice"*.

For the `--retain_unpaired` half it is not a choice. `test_passthrough_serial_parallel_parity`
(`src/parallel.rs:2474`) passes `Some(open_fq(&pt_path))` / `Some(&mut rdr_pt)` on both paths, and
`src/cli.rs:948-953` rejects `--passthrough` with `--retain_unpaired`. The test calls the library
directly so it *can* be forced, but doing so configures a combination the product refuses — and the
`Unpaired` arm's own comment (`src/parallel.rs:643-648`) leans on the exclusion to call its
`passthrough_records_dropped` accounting *"unambiguous"*. A test that enables both would be pinning
behaviour that cannot occur and would silently become the only coverage of an illegal state.

Since Validation 9 — which the plan calls *"the one that distinguishes a fix from a cover-up"* —
depends entirely on that `--retain_unpaired` assertion existing, this is load-bearing, not
stylistic. **The plan must commit: the rescue-on half goes in a non-passthrough test.** Confidence:
high.

### 1.8 IMPORTANT — the plan missed the second blind test, which is the right host

§2's *"Why no test catches it"* names `test_passthrough_serial_parallel_parity` as *the* test
covering these fields. There is a second one:

`test_parallel_serial_bp_after_trim_parity_paired` (`src/parallel.rs:1622-1729`) — the test PR #543
added for the sibling issue — asserts full `PairValidationStats` parity at `:1705-1710`, with the
identical assertion message. It is blind to #544 for a *different* reason: `parity_test_config()`
sets `length_cutoff: 0` (`src/parallel.rs:1433`), inherited by
`parity_test_config_post_adapter()` (`:1456-1464`), so `src/filters.rs:103` can never return
`TooShort` and the `Unpaired` arm is unreachable. Both counters stay 0 on both paths; the assertion
is vacuous for these two fields.

So there are **two** existing tests asserting these fields and both are blind, for two unrelated
reasons. The plan diagnoses one. That matters because the one it missed is a markedly better host
for Steps 5–7:

- **Non-passthrough** (`None` for the pt reader/writer on both paths), so `--retain_unpaired` is
  legal there — it resolves finding 1.7 directly.
- It already has the serial-vs-parallel harness at `--cores 4`.
- It already implements **Step 6** verbatim: `const EXPECTED_R1 / EXPECTED_R2` with the comment
  *"Hand-derived: R1 is 15 × 30 bp trimmed + 15 × 34 bp clean"* (`:1712-1717`).
- It already implements **Step 7** verbatim: *"config must clip both reads, or this test cannot
  fail"* (`:1721-1726`) — the same discriminating-power guard the plan proposes writing from scratch.

A sibling in the same `mod` following the same shape, with `length_cutoff` raised and the three-shape
fixture from 1.6, satisfies Steps 5, 6 and 7 together and needs no passthrough machinery at all.
I would make that the plan's recommendation rather than leaving it to the implementer.
Confidence: high.

### 1.9 IMPORTANT — §3's edge-case table is missing the both-rescued case

The table's row *"`--retain_unpaired`, both mates pass | Pair is not `TooShort`; arm not reached"*
conflates two different things: passing the **discard** cutoff (`--length`) and passing the
**rescue** cutoff (`--length_1` / `--length_2`). `src/filters.rs:100-107` sets the two flags
independently of each other and of which mate failed — and `src/filters.rs:88-99` spells out why, as
the #245 parity fix:

> *"when at least one mate fails the discard cutoff, each read is independently routed to its unpaired
> file based solely on its own length vs the per-side cutoff. There is no requirement that the OTHER
> mate must have failed; **both sides can land in unpaired in the same record**."*

So the genuinely interesting fourth case is absent: `length_cutoff = 40`, both mates 37 bp → both
short vs the discard cutoff, both ≥ 35 → `r1_ok == r2_ok == true`, **both rescued, one pair
incrementing both counters**. That is the shape that would catch a fix using `else if` instead of two
independent `if`s, or one that gated R2's increment on `!r1_ok` — and it is the exact regression
#245 fixed, so it is a live failure mode in this code, not a hypothetical. Add the row, and ideally a
fourth fixture shape for it. Confidence: high.

### 1.10 IMPORTANT — "every consumer gates on `> 0`" is false for the JSON

§2 states *"every consumer gates on the counters being non-zero"* and §7 repeats *"All gate on
`> 0`"*. Three of the four do (`src/report.rs:469`, `src/main.rs:1902`, `:2123` — all verified). The
JSON does not: `src/report.rs:1073-1074` sits inside `match pair_stats { Some(ps) => … }`
(`:1060-1076`) and calls `json_int` unconditionally, so the keys are always present.

Measured confirmation: `--cores 1` emits `"r2_unpaired": 0`, `--cores 4` emits `"r2_unpaired": 4` —
same keys, different value, in **both** per-mate JSON files.

Two consequences. (a) §2's framing — *"a spurious count does not merely print a wrong number"* — is
right about the text report and stderr but wrong as stated about the whole surface; for JSON
consumers it **is** exactly a wrong number. (b) It changes what Step 8's CHANGELOG sentence should
name: the user-visible change for anyone parsing the JSON (MultiQC, nf-core) is `4 → 0` in a field
that was always there, not the appearance/disappearance of a block. The plan's Validation 1 already
expects the right thing, so this is a prose-accuracy defect, not a validation gap.
Confidence: high (measured).

### 1.11 IMPORTANT — wrong, self-contradictory line references for the parity test

§2 cites the test at `src/parallel.rs:2360` with its `PairValidationStats` assertion at
`:2338-2342` — an assertion 20 lines *before* the function that contains it, which cannot be right,
so the plan is internally inconsistent here. Actual:

- `#[test]` `:2473`, `fn test_passthrough_serial_parallel_parity` `:2474`
- `PairValidationStats` assertion `:2532-2537`
- `:2338-2342` is inside `write_multiome_fixture`, unrelated

(`gen_multiome_rows` at `:2283` and its doc comment at `:2275-2282` are both correct, so this is a
localised slip, not a systematic one.) Given the brief's note that the sibling plan shipped three
false factual claims in its first revision, correct these before an implementer navigates by them.
Confidence: high.

### 1.12 Step 4's blast radius is nil — and §10's hedge points away from the fix

The brief asked how many callers `gen_multiome_rows` has and whether any depend on the 33%
both-short drop rate. Five callers; I read every assertion in each:

| Caller | Test | Assertions that could be affected |
|---|---|---|
| `:2396` | `test_passthrough_serial_smoke` (n=50) | `checked == 50`; `kept + dropped == 50`; `dropped > 0`; `kept > 0`; lockstep; pt seq verbatim |
| `:2477` | `test_passthrough_serial_parallel_parity` (n=12 000) | serial == parallel; record identity; lockstep; `checked == 12 000`; `kept + dropped == 12 000` |
| `:2589` | `test_passthrough_truncation_detected_mid_stream` (n=10 000) | error string only |
| `:2633` | `test_passthrough_id_mismatch_detected` (n=200) | error string only |
| `:2680` | `test_passthrough_dont_gzip_plain_output` (n=300) | `checked == 300`; `r1_recs.len() == kept`; lockstep |

**None depends on both mates being short.** They depend only on *which pairs are dropped*, and
lengthening R2 in the short rows does not change that: `src/filters.rs:100-103` returns `TooShort`
when `r1_short || r2_short`, and R1 at 16 bp is already below `length_cutoff = 20`. The drop set
stays "every third record" under every reshaping in 1.6. No absolute basepair or length assertion
exists on this fixture. So mutating the shared generator is safe, and §10's worry is unfounded.

More importantly, §10's conclusion — *"the right move is a second generator for the asymmetric case
rather than mutating the shared one"* — points at the option that does **not** fix the problem §2
identifies, unless the variant is also wired into the parity test. An implementer who adds a variant
consumed only by a *new* test leaves `test_passthrough_serial_parallel_parity` exactly as blind as
§2 says it is. The hedge is not adequate as written: it should either state plainly that the shared
generator is safe to mutate (it is — I checked), or require that the variant replace
`gen_multiome_rows` **in the parity test**. Confidence: high.

Note also that the smoke test's own comment at `:2427` already says *"every 3rd row has short R1"*,
contradicting the generator's doc comment two functions away. Step 4's cleanup should take that line
too.

---

## 2. Assumptions

| # | Plan's assumption | Status |
|---|---|---|
| 1 | Counters mean "records written" | **Verified, no conflict found.** `src/cli.rs:277-280` says *"Retain unpaired reads when the mate is too short"*; `grep` over `docs/`, `wiki/`, `README.md` for `Unpaired read` / `r1_unpaired` returns **nothing**, so no documentation defines the field at all and none contradicts the report wording. Measured: rescue-on counts equal the records in the files. Treating it as a fixed rule is right. |
| 2 | `writer_up_r*` is `Some` exactly when `--retain_unpaired` | **Verified** at `:108`, `:370-379`, `:430-438`. See 3.1 for a latent narrowing. |
| 3 | One caller, one consumer | **Verified**, plus the `process_paired_batch` route closed (§1.1). |
| 4 | `merge` is a field-wise sum | **Verified** `src/report.rs:148-149`. |
| 5 | `--passthrough` + `--retain_unpaired` rejected | **Verified true** — `src/cli.rs:948-953`, test `:2021`. Mark the §11 risk closed. But see 1.7: it is also a constraint, not only a relief. |
| 6 | Changing `gen_multiome_rows` affects other tests | **False, and safely so** — see 1.12. Nothing else asserts on the shape. |

**Unstated assumptions I surfaced:**

- *That the three-shape fixture is reachable under the existing configs.* It is not: both
  `parity_test_config()` (`length_cutoff: 0`) and `passthrough_test_config()` (`length_cutoff: 20`,
  with `UnpairedLengths { r1: 35, r2: 35 }` supplied at the call sites) need checking against
  whichever shapes the implementer picks. The plan never states the cutoff values the new assertions
  depend on, which is the main thing an implementer would have to invent (see §6).
- *That `r1_ok` and `r2_ok` are mutually exclusive-ish.* They are fully independent — §1.9.
- *That the serial path is the reference.* Correct here, but worth saying once: the fix makes the
  parallel path match the serial path, and the serial path matches Perl 0.6.11
  (`src/filters.rs:88-99` cites `0.6.11:trim_galore:2325-2343`). That strengthens §1's framing and is
  worth one sentence in §7, since the CI matrix cannot demonstrate it.

---

## 3. Lower-severity observations

### 3.1 `retain_unpaired` is derived from the R1 path alone

`src/parallel.rs:108` derives one flag from `unpaired_r1_path.is_some()` and uses it to gate **both**
writers, where the serial path gates R1 and R2 independently on their own `Option`s
(`src/trimmer.rs:583`, `:591`). Both callers of `run_paired_end_parallel` pass both-`Some` or
both-`None` (`src/main.rs:1706-1714`, `:2021-2029`), so this is unreachable today — but after the fix
the *counters* inherit the coupling too, so a future asymmetric caller would make them diverge again.
A `debug_assert!(unpaired_r1_path.is_some() == unpaired_r2_path.is_some())` beside `:108`, matching
the existing passthrough `debug_assert!` at `:109-113`, pins it for one line. Optional; confidence
high that it is currently unreachable.

### 3.2 Validation 10's "846+ pass"

`grep -c '#\[test\]'` over `src/*.rs` + `tests/*.rs` gives **845** today, so "846+" is consistent with
adding at least one test. Not a finding; noting it because the brief asked for the table's numbers to
be checked rather than trusted.

### 3.3 Validation 6 should name both reports

"diff whole text reports before/after for the rescue run" — the pair-validation section exists only
in the **R2** report (measured, §1.5). Diffing only one file catches half the surface; diffing the R1
report alone catches none of it. Spell out both filenames, and add the two JSON files, since the JSON
is where the value (not the block) changes.

### 3.4 Step 6's anchor must be derived, not captured

Step 6 says to *"assert the absolute expected count on a fixture with hand-derived short-mate
numbers"*. Worth one explicit sentence that the number is computed from the generator's index rule
(e.g. "n/3 rows are R1-short-only ⇒ `r2_unpaired == n/3`"), **not** read off a first run after the
fix — otherwise the anchor pins whatever the code produced and Step 6 degrades into the
`serial == parallel` tautology it exists to escape. `:1712-1717` in the sibling test is the model to
copy, down to the arithmetic in the comment.

### 3.5 Degenerate "fixes" the validation table would and would not catch

The brief asked whether Validation 9 is correctly specified and what else slips through. Validation 9
as stated ("revert both increments; the `--retain_unpaired` assertion must fail") is correctly
specified **provided** that assertion exists outside the passthrough test (1.7). Sweeping the
remaining degenerate shapes:

| Degenerate fix | Caught? |
|---|---|
| No change at all (increment on `r*_ok` only, no writer test) | ✅ Validation 1 — no-rescue run reports 0 |
| Both counters hard-zeroed | ✅ Validation 9, given 1.7 is resolved |
| Gate on `writer.is_some()` but drop the `r1_ok` test | ✅ on the R1 side — over-counts vs the file, fails Validation 3/4 |
| Gate on `writer.is_some()` but drop the `r2_ok` test | ❌ **not caught** on an all-dropped-pairs-are-R1-short fixture: `r2_ok` is `true` for every pair that reaches the arm, so removing the test changes nothing. Needs a pair where the writer is open and `r2_ok == false` — i.e. the both-short shape from 1.6 |
| Count pairs rather than mates (`+= 1` to both whenever either wrote) | ✅ on BS-seq, where R1 and R2 counts differ (0 vs 4) — Validation 3/4 |
| `else if` instead of two independent `if`s | ❌ **not caught** without the both-rescued shape from 1.9 |

Two of the six slip through on the fixture Step 4 specifies. Both are closed by the three-or-four
shape fixture; that is the single highest-value change to this plan.

---

## 4. Efficiency analysis

Nothing to report, and I agree with §6. Two `usize` increments move inside guards that are already
being evaluated for the writes a line away; no allocation, no new branch, no change to the per-record
critical path, nothing added to the merge. Per-worker `PairValidationStats` already exist
(`src/parallel.rs:350`) and already merge field-wise. Complexity, memory and scalability are all
unchanged. The only cost in this plan is test runtime, and Step 4 reuses fixtures that already run.
Confidence: high.

---

## 5. Validation sufficiency

**Sufficient once 1.6, 1.7 and 3.5 are addressed; not before.** The table covers the right axes —
core-invariance, counts-match-files, record identity, per-site mutation, the anti-cover-up check —
and Validation 9 is a genuinely good addition that most plans would omit. The gaps are all in the
fixture the validations run against, not in the validations themselves:

- **Validation 7 is unsatisfiable** on the Step 4 fixture (1.6). This is the blocking one.
- **Validation 9 has no legal home** as currently placed (1.7).
- **Two degenerate fixes pass everything** (3.5).
- **Validation 6 under-scopes the surface** (3.3).

Where could the code silently produce a wrong result with nothing catching it, after this plan ships
as written? Precisely one place: a future change to the `Unpaired` arm that drops the `r2_ok`
conjunct, or merges the two guards. The suite would stay green, the no-rescue report would stay
clean, and `--retain_unpaired` runs would over-count `r2_unpaired` relative to the file — the same
class of defect as #544, one conjunct over. The three-shape fixture closes it.

One thing the plan is right to say and I will underline: the CI validation matrix cannot help here.
It runs `--cores 1` (§7), which is the already-correct path. Steps 5–7 are the entire guard, which is
why their under-specification is the most consequential part of this review.

---

## 6. Under-specified for an implementer

Three places where an implementer would have to invent behaviour:

1. **Which test gets the `--retain_unpaired` assertions** (1.7, 1.8). Currently "implementer's
   choice" between two options, one of which is illegal. Name the target.
2. **The fixture's exact shapes and the cutoffs they depend on** (1.6, 1.9). Step 4 gives one shape
   and no `length_cutoff` / `UnpairedLengths` values; the counters only move for specific
   combinations of the two, and `parity_test_config()`'s `length_cutoff: 0` silently defeats all of
   them. Give the table.
3. **Whether to mutate `gen_multiome_rows` or add a variant** (1.12). §10 states a preference that
   does not resolve §2's problem. Decide it: the shared generator is safe to mutate.

Everything else — Steps 1, 2, 3, 8 — is specified tightly enough to code from without guessing.

---

## 7. Alternatives

**(a) Pass the writer state into `classify_paired`.** Add a `retain_unpaired: bool` parameter and
keep the increments where they are. Rejected, correctly by implication: it adds a parameter to a
function that already takes seven, duplicates a predicate that the caller must evaluate anyway, and
leaves two sources of truth about whether a rescue happened. The plan's choice is better.

**(b) Derive the counters from the writers instead of counting.** Not available — the per-worker
sinks are `Vec<u8>` / `GzEncoder` byte buffers with no record count.

**(c) Change the report's semantics instead of the counter** — redefine `Unpaired read N kept` as
"would have been kept". Worth one line in §10 to dismiss rather than passing over: it would make the
parallel path's current output correct by definition, requires no code change, and is wrong, because
the serial path and Perl 0.6.11 both mean "written", the field name is `*_kept`, and a run with no
unpaired file reporting a non-zero "kept" is incoherent whatever the gloss. Mentioning and rejecting
it makes §8's assumption 1 a decision rather than an inheritance.

**(d) Gate the report block on `--retain_unpaired` rather than on `> 0`.** Treats the symptom in the
reporter, leaves the JSON wrong (1.10) and the struct lying. Not competitive.

**(e) Assert on report *output* rather than on the struct.** The plan's Steps 5–7 all assert on
`PairValidationStats`. An end-to-end check that the text block is absent — essentially Validation 1
as a `tests/` integration test — would additionally pin the `> 0` gating and the
only-in-the-R2-report placement (1.5). Not necessary, but it is the cheapest way to make Validation 1
permanent rather than a one-off manual step; the plan currently leaves Validations 1–6 as manual
checks that nothing re-runs. Worth considering for one of them.

---

## 8. Action items

### Critical

1. **Give the fixture three shapes, not one** (§Step 4; `src/parallel.rs:2283`; `src/filters.rs:100-107`).
   R1-short-only, R2-short-only, and both-short. As written, Step 4 produces only `r1_ok == false`,
   so `r1_unpaired` stays 0 on both paths and **Validation 7 cannot be satisfied** — reverting the R1
   increment alone would pass the suite. Keeping both-short also preserves the edge case §3 says the
   current fixture pins. *Confidence: high; verified in source.*

2. **Commit the `--retain_unpaired` assertions to a non-passthrough test** (§Step 5, §10;
   `src/cli.rs:948-953`, test `:2021`; `src/parallel.rs:643-648`). `--passthrough` with
   `--retain_unpaired` is rejected at the CLI, so `test_passthrough_serial_parallel_parity` cannot
   legally host that half, and §10's "implementer's choice" offers an option that is not available.
   Validation 9 depends on this assertion existing. *Confidence: high; verified in source.*

### Important

3. **Point Steps 5–7 at `test_parallel_serial_bp_after_trim_parity_paired`** (`src/parallel.rs:1622`,
   assertion `:1705-1710`). §2 misses that this second test also asserts
   `PairValidationStats` parity and is likewise blind — not from fixture shaping but because
   `parity_test_config()` sets `length_cutoff: 0` (`:1433`), making the `TooShort` arm unreachable.
   It is non-passthrough (resolving item 2), and already implements Step 6 (`EXPECTED_R1`/`EXPECTED_R2`,
   `:1712-1717`) and Step 7 ("or this test cannot fail", `:1721-1726`). A sibling in its shape
   satisfies Steps 5–7 at once. *Confidence: high.*

4. **Add the both-rescued edge case to §3** (`src/filters.rs:88-107`). Both mates below `--length`
   but both ≥ `--length_1/2` → `r1_ok == r2_ok == true`, both rescued from one pair. The table's
   "both mates pass → arm not reached" row conflates the discard cutoff with the rescue cutoff. This
   is the #245 regression shape, so it is a live failure mode, and it is one of the two degenerate
   fixes that would otherwise slip through (§3.5). *Confidence: high.*

5. **Correct "every consumer gates on `> 0`"** (§2, §7; `src/report.rs:1060-1076`). The JSON emits
   `r1_unpaired` / `r2_unpaired` unconditionally — measured as `4` at `--cores 4` against `0` at
   `--cores 1`, same keys. For JSON consumers the fix changes a *value*, which is what Step 8's
   CHANGELOG sentence should name. *Confidence: high; measured.*

6. **Fix the parity-test line references** (§2). `:2360` / `:2338-2342` are wrong and mutually
   impossible; actual `:2474` with the assertion at `:2532-2537`. *Confidence: high.*

7. **Resolve §10's Step-4 hedge in favour of mutating the shared generator.** I checked all five
   callers (`:2396`, `:2477`, `:2589`, `:2633`, `:2680`) and no assertion depends on both mates being
   short — only on which pairs are dropped, which no reshaping changes because R1 at 16 bp is already
   below the cutoff. As written, the hedge points toward a variant that, if consumed only by a new
   test, leaves the parity test as blind as §2 says it is. *Confidence: high.*

8. **Close assumption 5 and §11 risk 3 as verified** — `src/cli.rs:948-953` plus the test at `:2021`.
   *Confidence: high.*

### Optional

9. Spell out in Steps 5–7 the `length_cutoff` and `UnpairedLengths` values each new assertion
   depends on — `parity_test_config()`'s `length_cutoff: 0` silently defeats every shape in item 1.
10. Scope Validation 6 to both text reports and both JSON files; the text block lands only in the R2
    report (measured).
11. Say in Step 6 that the anchor is derived from the generator's index rule, not captured from a
    post-fix run.
12. Note in §8 that the serial path's semantics are the Perl 0.6.11 semantics
    (`src/filters.rs:88-99` cites `0.6.11:trim_galore:2325-2343`) — it turns assumption 1 from an
    inheritance into a decision, which matters because the CI matrix cannot demonstrate it.
13. Add `debug_assert!(unpaired_r1_path.is_some() == unpaired_r2_path.is_some())` at
    `src/parallel.rs:108`, matching the passthrough assert at `:109-113`. Unreachable today
    (`src/main.rs:1706-1714`, `:2021-2029`), but the counters now inherit the coupling.
14. Mention and dismiss the "redefine the counter as *would have been* kept" alternative in §10
    (§7(c) above), so assumption 1 reads as chosen rather than assumed.
15. Fix the stale comment at `src/parallel.rs:2427` ("every 3rd row has short R1") while Step 4 is in
    the generator.

### Confirmed sound — no action

- Steps 1–3 and the whole call-graph argument (§1.1).
- Scope limited to the `TooShort` arm (§1.3).
- Count-after-write ordering (§1.4) — matches serial, and unobservable either way.
- §2's measurement table and its three conclusions (§1.5) — reproduced.
- §4 (no signature change), §6 (no efficiency impact) — both correct.
- Step 8's one-sentence CHANGELOG, mirrored into
  `docs/src/content/docs/reference/changelog.md` (both files verified present, both carrying the
  `### Unreleased` / `#### Bug fixes` #541 entry in the format to copy). One sentence is the right
  size per the house rule; item 5 only changes *what* it should say.

---

## Note on the CHANGELOG question

The brief asked whether a CHANGELOG sentence is sufficient for removing currently-emitted report
lines. I think yes, and that the house rule settles it: one sentence, with a `**Breaking:**` prefix
only if you consider this breaking. I would not — the removed lines were wrong, no documented
interface promised them, and the change brings `--cores >= 2` into line with `--cores 1` rather than
away from anything. The one thing worth getting right is *what* the sentence names: per item 5 the
load-bearing change for downstream tooling is the JSON value going `4 → 0` on no-rescue parallel
runs, not the text block disappearing. A reader deciding whether an upgrade affects them is a
MultiQC/nf-core consumer, and that is the sentence they need.

---

**Report path:** `/Users/fkrueger/Github/TrimGalore/plans/11112026_unpaired-counter-cores-parity/PLAN_REVIEW_B.md`
