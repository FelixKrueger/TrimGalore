# PLAN_REVIEW_A — `r1_unpaired` / `r2_unpaired` cores parity (#544)

**Reviewer:** A (independent; no contact with Reviewer B)
**Plan under review:** `plans/11112026_unpaired-counter-cores-parity/PLAN.md`
**Tree:** `dev` @ `86aa301`, tracked files clean
**Date:** 2026-11-11

**Verdict in one line:** the fix is correct and the call-graph argument holds — I re-verified every
load-bearing claim in source, plus two of the three measurement rows with the shipped binary. But the
**test plan does not do what the plan says it does**: as specified, Step 4 exercises only the R2
increment site, so Validation 7 cannot fail, and the plan has missed a second `PairValidationStats`
parity test (and a third shaped fixture) that change where Steps 4–7 should land.

---

## 0. What I verified, and how

Everything below marked ✅ I read in source at `86aa301` or reproduced with
`./target/release/trim_galore`. Nothing in this report is inferred from the plan's own text.

### The core argument — holds

| Plan claim | Status |
|---|---|
| `classify_paired` has exactly one caller | ✅ `src/parallel.rs:628`, the only call; grep shows no other. (The plan cites `:602`, which is `fn process_pairs`'s *definition* line, not the call. Defensible as "the caller is `process_pairs`", but the call site is `:628`.) |
| `PairOutcome::Unpaired` has one consumer | ✅ constructed at `:577`, matched only at `:642` |
| `process_pairs` has two call sites | ✅ `:386` (gzip) and `:415` (plain) |
| Both pass the unpaired writers `Some` iff `retain_unpaired` | ✅ `:370-379` (GzEncoders), `:430-438` (`Vec<u8>` sinks) |
| `retain_unpaired = unpaired_r1_path.is_some()` | ✅ `src/parallel.rs:108` |
| The caller's chain is `Some` iff `--retain_unpaired` | ✅ `src/main.rs:1706-1713` builds both paths only under `if cli.retain_unpaired`, else `(None, None)` |
| `Some(writer) && r*_ok` ≡ the serial predicate | ✅ serial is `if let Some(ref mut w) = unpaired_r1.as_deref_mut() && r1_ok` at `src/trimmer.rs:584-589`; parallel's arm is `if let Some(ref mut w) = writer_up_r1 && r1_ok` at `src/parallel.rs:649-652`. **The equivalence the whole fix rests on is real.** |
| `PairValidationStats::merge` is a field-wise sum | ✅ `src/report.rs:148-149` |
| Serial increments are only the two sites | ✅ `trimmer.rs:588`, `:594` |
| Parallel increments are only the two sites | ✅ `parallel.rs:572`, `:575` |
| `TooLong` / `TooManyN` never rescue on either path | ✅ `trimmer.rs:564-573` and `:600-607`; `parallel.rs:560-566` and `:579-585`. **The plan's scope is right** — rescue is a `TooShort`-only concept on both paths, so there is no fourth increment site hiding in another filter arm. |
| (not in the plan) The third paired path has no increment site | ✅ `trimmer::run_paired_end_to_bam` matches `PairFilterResult::TooShort { .. }` at `src/trimmer.rs:805` and discards — it never rescues and never increments, consistent with `--retain_unpaired` being rejected for `--output-format ubam` (`src/cli.rs:652-657`). So the two serial sites really are the only ones. |
| The `Unpaired` arm's comment is at `:643-648`, writes at `:649-658` | ✅ exact |
| `gen_multiome_rows` doc comment admits the shaping | ✅ `src/parallel.rs:2271-2282`; the words "so neither path increments … and the parity test passes" are there verbatim |

### Measurements — re-derived

I ran the shipped binary (current with `86aa301`) rather than trusting the table.

Row 1 (`--paired`, no rescue), `test_files/BS-seq_10K_R{1,2}.fastq.gz`:

- `--cores 1` → no `Unpaired read` lines in either report; JSON `"r1_unpaired": 0, "r2_unpaired": 0` ✅
- `--cores 4` → `Unpaired read 1 kept: 0` / `Unpaired read 2 kept: 4`; JSON `"r2_unpaired": 4`; **no
  unpaired files on disk** ✅

Row 3 (`--paired --retain_unpaired`):

- `--cores 1` and `--cores 4` both report `read 2 kept: 4` ✅
- Both write `*_unpaired_1.fq.gz` (empty) and `*_unpaired_2.fq.gz` containing **16 lines = 4
  records** ✅ — so "4 is the correct rescue count" and "the figure equals the records written" are
  both confirmed, not assumed.

The plan's conclusion that the parallel path reports *a true number for a question nobody asked* is
exactly right.

### Assumption 5 — **verified, and it is true**

The plan flags this as quoted-not-checked. It checks out: `src/cli.rs:948-952` bails with
`"--passthrough is incompatible with --retain_unpaired (passthrough requires strict pair semantics in
v1)"`, inside the `if let Some(ref pt) = self.passthrough` envelope at `:934`. There is a unit test,
`test_passthrough_rejects_retain_unpaired`, at `src/cli.rs:2021`. The `Unpaired` arm's reliance on
this is sound.

(It has a consequence the plan did not draw — see **I2**.)

### Assumption 1 — true, and the plan under-sells its evidence

The plan rests "records written" on the v2 report wording `Unpaired read N kept`. The upstream
authority is stronger. In Perl 0.6.11 (`git show 0.6.11:trim_galore`):

- `++$read_1_printed` sits immediately after `print UNPAIRED1 …`, **inside `if ($retain)`**
  (lines 2322–2342) — structurally identical to `src/trimmer.rs:584-595`.
- The report line is `"Number of unpaired read 1 reads printed: $read_1_printed"` (line 2396).

So the counter means "records printed to the unpaired file", by upstream definition, and the serial
Rust path is a faithful port. There is no conflicting definition in `--help` or `docs/` — I grepped
`docs/src/content/docs/` and `README.md` for `Unpaired read` and found nothing, so the report text
and Perl are the only definitions in play. **Assumption 1 should be promoted from "what the report
claims" to "what Perl 0.6.11 counts", and the CHANGELOG can say the fix restores Perl behaviour** —
a materially better justification for removing output than "intended".

---

## 1. Logic review

### C1 (Critical) — Step 4 as specified makes Validation 7 unsatisfiable

**§5 Step 4 vs §9 row 7.**

Step 4 says: *"Change it so some records have **only R1 short**, which produces `r1_ok == false,
r2_ok == true` and exercises the counter."* Singular. That shape rescues **R2 only**. With it:

- `r2_unpaired` becomes non-zero → the R2 increment site is covered.
- `r1_unpaired` stays **0 in both paths** → the R1 increment site is covered by nothing.

Therefore Validation 7 — *"Revert the R1 increment alone → a test fails"* — **will not fail**.
Step 6's absolute anchor for R1 would be `== 0`, which a missing increment also satisfies. Step 7's
guard, `r1_unpaired + r2_unpaired > 0`, passes on the R2 contribution alone and does not notice.

This is the exact failure mode in `feedback_harness_false_negatives` #8: an assertion that can only
ever pass. And it is not hypothetical — the sibling test `test_parallel_serial_bp_after_trim_parity_paired`
exists in its current form precisely because #541 learned that `stats_r1` and `stats_r2` accumulate
at separate sites and *"a fix applied to only one must not pass"* (`src/parallel.rs:1619-1621`). The
same reasoning applies here and the plan has not carried it across.

**Fix:** the fixture must produce **both** asymmetries — some rows with R1 short and R2 long
(`r2_ok`), and some rows with R2 short and R1 long (`r1_ok`) — with *different* counts for the two,
so a copy-paste transposition of the two anchors also fails. E.g. with `length_cutoff = 20`,
`UnpairedLengths { r1: 35, r2: 35 }`: `i % 5 == 0` → R1 16 bp / R2 48 bp; `i % 5 == 1` → R1 42 bp /
R2 16 bp; and keep a both-short class if the existing drop-rate assertions need one. The plan should
commit to this shape rather than leaving it to the implementer.

Confidence: **high** — read directly off Step 4's text against Validation 7's text; no execution
needed.

### I1 (Important) — "the one parity test that covers these fields" is false; there are two, and the other one is the better host

**§1 and §2 ("Why no test catches it").**

`grep "PairValidationStats must match"` returns **two** sites:

1. `src/parallel.rs:2534` — `test_passthrough_serial_parallel_parity` (fn at `:2474`), the one the
   plan knows about.
2. `src/parallel.rs:1707` — `test_parallel_serial_bp_after_trim_parity_paired` (fn at `:1622`),
   which the plan never mentions. It asserts `serial.2 == parallel.2` at `:1704-1709` — the same
   struct, the same invariant.

Why does (2) pass today? Not fixture shaping: `parity_test_config_post_adapter()`
(`:1455-1463`) inherits `length_cutoff: 0` from `parity_test_config()` (`:1424`), so
`filter_paired_end` never returns `TooShort`, the `Unpaired` arm is never reached, and every rescue
counter is structurally dead. Its hand-derived anchors cover `total_bp_after_trim` only.

This matters for three reasons:

- The plan's framing — *"the one parity test that covers these fields has a fixture deliberately
  shaped so they stay zero"* — is wrong on the count, and wrong on the mechanism for the test it
  missed. Given the lead's note that #541's first revision carried three false claims, this is the
  same class of defect recurring.
- **Test (2) is the better home for Steps 5–7.** It is non-passthrough, so it has no
  `--retain_unpaired` incompatibility (see I2); it already hand-derives absolute expectations, which
  is the #541 pattern Step 6 says it wants to copy; it already carries a *"config must clip both
  reads, or this test cannot fail"* guard, which is the Step 7 pattern; and it is paired, serial-vs-
  parallel, at `cores = 4`. Adding a `length_cutoff` + asymmetric-mate variant alongside it is a
  smaller, cleaner change than reworking the passthrough fixture.
- A reader of test (2) today would reasonably believe `PairValidationStats` parity is covered. It
  should get a one-line note that the rescue counters are exercised elsewhere, or it will drift back
  into the same false confidence.

Confidence: **high** — both assertions read in source; `length_cutoff: 0` confirmed in both config
helpers.

### I2 (Important) — Step 5's `--retain_unpaired` run must not have passthrough active; the plan leaves this to taste

**§5 Step 5 vs §3 edge-case table vs `src/cli.rs:948`.**

Step 5 says to extend `test_passthrough_serial_parallel_parity` so parity is asserted twice, *"once
with `--retain_unpaired`"*. But §3's own edge-case table says `--passthrough` + `--retain_unpaired`
is *"rejected by `Cli::validate`, so unreachable"* — and I confirmed it is. The library entry points
take the passthrough reader and the unpaired writers as independent parameters, so a unit test
*can* pass both; it would then pin behaviour for a configuration the CLI forbids, and drive
`process_pairs` down the `pt.is_some() && writer_up_*.is_some()` path that the arm's own comment at
`:644-648` declares illegal. The serial side's comment (`src/trimmer.rs:579-583`) states the same
invariant and uses it to call its drop accounting "unambiguous".

The plan does note in §10 that *"a sibling avoids coupling these counters to passthrough"* — but
frames it as implementer's choice and an aesthetic preference. It is a **constraint**: the
`--retain_unpaired` run must pass `reader_passthrough = None`, or be a non-passthrough test
entirely (which I1 argues for anyway). As written, an implementer following Step 5 literally will
write a test that encodes a rejected combination, and the next maintainer will read it as evidence
the combination is supported.

Confidence: **high**.

### I3 (Important) — "every consumer gates on `> 0`" is wrong for the JSON

**§2 table and §7.**

§2 says *"every consumer gates on the counters being non-zero"* and §7 says *"All gate on `> 0`"*.
Verified consumer by consumer:

| Consumer | Gated? |
|---|---|
| Text report block, `src/report.rs:469-481` | ✅ gated on `> 0` |
| stderr summary, `src/main.rs:1902-1910` | ✅ gated |
| stderr summary (uBAM), `src/main.rs:2123-2130` | ✅ gated |
| **JSON, `src/report.rs:1073-1074`** | ❌ **not gated** — the `pair_validation` object is emitted unconditionally in PE mode (`:1060-1080`), so both fields are always present |

So the fix does two different things: it *removes lines* from the text report and stderr, and it
*changes a number* in the JSON (`"r2_unpaired": 4` → `0` on my BS-seq run). The JSON is the
machine-read surface — MultiQC-adjacent tooling and anything diffing reports across versions — and
is arguably the more consequential half. §2's "a spurious count does not merely print a wrong number"
framing is half right and should be split.

Two consequences:
- The CHANGELOG sentence should be written to cover both surfaces, not just the vanished block.
- Validation 1 already expects JSON `0`, which is correct; Validation 6 ("diff whole text reports")
  should say **the R2 report** — the paired block lands only there (confirmed on my runs, and
  `tests/integration_passthrough.rs:176` reads `r2.fq_trimming_report.txt` for the same reason).

Confidence: **high** — read at `src/report.rs:1060-1080`.

### I4 (Important) — nothing automated covers the user-visible symptom, and a third shaped fixture is unaccounted for

**§5 Steps 5–7, §9 rows 1–4.**

Every test the plan proposes calls `run_paired_end` / `run_paired_end_parallel` directly and asserts
on returned structs. Validations 1–4, which assert on the *report text and the files on disk*, are
manual CLI runs — they verify the fix once and guard nothing afterwards. The defect in #544 is a
**report-output** defect. A fix that counts correctly but a future change that alters the `> 0`
gating (or re-routes which report the block lands in) would reproduce the symptom with every
proposed test still green.

The repo has 18 binary-driving integration tests in `tests/`. One of them —
`tests/integration_passthrough.rs` — contains a **third copy of the shaped fixture**, with the same
admission in its comment (`:50-54`): *"R2 also short on dropped rows (avoids the pre-existing
parallel/serial `r2_unpaired` accounting divergence noted in the lib parity test)"*. The plan does
not mention it. Two things follow:

- That comment becomes a stale reference to a fixed defect the moment this lands. Per the project's
  own comment rules ("current state only — never reference the state you are fixing"), the
  justification clause has to go regardless of whether the fixture changes.
- That test runs only at `--cores 1` (`:106`, `passthrough_cli_cores_1_smoke`), so de-shaping it
  alone does not cover the parallel symptom. But it is the obvious skeleton for the missing guard: a
  `--cores 4` paired run, no `--retain_unpaired`, asserting the trimming report contains no
  `Unpaired read` line and the JSON reads `0`. That is ~30 lines and it is the only thing that pins
  the actual user-visible behaviour.

Confidence: **high** on the stale comment and the cores-1 restriction; **medium** on how much the
integration test is worth relative to its cost — but given that CI's Perl matrix runs `--cores 1`
(the plan correctly notes this), the unit tests are the *only* guard, and they guard the counter
rather than the output.

### O1 (Optional) — two wrong line citations, one of them internally impossible

**§2 "Why no test catches it".** The plan cites `test_passthrough_serial_parallel_parity` at
`src/parallel.rs:2360` with its `PairValidationStats` assertion at `:2338-2342`. Actual: the fn is at
`:2474`, the assertion at `:2531-2538`. Note `:2338-2342` is *before* the claimed fn line, which
cannot be true of an assertion inside it — so the pair is self-inconsistent on its face. Everything
else in the plan's citations I checked was accurate (`:571-576`, `:602`, `:642-662`, `:643-648`,
`:649-658`, `:2271-2282`, `:2283`, `:386`, `:415`, `:108`, `trimmer.rs:584-595`, `report.rs:148-149`,
`:469-479`, `main.rs:1902-1909`, `:2123-2130`). Low impact — but given #541's plan shipped three
false claims, worth correcting before an implementer navigates by them.

### O2 (Optional) — Step 2's ordering rationale is vacuous on this path

**§5 Step 2.** *"placing the increment after the write mirrors `src/trimmer.rs:584-595` and keeps the
count false if the write fails."* The first half is right and is reason enough. The second half is
not reachable here: on the parallel path `W` is `&mut Vec<u8>` or `GzEncoder<&mut Vec<u8>>`
(`:370-379`, `:430-438`), both of which write into memory and cannot fail, so the `?` never fires and
count-before vs count-after is unobservable. (On the serial path `w.write_record` hits a real file, so
there the argument does hold.) Keep the ordering for symmetry; drop the failure justification — and
in particular keep it out of the source comment, which Step 3 is already trimming for the same reason.

### O3 (Optional) — Perl gates the report line on `$retain`, not on `> 0`

Adjacent divergence, found while verifying Assumption 1, **out of scope for this fix**: Perl 0.6.11
prints both unpaired lines whenever `--retain_unpaired` is set (`if ($retain)`, line 2394), including
as `0`. Rust gates on `> 0` instead. So after this fix, a `--retain_unpaired` run that rescues
nothing prints nothing where Perl prints two zero lines. The *more Perl-faithful* shape would be
"gate the block on `retain_unpaired`, count on write" rather than "gate on `> 0`, count on write".
That is a bigger change and a separate decision — worth a follow-up issue, not scope creep here. Flag
it so the CHANGELOG does not over-claim Perl parity.

### O4 (Optional) — the affected surface includes every `--clumpify` run

**§3 edge-case table, last row.** *"`--cores 1` → worker pool bypassed; serial path already correct"*
is true but incomplete as a statement of blast radius. `src/main.rs:1736` dispatches on
`if cli.cores > 1 || cli.clumpify`, and only `--clump_only` rejects `--retain_unpaired`
(`src/cli.rs:921-925`) — plain `--clumpify --paired` is legal and always routes to the parallel path.
So a user who never typed `--cores` still sees the spurious lines under `--clumpify`. #541's
CHANGELOG entry named `--clumpify` explicitly for exactly this reason; this one should too.

### O5 (Optional) — Step 3's rewritten comment drops a live fact

**§5 Step 3.** The existing comment at `:643-648` carries two facts: (a) the rescue writes are no-ops
without `--retain_unpaired`, and (b) validation rejects `--passthrough` + `--retain_unpaired`, which
is what makes the `passthrough_records_dropped += 1` three lines below unambiguous. Replacing the
whole block with one line about counting deletes (b), which is still load-bearing for the code it
sits above. Suggest two lines, e.g. the proposed invariant line plus *"`--retain_unpaired` is
rejected with `--passthrough`, so a passthrough record here is always a drop."* The serial path keeps
its version of (b) at `trimmer.rs:579-583`; dropping it on one side only is how the two comments
start to disagree.

---

## 2. Assumptions

| # | Plan's assumption | Verdict |
|---|---|---|
| 1 | Counters mean "records written" | ✅ **True, and stronger than claimed** — Perl 0.6.11 counts after the print, inside `if ($retain)`, and calls it "printed". No conflicting definition in `--help` or `docs/`. Promote to a cited fixed rule. |
| 2 | `writer_up_r*` is `Some` iff `--retain_unpaired` | ✅ Verified end-to-end: `main.rs:1706-1713` → `parallel.rs:108` → `:370-379` / `:430-438`. |
| 3 | One caller, one consumer | ✅ Verified (caller at `:628`, not `:602`). |
| 4 | `merge` is a field-wise sum | ✅ `report.rs:148-149`. |
| 5 | `--passthrough` + `--retain_unpaired` rejected | ✅ **True** — `cli.rs:948-952`, test at `:2021`. The plan's own hedge can be closed. But see **I2**: being unreachable is exactly why Step 5 must not test it. |
| 6 | Changing `gen_multiome_rows` will affect other tests | ❌ **Over-cautious — I closed this.** See below. |

### Assumption 6 is resolvable now, and the answer is "it doesn't spread"

`gen_multiome_rows` has **five** call sites: `:2396` (serial smoke), `:2477` (the parity test),
`:2589` (truncation), `:2633` (id mismatch), `:2680` (plain-output). I read the assertions in all
four non-parity tests:

- `:2394` smoke — asserts `passthrough_records_checked == 50`, `kept + dropped == 50`, `dropped > 0`,
  `kept > 0`, row-by-row lockstep, and the verbatim 16 bp barcode. **None depends on which mate is
  short.**
- `:2587` truncation and `:2633` id mismatch — assert an error fires and the message names the row.
  Fixture content is irrelevant beyond record count.
- `:2678` plain output — asserts `passthrough_records_checked == 300`, `r1_recs.len() ==
  passthrough_records_kept`, and lockstep. Parallel-side only; no serial comparison.

And mechanically: lengthening R2 on a short-R1 row does not change the pair's verdict. With
`length_cutoff = 20`, R1 at 16 bp still fails, so the pair is still `TooShort`, still counted in
`pairs_removed`, still dropped from the passthrough stream, and still written nowhere (no unpaired
writer). **Only `r2_ok` changes.** The 33% drop rate is therefore not disturbed at all.

So §10's *"the real unknown"* and §11's *"Step 4's blast radius is the real unknown; everything else
is mechanical"* can both be struck: the blast radius is nil, and the plan should say so rather than
hedging to "check first, prefer a variant". (The variant is still the right call — but for I1's
reason, that the non-passthrough test is a better host, not because mutating this one is risky.)

The two things the plan's self-review names as its two biggest risks are therefore both already
settled — Assumption 5 is true and Step 4 is inert — while the risk it does not name (C1: the test
cannot cover the R1 site) is the one that matters.

### Implicit assumptions the plan does not surface

- **That `r1_ok` / `r2_ok` are computed identically on both paths.** They come from the same
  `filters::filter_paired_end` call with the same `UnpairedLengths` on both paths
  (`trimmer.rs:540-547`, `parallel.rs:545-552`), so this is true — but it is the second half of the
  equivalence argument and the plan only states the writer half.
- **That no `PairValidationStats` field is derived from the unpaired counters.** Verified: `merge`
  sums each field independently and no consumer computes a percentage against them (unlike
  `pairs_removed`, which `report.rs` divides by `pairs_analyzed`). So zeroing them cannot shift
  another printed figure.
- **That the worker's count survives a later writer-thread failure.** It does, on both paths, and it
  does not matter — a write failure aborts the run before any report is written.

---

## 3. Efficiency

Nothing to flag, and the plan's own assessment is right: two `usize` increments move into guards
that are already evaluated for the writes that sit beside them. No allocation, no new branch, no
signature change, no change to what crosses the channel between worker and writer thread. Per-worker
accumulation plus a field-wise `merge` is already the established shape for every other counter in
this struct.

The one thing worth noting is what the move *buys* on the hot path: in a no-rescue run, the
`TooShort` arm now does strictly less work per short pair (two branch-and-increments fewer), since
`writer_up_r1` is `None` and the `if let` short-circuits before `r1_ok` is even tested. Immeasurable,
but it is a strict improvement, not a cost.

---

## 4. Validation sufficiency

**Where the plan is strong:** Validation 9 (revert both increments → the `--retain_unpaired`
assertion must fail) is exactly the right instinct and is the test most review processes omit.
Validation 5 (SHA-256 of trimmed output unchanged) correctly pins that this is a stats-only change.
Validation 3 cross-checks the reported figure against the records on disk, which is the only way to
confirm the *semantics* rather than just the parity.

**Where it is not sufficient:**

1. **Validation 7 cannot pass** as long as Step 4 produces only R1-short rows — see **C1**. This is
   the gap that matters: it is a validation whose failure mode is "silently vacuous", which is worse
   than a missing validation because it reports success.
2. **Validation 9 is correctly specified but incompletely load-bearing.** Reverting both increments
   leaves the counters at 0 everywhere. That fails only if some assertion demands a **specific
   non-zero** value. Step 7's guard (`sum > 0`) delivers that, and Step 6's absolute anchor delivers
   it better. But Step 6 is phrased as *"assert the absolute expected count"* without saying on
   which test or with what numbers, and Step 7 places its guard *"in the `--retain_unpaired` run"*.
   If an implementer does Step 5 and Step 9 but treats Steps 6–7 as optional polish, Validation 9
   passes vacuously. Steps 6 and 7 should be stated as **prerequisites of Validation 9**, not as
   separate steps.
3. **No automated coverage of the symptom** — see **I4**. Every proposed assertion is on a returned
   struct; the report text and JSON are checked once, by hand.
4. **Validation 10's test count (`846+`)** is a floor, not an assertion. If Step 4 is done as a
   variant generator plus new tests, the number rises; if a test is accidentally made vacuous (C1),
   the number also rises. The count cannot distinguish those.

**Degenerate "fixes" the validations as written would not catch:**

| Degenerate fix | Caught? |
|---|---|
| Delete both increments, add nothing | Only by Step 6/7's absolute anchors — see point 2 |
| Add the increment inside the `Unpaired` arm but **outside** the `if let` (`if r1_ok { … }`) — i.e. relocate the bug | Yes, by the no-rescue assertion in Step 5 — **provided Step 5's no-rescue run is automated**, not just Validation 1's manual run |
| Gate on `writer_up_r1.is_some()` but drop the `r1_ok` test | Over-counts on the rescue run; caught by Step 6's absolute anchor, **not** by `serial == parallel` (serial would differ) — actually caught, since serial is unchanged |
| Fix R2, forget R1 | **Not caught** — C1 |
| Fix both, but transpose the two counters | **Not caught** if both expected values are equal; fixed by C1's "different counts for the two classes" |

**One more validation worth adding:** assert that the records in the unpaired *files* equal the
counters inside a test, not only in Validation 4's manual `gzip -dc | wc -l`. The library tests
already write real files (`FastqWriter::create` / `run_paired_end_parallel` to paths), so
`read_records(&unpaired_r2).len() == stats.2.r2_unpaired` is one line and it pins the semantics
rather than the parity. That is the assertion that would have caught this defect from the start.

---

## 5. Alternatives

**A. Pass `retain_unpaired: bool` into `classify_paired` and keep the increments there.** Correct,
but adds a parameter to a function already carrying `#[allow(clippy::too_many_arguments)]`, and
leaves the count one function away from the write it is supposed to describe. The plan's choice is
better and for the right reason. Worth one line in the plan so nobody re-litigates it.

**B. Return the counts from `process_pairs` / derive them from the sinks.** Counting records written
into `buf_up_r1` would make the invariant structural rather than conventional, but requires either
counting records in a compressed buffer or threading a counter through the writer — strictly worse.
Reject.

**C. Gate the report block on `retain_unpaired` instead of on `> 0`, and leave the counters alone.**
This is the Perl-faithful shape (O3) and would fix the reported symptom — the block would vanish
without `--retain_unpaired` at any core count. But it leaves `PairValidationStats` divergent between
paths, so the JSON still disagrees with itself across `--cores`, and the parity tests still can't be
made honest. Reject as *the* fix; reconsider as a follow-up for the zero-lines-under-rescue case.

**D. Scope the change to also fix the sibling fixtures.** The plan touches one of three shaped
fixtures (`gen_multiome_rows`). The other two are the non-passthrough parity test's dead
`length_cutoff: 0` (I1) and `tests/integration_passthrough.rs:50-54` (I4). Leaving the other two
means the next person investigating "is this covered?" finds two comments that say the divergence is
known and tolerated, describing a defect that no longer exists. Recommend folding both in — the
integration-test comment is mandatory (it will be factually wrong), the parity-test note is cheap.

---

## 6. Under-specified for an implementer

- **Step 4** — which rows get which shape, and the resulting expected counts. Must be pinned (C1).
- **Step 5** — whether the `--retain_unpaired` run has passthrough active. Must be "no" (I2).
- **Step 6** — *which* test, and the hand-derived numbers. "As in #541, assert the absolute expected
  count" names a pattern, not a value. #541's equivalent is `const EXPECTED_R1: usize = 15 * 30 + 15
  * 34;` at `src/parallel.rs:1712-1714`, derived in a comment. This plan should carry the arithmetic
  the same way, for both counters.
- **Step 8** — the CHANGELOG sentence itself is not drafted. `### Unreleased` exists in
  `CHANGELOG.md` with one #541 entry already in it, so the slot is there; the sentence should name
  both surfaces (I3) and `--clumpify` (O4). One sentence, per the project's rule.
- **Step 9** — "confirm a test fails" should name *which* test, so a vacuous pass is visible.

---

## 7. Action items

### Critical

1. **Reshape Step 4 to exercise both increment sites, with unequal expected counts.** As written the
   fixture rescues R2 only, so Validation 7 (revert the R1 increment → a test fails) cannot fail and
   Step 6's R1 anchor is `0`, which a missing increment satisfies. Specify the row classes and the
   two expected values in the plan. *§5 Step 4 vs §9 row 7; `src/parallel.rs:2283`. Confidence:
   high.*

### Important

2. **Correct §1/§2: there are two `PairValidationStats` parity tests, and retarget Steps 4–7 at the
   other one.** `test_parallel_serial_bp_after_trim_parity_paired` (`src/parallel.rs:1622`, assertion
   `:1704-1709`) asserts the same struct and is vacuous for these fields because
   `parity_test_config_post_adapter` inherits `length_cutoff: 0` (`:1455-1463` → `:1424`). It is
   non-passthrough, already hand-derives absolutes, and already carries a "this test cannot fail
   otherwise" guard — the better host. *Confidence: high.*
3. **Make "the `--retain_unpaired` run has no passthrough" a constraint, not implementer's choice.**
   `src/cli.rs:948-952` rejects the combination (Assumption 5 verified true), so Step 5 as phrased
   would pin behaviour the CLI forbids and contradict §3's own edge-case row. *§5 Step 5, §10.
   Confidence: high.*
4. **Fix §2/§7: the JSON consumer is not gated on `> 0`.** `src/report.rs:1060-1080` emits
   `pair_validation` unconditionally in PE mode, so the fix changes `"r2_unpaired": 4` → `0` for
   existing `--cores >= 2` users on the machine-read surface. Split the "removes lines" framing from
   "changes a number", and write the CHANGELOG sentence to cover both. Also: Validation 6 should name
   the **R2** trimming report, which is where the paired block lands. *Confidence: high.*
5. **Add one automated test of the user-visible symptom, and update the third shaped fixture.** All
   proposed tests assert on returned structs; Validations 1–4 are manual. `tests/integration_passthrough.rs:50-54`
   contains a third copy of the shaping with the same admission in its comment — that comment becomes
   a stale reference to a fixed defect and must be rewritten regardless (project comment rule:
   current state only). That file runs only at `--cores 1` (`:106`); a `--cores 4` paired no-rescue
   case asserting the report has no `Unpaired read` line is the missing guard. *Confidence: high on
   the stale comment; medium on the integration test's cost/benefit.*
6. **State Steps 6 and 7 as prerequisites of Validation 9**, not as independent steps. Without a
   specific non-zero anchor, "revert both increments" passes vacuously — the exact cover-up
   Validation 9 exists to prevent. *§5 Steps 6–7, §9 row 9. Confidence: high.*

### Optional

7. **Promote Assumption 1 to a cited fact and close Assumption 5.** Perl 0.6.11 increments
   `$read_1_printed` after the `print`, inside `if ($retain)` (tag `0.6.11`, `trim_galore:2322-2342`),
   and reports it as "printed" (`:2396`). The fix restores Perl behaviour — a better CHANGELOG
   justification than "intended". Assumption 5 is verified true (`src/cli.rs:948-952`, test `:2021`);
   strike the hedge from §8 and §11.
8. **Strike Assumption 6 and the §10/§11 "blast radius" risk.** All five `gen_multiome_rows` callers
   checked; none depends on the both-short shape, and lengthening R2 on a short-R1 row changes only
   `r2_ok` — the pair verdict, drop rate, and written records are untouched. The variant generator is
   still right, for I1's reason rather than for risk.
9. **Fix the two wrong citations** (`:2360` → `:2474`; `:2338-2342` → `:2531-2538`; the latter is
   before the claimed fn line and so self-inconsistent).
10. **Drop the "keeps the count false if the write fails" rationale from Step 2.** Unreachable on
    this path — the sinks are `Vec<u8>` / `GzEncoder<&mut Vec<u8>>`. Keep the ordering for symmetry
    with `trimmer.rs:584-595`.
11. **Keep the passthrough-rejection clause in Step 3's rewritten comment.** It is still load-bearing
    for the `passthrough_records_dropped += 1` below it; the serial path keeps its copy at
    `trimmer.rs:579-583`.
12. **Name `--clumpify` in the CHANGELOG.** `src/main.rs:1736` routes every `--clumpify` run to the
    parallel path, so users who never typed `--cores` are affected — as #541's entry noted for the
    same reason.
13. **File a follow-up for the `> 0` vs `$retain` gating divergence** (O3): after this fix, a
    `--retain_unpaired` run that rescues nothing prints nothing where Perl prints two zero lines.
    Out of scope; worth not over-claiming Perl parity in the meantime.
14. **Add `read_records(&unpaired_r2).len() == stats.2.r2_unpaired` inside the rescue test.** One
    line, pins the semantics rather than the parity, and is the assertion that would have caught this
    defect originally.

---

## 8. Summary

The fix itself is right, and the plan's central argument — that `Some(writer) && r*_ok` in
`process_pairs` is identical to the serial predicate, so the two increments can simply move — I
verified end to end in source. The measurements reproduce. Both of the plan's self-declared open
risks turn out to be non-risks: Assumption 5 is true, and Step 4's blast radius is nil.

What the plan gets wrong is the part it is explicitly there to get right. Its stated reason for
existing is that the one test covering these fields was shaped to hide the divergence — but there
are two such tests and three such fixtures, and the test plan as specified leaves the R1 increment
site uncovered, which makes its own Validation 7 unsatisfiable. Address C1 and items 2–6 and this is
a clean, small, well-targeted change.

---

*Reviewer A. Independent review; no coordination with Reviewer B. Source read at `86aa301`; no files
in the repository modified other than this report. Measurements re-derived with
`./target/release/trim_galore` and `git show 0.6.11:trim_galore` — both read-only.*
