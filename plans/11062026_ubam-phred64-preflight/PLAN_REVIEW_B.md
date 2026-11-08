# Plan Review B — `plans/11062026_ubam-phred64-preflight/PLAN.md`

**Reviewer:** B (independent, fresh context)
**Date:** 2026-11-07
**Target:** `/Users/fkrueger/Github/TrimGalore/plans/11062026_ubam-phred64-preflight/PLAN.md`
**Tree:** `dev` @ `49152c1`, release binary from the same commit (`--version` confirms `49152c1 — macos/aarch64`).

*Status: complete.*

---

## 0. Verification log (what I re-derived rather than trusted)

### 0.1 Line citations in §2 — all correct

| Plan claim | Verified |
|---|---|
| `src/main.rs:369` uBAM-input guard | ✅ `if cli.phred64 && any_bam {` is line 369 |
| `src/main.rs:386` the condition being widened | ✅ `if cli.phred64 && cli.fastqc_requested() {` is line 386 |
| `src/main.rs:26` `PHRED64_PREFLIGHT_READS = 10_000` | ✅ |
| `src/main.rs:36` `lowest_quality_char` | ✅ |
| `src/main.rs:553` startup NOTE | ✅ (NOTE body spans 553–564) |
| `src/cli.rs:605` `fastqc_requested()` | ✅ `self.fastqc \|\| self.fastqc_args.is_some()` |
| `src/cli.rs:630` `--passthrough` + uBAM rejection | ✅ exact line |
| `naming::ensure_output_dir` at `src/main.rs:546` | ✅ |
| §3.4a exclusion set | ✅ `--dont_gzip`, `--clumpify`, `--passthrough`, `--clock`, `--implicon`, `--demux`, `--retain_unpaired`; `--hardtrim5/3` and `--clump_only` are indeed **absent** |

### 0.2 Ordering — the siting claim holds

Every dispatch branch that can write output begins **after** line 386 *and* after `ensure_output_dir`
(546): `--hardtrim5` at 597, `--hardtrim3` at 628, `--clock` at 668, `--implicon` at 690,
`--clump_only` at ~761, `run_ubam_output` at 982. So "refused before anything is written" holds for
every mode the widened gate now covers, including the early-returning branches the brief singled out.

### 0.3 The §2 measurement table — reproduced

Own fixture (4 × 40 bp, all `I` = ASCII 73 = Q40 under +33), own runs, `samtools` read-back:

| Invocation | Measured exit | Measured result | Plan says |
|---|---|---|---|
| `--output-format ubam --phred64` (default `-q 20`) | 0 | `.bam` written, **0 records** | matches |
| `--output-format ubam --phred64 -q 5` | 0 | 4 records, QUAL = **9** | matches (Q40→Q9) |
| `--output-format ubam --phred64 -q 0` | 0 | 4 records, QUAL = **9** | matches |
| `--hardtrim5 20 --output-format ubam --phred64` | 0 | `in.20bp_5prime.bam`, 4 recs, QUAL = **9** | matches |
| `--clump_only --output-format ubam --phred64` | 0 | `in_clumped.bam`, 4 recs, QUAL = **9** | matches |
| `--hardtrim5 20 --phred64` (FASTQ out) | 0 | `IIIIIIII…` verbatim | matches |

I additionally checked `--clump_only --phred64` to FASTQ: also verbatim `IIII…`. The plan's §2
conclusion #1 ("the exposed modes are the ones that do no quality arithmetic … inert on FASTQ
output") is sound.

### 0.4 The mixed-quality flattening result — mapping reproduces exactly, one figure does not

Fixture of 48 bases cycling Q0,5,10,15,20,25,30,31,32,35,38,40 (+33), `--output-format ubam
--phred64 -q 0`:

```
in : 0 5 10 15 20 25 30 31 32 35 38 40
out: 0 0  0  0  0  0  0  0  1  4  7  9
```

Identical to the plan. 12 distinct input levels → 5 distinct output levels, i.e. **7 destroyed** —
the plan's figure, exactly.

**But "48% of bases flattened to Q0" is not reproducible from the plan as written.** My fixture
gives **66.7%**. The percentage is a property of the fixture's quality *composition*, which the plan
never states, so the number cannot be checked or re-derived by a reader. The load-bearing claims
(the mapping, the floor-at-0, the 7-of-12) all hold; only this one figure is unanchored. See
Optional O-1.

### 0.5 #362 regression guard — the risk is real and currently passing

Genuine Phred+64 fixture (lowest char `B` = ASCII 66 = Q2) through `--output-format ubam --phred64
-q 0`: exit 0, 3 records, `QUAL` = `2 2 3 3 4 4 5 5 6 6` against expected `ASCII − 64` =
`2 2 3 3 4 4 5 5 6 6`. Correct today, and the widened condition leaves it correct because the lowest
byte is 66 ≥ 64. Validation #4 is guarding something that genuinely works and genuinely could break.

### 0.6 `--passthrough` + uBAM — the plan's assumption checks out

`src/cli.rs:630` rejects `--passthrough` under `--output-format ubam` unconditionally, and
`validate()` runs at `src/main.rs:322`, well before the scan. So in the uBAM case the `.chain()`
arm is genuinely dead, exactly as §8 claims. (It stays load-bearing for the `--fastqc` arm, where
`preflight_scans_the_passthrough_carrier` covers it.)

### 0.7 Prior-art reading — one quote exact, one framing wrong

The #536 commit body says, verbatim: *"because refusing there would reject `--hardtrim5/3` and
`--clump_only` runs that pass the scores through untouched. #532 keeps that wider question."* The
plan quotes this accurately and its reading of #358 → #362 matches the code. But see **I-1** below:
the §2 sentence about which configuration "already had a backstop" is wrong, and wrong in a way
#536 itself already corrected once.

---

## 1. Critical

### C-1 — §2 conclusion #1 is refuted: FASTQ output is *not* inert, and the plan's whole scope rationale rests on it

The plan states:

> **The exposed modes are the ones that do no quality arithmetic.** The plain trim path at default
> `-q 20` largely self-defends — the shift pushes everything under the cutoff and the run empties
> itself visibly. […] Gating on **output format** protects the corrupting invocations and leaves
> the inert FASTQ ones working.

The plan measured `-q 0` and `-q 20` and interpolated between them. There is a band in between
where **FASTQ output silently deletes good bases**. Measured, this tree, this binary, Phred+33
fixture of 60 bp (40 × `I` = Q40, then 20 × `D` = Q35), `--length 0`:

| `--phred64 -q N`, FASTQ output | Exit | Reads kept | First read length |
|---|---|---|---|
| `-q 1` | 0 | 3 | 60 — intact |
| `-q 5` | 0 | 3 | **40 — 20 good bases silently deleted** |
| `-q 9` | 0 | 3 | **40 — 20 good bases silently deleted** |
| `-q 10` | 0 | 3 | 0 — emptied |
| `-q 20` | 0 | 3 | 0 — emptied |

Control run (`-q 5`, no `--phred64`): 60 bp retained. So `--phred64` alone, on FASTQ output, exit 0,
removed a third of every read.

The mechanism: a surviving base needs `ascii − 64 ≥ q`. Real Phred+33 Illumina tops out around
ASCII 70–73 (Q37–Q40), so the shifted ceiling is 6–9. That gives a **silent band of roughly
`2 ≤ -q ≤ (max_ascii − 64)`** — `[2, 9]` for Q40-capped data, `[2, 6]` for NovaSeq-binned Q37 data.
Above the band everything shifts below the cutoff and the run empties (genuinely loud — and this
part of the plan's reasoning *is* sound, and generalises: no real Phred+33 byte reaches ASCII 84,
so `-q 20` always empties). Below the band nothing trims. Inside it, exit 0 and truncated reads.

Why this is Critical rather than Important:

1. **It falsifies the premise the plan uses to choose its scope.** §2 conclusion #1 and the §5.1
   condition are the same argument. The claim "FASTQ output passes scores through verbatim" is true
   of the quality *bytes* and false of the trim *decisions*. On uBAM the false declaration corrupts
   QUAL values; on FASTQ it corrupts read lengths. Both are the same bug — a false declaration
   feeding arithmetic — and the plan counts only one of them as arithmetic.
2. **It makes §10's NOTE recommendation actively harmful.** §10 says "with the dangerous uBAM paths
   refused, most of the NOTE's second half describes an outcome that is now unreachable for
   mislabelled input. Recommend trimming it." That premise is false. The current NOTE at
   `src/main.rs:555-558` is the *only* warning a user in the silent band ever sees, and it is
   already correct about it — *"Wherever quality scores are read (**quality trimming**, and BAM QUAL
   on uBAM output), Phred+33 data taken as Phred+64 loses 31 from every score, floored at 0"*.
   Trimming the NOTE to "the FASTQ-output case plus one clause on the uBAM refusal", as §10 directs
   ("Taken as: trim, minimally"), risks deleting the only live warning for the one path that stays
   silent after this change.
3. **It changes what §10's second open question is.** §10 dismisses FASTQ output as "`-q 0` still
   passes through with scores intact but a wrong `Quality encoding type selected` line … Harmless to
   the data; out of scope here." At `-q 0` that is accurate. The question is not about `-q 0`.

**What I'd ask for.** Not necessarily a wider gate — refusing `--phred64` on all FASTQ output would
reject the `--hardtrim5/3` runs #536 deliberately preserved, and that trade-off is still the user's
call. But the plan must (a) correct §2 conclusion #1 so the record is right, (b) withdraw §10's
"trim the NOTE" direction or re-derive it, keeping the `quality trimming` clause intact, and
(c) state the residual explicitly: *after this change, `--phred64` + FASTQ output + `-q` inside the
band still silently truncates, and nothing catches it.* 

**On whether to widen the gate — I checked the obvious extension and it does not work.** The
tempting condition is `… || cli.quality > 0`. It is wrong: `src/cli.rs:49` declares
`#[clap(short = 'q', long = "quality", default_value = "20")] pub quality: u8`, so `cli.quality > 0`
is true on *every* run that does not explicitly pass `-q 0` — including every `--hardtrim5/3` and
`--clump_only` run, which is exactly what #536 preserved. Nor can the implementation ask whether
`-q` was user-set: `src/cli.rs:782-784` already records that `-q` is *silently accepted* under
`--clump_only` precisely because "non-Option clap defaults … can't distinguish user-set from
default without threading `ArgMatches::value_source`".

The condition that would actually work is mode-based — "quality trimming is live", i.e. not
hardtrim, not `--clump_only`, not `--clock`, not `--implicon`. But that is verbatim the rule the
`src/main.rs:358-363` comment already rejected on durability grounds: *"a mode-dependent rule
('rejected unless you're in hardtrim, or clump-only-to-FASTQ, …') is worse to document, worse to
test, and one refactor away from being wrong."* So I am **not** asking for a wider gate. I am asking
the plan to stop claiming the residual does not exist. Option 1 as scoped is a defensible increment;
its justification needs rewriting, and the NOTE needs to survive intact.


### C-2 — §5.5's helper list omits `assert_preflight_refused`, and test 3 then passes for the wrong reason (demonstrated)

§5.5 says: *"Per §9. Same file, same helpers (`write_phred33_fixture`, `run_raw`, `fixture_copy`)."*
That list omits the one helper that makes these tests mean anything:
`assert_preflight_refused(stderr, out_dir)` at `tests/integration_fastqc_phred64.rs:536`, whose own
doc comment states why it exists — *"FastQC refuses the same inputs later, naming the same files, so
a non-zero exit alone does not distinguish the two."*

For the new tests the hazard is worse than for #536's, because `--clump_only` and `--hardtrim5/3`
carry large rejection matrices of their own. I demonstrated it. The existing pre-flight tests all
pass `-q 0 --length 0`; an implementer following §5.5 ("same helpers") will copy that pattern. On
test 3 it does this:

```
$ trim_galore --clump_only --output-format ubam --phred64 -q 0 --length 0 in.fastq
exit=1
Error: --clump_only does not filter; --length is not compatible
```

Exit 1, no output directory, no BAM — so a test asserting "non-zero exit, nothing written" (which is
literally what §9 row 3 prescribes) **passes green today, before the feature is implemented at all**,
and keeps passing under every mutant in §9's mutation list. `--clump_only` rejects `--length` at
`src/cli.rs:826`.

I checked the other two for the same trap and they are clean: `--hardtrim5 20 … -q 0 --length 0`
exits 0 (both flags inert in that mode), and plain trim `… -q 0 --length 0` exits 0. **Only test 3
is affected**, which is exactly the test covering one of the two headline silent paths.

**Fix:** add `assert_preflight_refused` to §5.5's helper list and require every new refusal test to
call it, and drop `--length 0` from the `--clump_only` case. The plan's own §9 "Mutation check"
preamble invokes the handoff lesson about tests that "passed for the wrong reason"; this is that
lesson recurring inside the plan that cites it.

---

## 2. Important

### I-1 — §2's "the one configuration that already had a backstop" is wrong, and #536 already corrected this exact claim

§2 says the `--fastqc` gate *"landed on the one configuration that already had a backstop
(`src/fastqc.rs:79` refuses any file holding a character below ASCII 64), and left uBAM output —
which has none, and whose wrong answer persists on disk — unchecked."*

The `--fastqc` gate does not pick out one configuration. `fastqc_requested()` is orthogonal to
`output_format`, so the gate already covers `--fastqc` **plus uBAM output** — and on that half there
is no backstop. All three uBAM-output FastQC call sites pass `phred64: false`:

- `src/main.rs:2436` (`run_ubam_output_single`)
- `src/main.rs:2595` (`run_ubam_output_paired_two_files`)
- `src/main.rs:2715` (`run_ubam_output_paired_single_file`)

each as `fastqc::run(&output_path, …, cli.cores, false)`. With `phred64: false` the engine reads the
BAM as Phred+33 and accepts it, so nothing contradicts the declaration.

This is not a new discovery — it is in #536's own commit body: *"'FastQC rejects such a file either
way' was false on every BAM-output path, where the call sites pass `phred64: false` and FastQC
accepts the file. Gone."* The plan reinstates the same mistaken premise in a different sentence.

No functional consequence (the gate fires on those runs regardless of which arm matches), but §2 is
the written record for this decision and the claim is load-bearing in its argument. Rewrite as: the
`--fastqc` gate covers FASTQ output, where FastQC is a late backstop, and uBAM output, where it is
not; what is left unchecked is uBAM output *without* `--fastqc`.

### I-2 — the third prescribed mutation cannot fail the test it names

§9 prescribes: *"flip `lowest < 64` to `lowest <= 64` → test 4 must fail."*

Test 4's fixture is `good_p64.fastq`, which `fixture_copy` sources from
`test_files/phred64_test.fastq`. I checked its contents: exactly two distinct quality characters,
`B` (66) and `h` (104); **lowest = 66**, as the plan says. Under the mutant the predicate becomes
`66 <= 64`, which is false, so the gate still does not bail and **test 4 passes under the mutant**.
The mutation and the assertion point the same way — the precise failure the plan's own mutation-check
preamble warns about two paragraphs earlier.

The mutant *is* caught, by a test that already exists:
`preflight_accepts_ascii_64_which_is_phred64_q0` (`tests/integration_fastqc_phred64.rs:508`), whose
fixture's low character is `'@'` = ASCII 64 exactly. Name that test instead. No new test is needed
for this mutant.

### I-3 — the plan surveys no uBAM test file, and §9 row 4 duplicates four existing guards

§2's "Existing code this must align with" table and §11's "Edge cases checked" never mention
`tests/integration_ubam_out.rs` or `tests/integration_clump_only_ubam.rs` — the two files holding
every existing test of the exact combination being changed. There are **four** genuine-Phred+64 →
uBAM guards already in the tree, all on `test_files/phred64_test.fastq` (lowest 66), all asserting
exit 0 plus raw-Phred `QUAL` values:

| Test | Path | Covers |
|---|---|---|
| `phred64_ubam_out_se_stores_true_phred` | `tests/integration_ubam_out.rs:632` | SE trim → uBAM |
| `phred64_ubam_out_pe_two_file_stores_true_phred` | `tests/integration_ubam_out.rs:799` | PE two-file trim → uBAM |
| `phred64_ubam_out_specialty_stores_true_phred` | `tests/integration_ubam_out.rs:666` | `--hardtrim5 20` → uBAM |
| `phred64_clump_only_ubam_out_stores_true_phred` | `tests/integration_clump_only_ubam.rs:779` | `--clump_only` → uBAM |

I ran all three files at `49152c1`: **75 passed, 0 failed** (44 + 17 + 14). Every one of these stays
green under the proposed change, because 66 >= 64.

Three consequences:

1. **§9 row 4 is redundant.** `phred64_ubam_out_se_stores_true_phred` *is* row 4 — same fixture,
   same flags (`--phred64 --output-format ubam -q 0 --length 0`), and it already asserts the exact
   raw-Phred vector (24 x Q40 then 10 x Q2) rather than just "exit 0". Row 4 should become "confirm
   these four stay green", not "add a test".
2. **§5.5 points the new tests at the wrong file for row 4 as written.**
   `tests/integration_fastqc_phred64.rs` has no BAM reader at all — its helpers (`basic_stat`,
   `per_base_mean`) read FastQC zips. Row 4's "`samtools view` shows true Phred" is not
   implementable there without importing `noodles`. `tests/integration_ubam_out.rs:616` already has
   `fn quals(path: &Path) -> Vec<Vec<u8>>`, and `:601` already has an `issue #358: --phred64
   quality-encoding guards` section that is the natural home.
3. **It strengthens the plan's §8 assumption.** §8 argues a blanket `--phred64` + uBAM refusal "would
   take back #362". That is not just an argument — three of the four tests above would fail outright
   under a blanket refusal. Worth saying, because it is the plan's best defence of its own scope.

Recommend §9 gain a row: "the four existing genuine-Phred+64 → uBAM guards stay green", naming them,
and that the new refusal tests for the hardtrim and clump_only arms sit next to their positive
counterparts rather than in the FastQC file.

### I-4 — §11's reading of the BAM skip is backwards; it is load-bearing exactly where the plan calls it decorative

§11 says: *"The guard at 369 removes uBAM inputs before the scan, so the BAM skip inside the loop is
belt-and-braces for the `--passthrough` arm rather than load-bearing."*

The opposite holds. For `cli.input` entries the skip is indeed redundant (369 handles them). For the
`--passthrough` carrier it is the only thing standing between the scan and a misattributed error: the
carrier's own uBAM refusal lives at **`src/main.rs:456-462`**, i.e. *seventy lines after* the scan at
386. Without the skip, a uBAM carrier under `--fastqc` with FASTQ output gets blamed by the scan for
a quality byte `BamReader` synthesised rather than one the file holds, and never reaches its real
refusal.

Again this is #536's own finding, verbatim from its commit body: *"A uBAM `--passthrough` carrier was
blamed for a quality character the reader synthesises … because the scan ran above the carrier's own
format check. uBAM inputs are now skipped, so the carrier reaches its real refusal."*

Behaviourally harmless here — the plan retains the skip unchanged — but the recorded rationale is
what a future refactor will read before deleting it. Note also that no test pins this:
`preflight_scans_the_passthrough_carrier` uses a *FASTQ* carrier. (Pre-existing gap, not introduced
by this plan; worth one line in §9's "Not covered".)

### I-5 — §5.3 and §11 overstate the NOTE defect: it is ambiguous, not false

§5.3 asserts the NOTE *"qualifies `--clump_only` with 'writing FASTQ' and leaves `--hardtrim5/3`
**unqualified**"*, concludes *"on the exact invocation that corrupted silently, the program printed a
reassurance that was false"*, and §11 records this as *"It is **wrong** for `--hardtrim5/3` + uBAM,
so §5.3 is now a required correctness fix in its own right."*

The actual text (`src/main.rs:559-561`) is:

> Modes that do no quality arithmetic, such as --hardtrim5/3 and --clump_only writing FASTQ, pass
> the scores through untouched.

Two readings are grammatically available, and the plan commits to the weaker one. On the natural
reading, "writing FASTQ" is a reduced relative on the whole coordinated phrase ("--hardtrim5/3 and
--clump_only, writing FASTQ"), which is *correct* — I verified both: `--hardtrim5 20 --phred64` to
FASTQ emits `IIII…` verbatim, and `--clump_only --phred64` to FASTQ likewise. The sentence is also
self-limiting: its head is "Modes that do no quality arithmetic", and hardtrim-to-uBAM *does* do
quality arithmetic (`BamWriter` subtracts the offset), so it is not in the described set.

Clarifying an ambiguous safety warning is still worth doing. But ranking it a standalone correctness
defect "required regardless of the rest of this plan" is not supported, and §11 lists it as a finding
of the plan's own self-review, which gives it more weight than it earned. Downgrade the framing to
"ambiguous; make the qualifier distributive".

This matters more than a wording quibble because of **C-1**: §5.3 and §10 together move toward
shortening the NOTE, and the clause C-1 shows is load-bearing (`quality trimming`) sits in the half
§10 proposes to cut.

---

## 3. Optional

### O-1 — "48% of bases flattened to Q0" is not reproducible from the plan

§2's mixed-quality block gives three figures. Two reproduce exactly on my own fixture: the mapping
`0 0 0 0 0 0 0 0 1 4 7 9` is byte-identical to the plan's, and 12 distinct input levels collapse to
5, i.e. **7 destroyed** as stated. The third — "48% of bases flattened to Q0" — gave **66.7%** on
mine. That is not a contradiction: the percentage is a property of the fixture's quality
*composition*, which the plan never states. Either name the fixture's composition or drop the
percentage; as written it is the one number in §2 a reader cannot check, and it sits next to two that
are exactly right.

### O-2 — `src/fastqc.rs:79` is cited as the check; it is the error-message mapper

§2 cites `src/fastqc.rs:79` as the site that "refuses any file holding a character below ASCII 64".
Line 79 is `runner::run(&config, …).map_err(|code| …)` — it rewrites a non-zero exit into a
`--phred64`-flavoured message. The refusal itself happens inside `fastqc-rust`'s own `detect()`,
which trim_galore does not own. Minor, but it matters for I-1: nothing in `src/` enforces this, so
the "backstop" is a property of the pinned dependency and of `phred64: true` being passed.

### O-3 — genuine Solexa input + uBAM output becomes a dead end with a misdirected remedy

The message's closing clause is unconditional: `"Drop --phred64 if the data is Phred+33."` For a real
Solexa / Illumina 1.0 file (ASCII 59-63), dropping `--phred64` is *wrong* — it would read the data as
Phred+33 and be off by 31 the other way. Measured on a `;` (ASCII 59) fixture, the full message is:

> The input is therefore Phred+33 or Solexa / Illumina 1.0, which allows ASCII 59-63 for negative
> scores and which Trim Galore does not support, flooring them to Q0 — read as Phred+64 every score
> would drop by 31. Drop --phred64 if the data is Phred+33.

Today such a user can get uBAM output (negatives floored at 0, which is what the message itself says
Trim Galore does). After this change they are hard-blocked, with no stated next step. The real
workaround — FASTQ output, then `samtools import` — exists but is not mentioned. Affected population
is tiny (pre-2011 data) and refusing is defensible; worth one clause in the message or in §3's edge-case
table so the decision is deliberate rather than incidental.

### O-4 — §9's second mutation adds nothing over the first for tests 1-3

`||` → `&&` makes the condition `phred64 && (fastqc_requested() && ubam)`. Tests 1-3 pass
`--output-format ubam --phred64` with no `--fastqc`, so the gate is off — the same outcome as
mutation 1 (dropping the `UBam` arm). The two mutants are indistinguishable from tests 1-3's point of
view. The `&&` mutant is worth keeping, but its value is that it fails the *existing* `--fastqc`
tests (`preflight_refuses_mislabelled_input_before_writing_anything` and friends, which use FASTQ
output, so `ubam` is false). Say that, so the loop's output is interpreted correctly.

### O-5 — CHANGELOG section is ambiguous; `### Unreleased` sits empty above `### Version 2.4.0`

§5.6 says "Add one entry under `### Version 2.4.0`". `CHANGELOG.md:4-6` has an **empty**
`### Unreleased` heading immediately above `### Version 2.4.0 (Release on 3 November 2026)`. Given
§10's own observation that nothing shipped on the 3rd and master is behind dev, placing a new entry
under a dated-but-unreleased heading while an empty `Unreleased` heading sits above it will read
oddly to the next person. State the choice and why in one clause, or the implementer will guess.

### O-6 — a second CHANGELOG sentence is also affected, not just the #532 one

§5.6 correctly identifies that the #532 entry's closing sentence ("`--phred64` without `--fastqc` is
unaffected") becomes false. The **#500** entry two paragraphs below also claims: *"A Phred+33 file
mistakenly run with `--phred64` and `--fastqc` now fails at the FastQC step rather than producing a
report."* Per I-1 that is already false on every uBAM-output path (`phred64: false` there), and after
this change the failure comes from the pre-flight, not the FastQC step, on both. Worth one word
("fails up front") while the file is open.

### O-7 — §6 reasons about a cost the implementation need not pay

§6 accounts for "one extra `detect_input_format` + bounded scan per input". The extra detection is
avoidable for the `cli.input` arm: `input_formats` is already computed at `src/main.rs:343` over
exactly `cli.input`, in order, and the loop re-detects at `src/main.rs:392`. Only the `--passthrough`
carrier genuinely needs a fresh call. Pre-existing and genuinely negligible — raised only because §6
chose to reason about it, and because zipping against `input_formats` would also make the
"already detected" relationship visible to the next reader.

### O-8 — "option 1 of 4" names no other option

The **Decision** line says "option 1 of 4". The plan never states what options 2-4 were, so the
decision cannot be audited from the plan, which is the artifact meant to carry it. One line each
would fix it. (My C-1 note on why `cli.quality > 0` fails is offered in that spirit — if a "verify
whenever quality trimming is live" option was among the four, the plan should record that it was
considered and why it was rejected.)

---

## 4. Review-area summary

**Logic.** The core claim — site the widened condition at `src/main.rs:386`, refuse before anything is
written — holds, and I verified it rather than taking it: every writing dispatch branch
(`--hardtrim5` 597, `--hardtrim3` 628, `--clock` 668, `--implicon` 690, `--clump_only` ~761,
`run_ubam_output` 982) begins after both 386 and `ensure_output_dir` at 546, including the
early-returning ones the brief flagged. The set of affected modes is complete, and the
`--passthrough` exclusion (§8) is exactly right. What does not hold is the *rationale*: §2's
conclusion #1 is refuted by measurement (C-1), §2's backstop sentence is wrong (I-1), and §11's
reading of the BAM skip is inverted (I-4).

**Assumptions.** The two "fixed rules" in §8 are sound — Phred+64 genuinely cannot hold a byte below
64, and BAM genuinely stores raw Phred. The `--passthrough` assumption is verified and correct. The
unstated assumption that does the damage is "FASTQ output performs no quality arithmetic", which is
true of quality *bytes* and false of trim *decisions*.

**Efficiency.** Fine and correctly reasoned; the only note is O-7. A refused run does strictly less
I/O than today, as claimed.

**Validation sufficiency.** The weakest part of the plan. One prescribed mutation cannot fail the
test it names (I-2); one prescribed test passes green today, before implementation, because of a flag
incompatibility (C-2); one prescribed test duplicates four existing guards the plan never surveys
(I-3); and nothing covers the paired uBAM-output refusal path, although paired uBAM output has its
own two drivers and its own collision pre-flight. Baseline measured green: 75 tests across the three
relevant integration files at `49152c1`.

**Alternatives.** I looked for a wider gate and could not find a clean one — see the note under C-1.
`cli.quality > 0` is wrong (default 20, and silently accepted under `--clump_only`); a mode-based
"quality trimming is live" predicate works but is the exact rule `src/main.rs:358-363` already
rejected on durability grounds. So **option 1 as scoped is the right increment**; my objection is to
the justification and the validation, not the change.

---

## 5. Action items

**Critical**
1. **C-1** — correct §2 conclusion #1: FASTQ output silently truncates reads when `-q` falls in
   roughly `[2, max_ascii - 64]` (measured: `-q 5` and `-q 9` delete 20 of 60 bases, exit 0).
   Withdraw or re-derive §10's "trim the NOTE" direction and keep the `quality trimming` clause of
   `src/main.rs:555-558` intact. State the residual explicitly in §8.
2. **C-2** — add `assert_preflight_refused` to §5.5's helper list, require it in every new refusal
   test, and drop `--length 0` from the `--clump_only` case (it is rejected at `src/cli.rs:826`, so
   §9 row 3 as prescribed passes green before the feature exists).

**Important**
3. **I-1** — rewrite §2's backstop sentence; all three uBAM-output FastQC call sites pass
   `phred64: false`, so uBAM output has no backstop whether or not `--fastqc` is given.
4. **I-2** — §9's `<` -> `<=` mutation cannot fail test 4 (lowest byte is 66). Point it at the
   existing `preflight_accepts_ascii_64_which_is_phred64_q0`.
5. **I-3** — survey `tests/integration_ubam_out.rs` and `tests/integration_clump_only_ubam.rs` in §2;
   turn §9 row 4 into "these four stay green" and site the new hardtrim/clump_only refusal tests
   beside their positive counterparts, where `quals()` already exists.
6. **I-4** — fix §11: the BAM skip is load-bearing for the `--passthrough` arm (carrier refusal is at
   `src/main.rs:456`, seventy lines downstream), not belt-and-braces.
7. **I-5** — downgrade §5.3/§11's "the NOTE is wrong" to "the NOTE is ambiguous"; make the qualifier
   distributive rather than recording a correctness defect that the text does not clearly contain.

**Optional**
8. O-1 anchor or drop the "48%" figure · O-2 re-cite `src/fastqc.rs:79` · O-3 decide the Solexa +
   uBAM dead end deliberately · O-4 state what the `&&` mutant actually proves · O-5 pick a CHANGELOG
   section explicitly · O-6 the #500 entry needs a word too · O-7 reuse `input_formats` ·
   O-8 name options 2-4.

---

## 6. Scope of my verification

Everything above was run against `dev` @ `49152c1` with the pre-built release binary from that
commit (`--version` prints `49152c1 — macos/aarch64`). I made no edits under `src/`, `tests/`,
`docs/` or `CHANGELOG.md`; the only files I wrote are this report and fixtures under
`/tmp/claude-501/revB_phred64_20261107/`. I did **not** apply any source mutation — I-2 is derived
arithmetically (`66 <= 64` is false) rather than measured, which for that claim is stronger than a
mutation run. I did not run the full `cargo test` suite, `fmt`, or a cold clippy; I ran the three
integration files that touch `--phred64` (75 tests, all green). CI's rustc 1.99 lints remain
unverifiable locally, as §9 notes.
