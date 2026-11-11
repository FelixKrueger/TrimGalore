# PLAN_REVIEW_B — #421 `--fastqc` silent no-op (Reviewer B)

**Reviewed:** `plans/08122026_fastqc-silent-noop/PLAN.md` v1
**Tree:** `dev` @ `3ffcef4`, read-only. No build, no test run, no edits outside this file.
**Verdict:** the refusal design is sound and arm 5's condition is **exactly right** — I re-derived it
from `main.rs` dispatch without using the plan's prose and got the same four terms. But §2.3 / A4's
central claim ("the gate at all 12 call sites is `fastqc || fastqc_args.is_some()`") is **false for 5
of the 12 sites**, and that falsity hides a sixth arm with the same defect the issue is about. One
Critical item, four Important, three Optional.

---

## 1. Logic review

### 1.1 The highest-risk item: arm 5's condition — CONFIRMED CORRECT

Derived from the dispatch, not the plan. `main()` early-returns in this order:

| Order | Site | Condition | FastQC? |
|---|---|---|---|
| 1 | `main.rs:464` | `hardtrim5.is_some()` → `specialty::hardtrim5{,_to_bam}` | **no** |
| 2 | `main.rs:495` | `hardtrim3.is_some()` | **no** |
| 3 | `main.rs:526` | `clock` → `run_specialty_paired` | **no** |
| 4 | `main.rs:548` | `implicon.is_some()` | **no** |
| 5 | `main.rs:580` | `clump_only` (4 drivers) | yes — but see §1.3 |
| 6 | `main.rs:845` | `output_format == UBam` → `run_ubam_output` (3 sub-arms) | yes (`:2242`, `:2399`, `:2522`) |
| 7 | `main.rs:853` | `paired && input.len() == 1` → `run_paired_ubam_single_file` | **no** |
| 8 | `main.rs:867` | `paired` → `run_paired` | yes (`:1778/:1784/:1796`) |
| 9 | SE loop | → `run_single_file` | yes (`:1459`) |

Arm 7 is reachable only when 1–6 have all fallen through, so its CLI-visible predicate is
`paired && input.len() == 1 && output_format != UBam && !clump_only && !hardtrim5 && !hardtrim3 &&
!clock && !implicon`. The plan's four terms are the four that are *not* already refused by guards
1–4, so the condition is exact — neither over- nor under-inclusive. `grep -n "fastqc"
src/specialty.rs` returns nothing; `run_paired_ubam_single_file` (`main.rs:1810-2003`) contains
nothing. Confirmed.

**A3 confirmed, and it is guarded better than the plan thinks.** `main.rs:845` sits *above* `:853`,
so `--paired --output-format ubam <one.bam>` reaches `run_ubam_output` → `run_ubam_output_paired_
single_file` (`main.rs:2006-2017`) → `fastqc::run` at `:2522`. A three-term condition would refuse
it. But that is not only a CI matter: **`tests/integration_ubam_out.rs:572`
(`ubam_out_pe_fastqc_produces_exactly_one_report`) runs exactly `--paired --output-format ubam
--fastqc test_files/ubam_paired_test.bam` through the binary and asserts one `_fastqc.zip`.** So a
dropped `output_format` term goes red under plain `cargo test`, before CI. Same for
`tests/integration_ubam_out.rs:546` (SE uBAM-out) and
`tests/integration_clump_only_ubam.rs:701` (`--clump_only --output-format ubam --fastqc`).
See §5 — this corrects the plan's Remaining-risk #1.

### 1.2 A1 — the five-arm set is complete **for `--fastqc`**

Traced all nine arms above; also checked the shapes the brief named:

- **Multi-pair specialty.** `--hardtrim5/3` loop over `cli.input` inside one arm; `--clock`/
  `--implicon` go through `run_specialty_paired`. All multi-pair shapes live *inside* arms 1–4, so
  `hardtrim5.is_some()` etc. cover them. No new arms.
- **`--demux`.** Called at `main.rs:1468`, *after* the FastQC block at `:1458`, inside
  `run_single_file`. Capable. `--demux --paired` is refused at `cli.rs:1010`, so it cannot reach
  arm 7.
- **`--passthrough`.** Carrier FastQC at `main.rs:1796` inside `run_paired`. Capable. Refused with
  uBAM input (`main.rs:318`) and with uBAM output (`cli.rs:605`), so it cannot reach arm 7 either.
- **`--basename`.** Naming only; no dispatch effect.
- **`--clump_only` shapes.** All four drivers call `fastqc::run` (`clump_only.rs:374`, `:554`+`:555`,
  `:847`, `:1095`). Capable for `--fastqc`. **Not capable for `--fastqc_args` alone** — §1.3.
- **`--retain_unpaired`.** FastQC runs on the `_val_` files only, never on `_unpaired_`. Not the
  #421 defect (a report *is* produced), but it makes `outputs.md`'s "once per output file" wrong in
  a second way — see §1.6.

So A1 holds as stated. It does **not** hold for the disjunction the plan itself insists on.

### 1.3 CRITICAL — §2.3 and A4 are false for 5 of the 12 call sites, and there is a sixth arm

A4: *"The refusal must key on `fastqc || fastqc_args.is_some()` — **Verified** — that is the gate at
all 12 call sites."* It is not. The four `--clump_only` drivers take `fastqc: bool` as a parameter
and gate on it bare:

```rust
// src/clump_only.rs:373, :553, :844, :1094 — all four identical
if fastqc {
    fastqc::run(&output_path, fastqc_args, output_dir, cores.max(1))?;
}
```

and every caller passes `cli.fastqc`, never the disjunction (`main.rs:636`, `:671`, `:715`, `:778`,
`:819`). I checked for a normalisation step that would rescue this — `grep -rn "\.fastqc\s*=" src/`
returns **nothing**, and `cli.rs` has no `requires`/`default_value_if` on `fastqc_args` (`cli.rs:337`
is a bare `#[clap(long = "fastqc_args", allow_hyphen_values = true)]`). So:

```
trim_galore --clump_only --fastqc_args "--quiet" a.fastq   # exit 0, no report, no message
```

is **the #421 defect, on a sixth arm, reachable today** — and it is precisely the failure mode §2.3
warns about ("A refusal keyed on `cli.fastqc` alone would leave `--fastqc_args` silently ignored —
the same defect, one flag over"). The plan diagnoses the shape correctly and then applies the
diagnosis to only one side of the code.

Consequences if this ships as drafted:

1. **The issue is not closed.** `--clump_only --fastqc_args` remains a silent no-op, so the stated
   goal ("the flag either works or says why not") is unmet on one arm.
2. **Two docs lines stay false.** `docs/src/content/docs/modes/clump-only.md:77` — "`--fastqc` —
   runs on the reordered output" — and the `--fastqc_args` help text "Implies --fastqc"
   (`cli.rs:333`) are both untrue on this path.
3. **The plan's own edge-case table is incomplete.** It lists `--clump_only --fastqc a.fq` as
   "accepted — capable" (true) but has no row for the args-only twin, which is accepted *and*
   ignored.
4. **The tripwire deletion leaves a mis-signposted landmine.** After A6 collapses the match to
   `active && allowed.is_empty() → Err`, a future contributor who adds a case
   `--clump_only --fastqc_args "--quiet" a.fastq` gets a red test reading *"a FastQC flag is set but
   no FastQC artefact was written; the allowance is masking a broken case"* — correct diagnosis,
   but the module doc will by then say the incapable arms are all refused at the CLI, so the obvious
   reading is "the tripwire is wrong", not "clump_only's gate is wrong".

**Fix (recommended, small):** add the predicate to `Cli` once and use it everywhere.

```rust
// src/cli.rs — put it next to validate()
/// `--fastqc_args`' "Implies --fastqc" is doc text, not a clap `requires`.
pub fn fastqc_requested(&self) -> bool {
    self.fastqc || self.fastqc_args.is_some()
}
```

Then pass `cli.fastqc_requested()` at `main.rs:636`, `:671`, `:715`, `:778`, `:819` in place of
`cli.fastqc`, and use the same method for the five new guards. Five one-word edits, no signature
change in `clump_only.rs`. This makes A4 true as written, retires the sixth arm, and gives the
refusal and the gates one source of truth — which is the argument the tripwire's own
`fastqc_active` doc comment already makes (`tests/integration_preflight_tripwire.rs:218-222`:
*"Deriving the predicate here, once, from the args keeps one source of truth for it"*).

If Felix would rather keep this PR to the refusal alone, that is a defensible call — but it needs to
be a **stated deferral with a filed issue**, not an unnoticed gap behind an A-row marked "Verified".
Note §7's "Not touched: … every `fastqc::run` call site" would have to change either way.

Minor, same paragraph: §2.1's table claims 12 call sites but lists **11 rows** — `clump_only.rs:555`
(the R2 FastQC in `clump_only_paired`) is missing. The total of 12 is right.

### 1.4 IMPORTANT — where the block is sited changes existing error messages

The plan sites the guards "after the §3.4a `output_format` block", i.e. after `cli.rs:634`. Several
per-flag structural and range checks live *below* that point, and the new guards would pre-empt all
of them:

| Invocation | Message today | Message after (at line 634) |
|---|---|---|
| `--hardtrim5 0 --fastqc a.fq` | `--hardtrim5 must be between 1 and 999, got 0` (`cli.rs:1002`) | FastQC refusal |
| `--hardtrim5 20 --hardtrim3 20 --fastqc a.fq` | cannot combine the two (`cli.rs:~985`) | FastQC refusal |
| `--clock --paired --fastqc r1.fq r2.fq r3.fq` | `--clock mode requires an even number of input files, got 3` (`cli.rs:546` via `:1003`) | FastQC refusal |
| `--implicon --paired --fastqc r1.fq r1.fq` | `Read 1 and Read 2 appear to be the same file` (`cli.rs:556`) | FastQC refusal |

Reporting an unsupported *combination* ahead of a malformed *value* or a structurally invalid input
list inverts the useful order, and the codebase argues the other way in exactly these words at
`main.rs:344-351` ("Sited here … so those report their own defect first"). None of the existing
tests break (`test_hardtrim5_and_hardtrim3_together_rejected`, `test_validate_clock_odd_count_
rejected`, `test_validate_implicon_odd_count_rejected` pass no FastQC flag), so this regresses
silently.

**Fix:** site the block after the `check_restartable_input` loop (`cli.rs:~1017-1020`) — i.e. after
hardtrim range checks, `validate_paired_input` for all three modes, and the `--demux` single-end
check. Nothing in the FastQC guards depends on running early. Cost: zero. Add one V-row asserting
`--hardtrim5 0 --fastqc a.fq` still reports the range error.

### 1.5 IMPORTANT — arm 5's message asserts something that may be false (O1)

O1 says both statements are true for `--paired --fastqc single.fastq`. Only one is. `validate` runs
at `main.rs:231`, format detection at `:249`, and the "only legal if that file is a uBAM" guard at
`:334` — so the ordering claim is right, but the proposed wording

> `--fastqc is not supported when paired FASTQ output is produced from a single interleaved uBAM…`

tells a user who passed a plain FASTQ that their input is an interleaved uBAM. That is a wrong
statement, not a mis-ordered one, and the suggested remedy (`--output-format ubam`) does not work
for them either — `--paired --output-format ubam single.fastq` still fails at `main.rs:334`.

I do **not** think this justifies O1's alternative of splitting the check into `main.rs` per §3.4b.
Reword instead, so the message is true under both readings:

> `--fastqc is not supported with --paired and a single input file: that shape is the interleaved-uBAM
> paired path, which produces FASTQ output with no QC report. Use --output-format ubam (which does
> produce one), run FastQC separately, or pass R1 and R2 as two files.`

Naming the *shape* rather than asserting the *format* costs one clause and removes the false claim.

### 1.6 A6 — the tripwire deletion is safe; A6's justification is understated, its scope is incomplete

Both halves of A6 verify, and one verifies more strongly than claimed:

- **Eight markers, all in scope.** `.fastqc_incapable()` at lines 813, 824, 835, 846, 866, 882,
  1028, 1043 — hardtrim5 ×2, hardtrim3 ×2, clock, implicon, `--paired ip.bam` ×2. Exactly the five
  arms, nothing else. ✓
- **The `(true, false)` branch has never fired — and cannot be made to fire even by the meta-tests.**
  All nine direct `compare(...)` callers (lines 1076, 1097, 1104, 1120, 1182, 1197, 1210, 1255)
  pass `fastqc_capable: true`. So the branch is unexercised *and* unexercisable today, which is a
  stronger argument for deletion than "never fires".
- **The surviving rule is not weakened — it is widened.** Collapsing to
  `if active && allowed.is_empty()` makes the "capable arm with a FastQC flag must produce
  artefacts" rule unconditional across all 53 cases. After the refusal lands, any case with a
  FastQC flag must be on a capable arm by construction, so the widening is exactly right. ✓
  (Caveat: sound only once §1.3 is resolved — see landmine point 4 above.)

**Incomplete in step 3:** the plan lists "the field, the builder, its eight call sites, and the
`(true, false)` match arm". It omits (a) the `fastqc_capable: bool` **parameter of `compare()`**
(line 348), (b) the argument at `Tripwire::run` line 311, and (c) the nine literal `true` arguments
in the `comparison_*` / helper callers. The compiler catches an arity mismatch, so this is a
completeness note rather than a risk — but under `-D warnings` a half-done deletion (field removed,
parameter kept) fails on an unused binding, so spell it out.

V8's "53 cases pass" checks out — `grep -c "#\[test\]"` = 53, and the deletion changes no count.
`the_site_inventory_is_complete` (line 1159) counts `preflight_output_collisions(` occurrences, not
FastQC sites, so it is untouched by this change. ✓

### 1.7 Over-rejection hunt — I could not break a working invocation

I looked specifically for an accepted-today invocation the five guards would refuse:

- Every capable shape (`--fastqc` on SE trim, PE trim, `--passthrough` carrier, `--demux`,
  `--clump_only` ×4, uBAM-out ×3) fails all five conditions. ✓
- `--clump_only --paired --fastqc one.bam`: `!clump_only` excludes it, so `main.rs:597`'s
  better message survives. Verified that `main.rs:597-602` is FASTQ-output-only and that
  `--clump_only --paired --output-format ubam one.bam` goes to
  `clump_only_paired_to_bam_one_pair` (`main.rs:707`), which is capable. ✓
- Deprecated/hidden flags (`--gzip` at `cli.rs:419` etc.) do not affect dispatch. ✓

The only behaviour change beyond the intended refusals is the message-precedence shift in §1.4.

---

## 2. Assumptions

| | Plan's status | My finding |
|---|---|---|
| A1 | Verified | **Holds for `--fastqc`.** Fails for `--fastqc_args` alone: the four `--clump_only` drivers are a sixth arm (§1.3). The grep methodology is the cause — absence of `fastqc` finds missing calls, not weaker gates. |
| A2 | Verified | **Confirmed** independently from dispatch order (`main.rs:845` above `:853`), not from `validate_paired_input`'s comment. |
| A3 | Verified | **Confirmed, and stronger** — also guarded by `tests/integration_ubam_out.rs:572` under plain `cargo test`, not only `ci.yml:1173`. |
| A4 | Verified | **FALSE as stated.** 5 of 12 sites (`clump_only.rs:374/554/555/847/1095`) gate on bare `cli.fastqc`. The clap half is right: no `requires` on `fastqc_args` (`cli.rs:337`). |
| A5 | Verified | **Half true.** All four drivers reach `fastqc::run`, so no refusal is needed — but they are capable only for `--fastqc`, so "needs no refusal" is right while "capable" is not unconditionally. |
| A6 | Verified | **Confirmed, and the branch is unexercisable even by the meta-tests.** Deletion scope is incomplete (§1.6). |
| A7 | Decided | Felix's call; not re-litigated. |

**Implicit assumptions the plan does not surface:**

1. *That the CLI layer can see everything needed.* True for arms 1–4 and, structurally, for arm 5 —
   but arm 5's condition is a proxy for a format fact, and the message must not pretend otherwise
   (§1.5).
2. *That grepping for the callee name enumerates the gates.* This is the assumption that let §1.3
   through. The check that would have caught it is "for each `fastqc::run` call site, what is the
   guard expression?" — 12 sites, five with a different answer.
3. *That `fastqc_args` implying `fastqc` is only a documentation claim about clap.* It is also a
   claim about the code, and the code contradicts it on one path.

---

## 3. Efficiency

Nothing to add — five booleans at argument-validation time, no I/O, messages built only on the error
path. Correct as written. The `fastqc_requested()` method in §1.3 is likewise free (inlined bool
disjunction).

One structural note in the plan's favour: §4's refusal to add a helper is right for the *guards*
(the §3.4a precedent is eight flat `if cond { bail!() }`), but that argument does not extend to the
*predicate*. `fastqc || fastqc_args.is_some()` is currently spelled out at seven sites in `main.rs`
and derived a second time in the tripwire; a method is one source of truth for a fact that has
already been got wrong once.

---

## 4. Alternatives

1. **Fix the `clump_only` gate in this PR (recommended).** Five one-word substitutions plus one
   method. Turns A4 from false to true, closes the sixth arm, and makes `clump-only.md:77` and the
   `--fastqc_args` help text accurate. Cost: `main.rs` enters the diff, contradicting §7's
   "Not touched".
2. **Refuse `--clump_only --fastqc_args` without `--fastqc`.** Consistent with the refusal
   philosophy, but strictly worse: the arm *is* capable, so refusing removes a working feature that
   only needs one flag wired through. Reject.
3. **Defer the sixth arm to a follow-up issue.** Acceptable if the deferral is explicit, the plan's
   A4/A5 rows are corrected, and the tripwire module doc says the clump_only args-only case is a
   known gap (otherwise §1.3's landmine 4 fires). Second choice.
4. **O1's alternative (arm 5 checked in `main.rs` after detection).** Not worth it — §1.5's rewording
   buys the accuracy for one clause instead of a two-file split. Reject, but reject for the wording
   reason, not the ordering one.
5. **Site the guards late in `validate` rather than at line 634.** Strictly better; no downside
   found (§1.4).

---

## 5. Validation sufficiency

**Gaps, highest-risk first.**

- **V-gap 1 (Critical): the disjunction is only tested where it refuses.** V4 tests
  `--fastqc_args` + `--hardtrim5` → refused. Nothing tests `--fastqc_args` on a *capable* arm
  producing a report. That is the axis §1.3's defect lives on. Add:
  - `--clump_only --fastqc_args "--quiet" a.fastq` → an `_fastqc.zip` appears (**red today**);
  - `--paired --output-format ubam --fastqc_args "--quiet" one.bam` → report appears (green today,
    pins `main.rs:2521`).
- **V-gap 2 (Important): no row asserts existing messages keep precedence.** §1.4. Add
  `--hardtrim5 0 --fastqc a.fq` → range error, and `--clock --paired --fastqc r1 r2 r3` →
  odd-count error.
- **V-gap 3 (Important): binary-level refusal coverage collapses from 8 marked shapes to 1.**
  Step 4 proposes one integration case (`--hardtrim5 --fastqc`). Arm 5 is the one whose condition
  is non-trivial, and its *accept* twin already has a binary-level test
  (`tests/integration_ubam_out.rs:572`) while its *refuse* case would have none. Add
  `--paired --fastqc <ip.bam>` → non-zero + message.
- **V-gap 4 (Optional): step 5's `--help` clause has no V-row.** CI greps `--version`, not
  `--help`, so nothing checks it. A `contains` assertion on the long help output is cheap.

**V10 is not runnable as written, and matters less than the plan says.**

- `ci.yml:584` sits inside a step that begins `! command -v fastqc` (line 583). On any machine with
  FastQC on `$PATH` — plausible here — the step fails for reasons unrelated to this change. The
  runnable part is the two lines after it.
- Every step uses `/tmp` and `rm -rf`; in this sandbox `rm -rf` is blocked and `$TMPDIR` is the
  writable root, so they need transcribing rather than copy-pasting.
- More importantly, **three of the four are already covered by `cargo test`**:
  `ci.yml:1146` ≈ `tests/integration_ubam_out.rs:546`; `ci.yml:1173` ≈
  `tests/integration_ubam_out.rs:572`; `ci.yml:584` ≈ tripwire `se_trim_fastqc`
  (line 527) + `se_trim_fastqc_args_outdir` (532) + `se_trim_fastqc_args_only` (549). Only
  `ci.yml:1191` (`--paired --output-format ubam --fastqc <two FASTQ>`) has no local equivalent.

So **Remaining risk #1 is wrong on the facts**: CI is not "the only place three of the four capable
shapes are exercised end-to-end" — V9 (`cargo test`) already catches a mis-specified arm 5. Restate
V10 as "transcribe `ci.yml:1191` only; the rest are covered by V9", and the plan gets a cheaper,
more reliable guard. `cargo test --test integration_ubam_out --test integration_clump_only_ubam`
is the honest local proxy.

**Where wrong results could still pass unnoticed:** nowhere in the trimming path — this change
touches no output bytes, and the Perl validation matrix passes no `--fastqc`, so the byte-identity
invariant is genuinely unaffected (§7 is right). The residual risk is entirely
"a capable arm stops being exercised, or an incapable arm stops being refused", which V-gaps 1 and 3
close.

---

## 6. Docs

`guide/outputs.md:145` is **not** the only line needing the exception. Full sweep of `README.md`,
`docs/src/content/docs/**`, and `--help`:

| Location | Text | Needs |
|---|---|---|
| `guide/outputs.md:143` | "`--fastqc` runs the bundled … library on the trimmed output files … Works on both FASTQ output (default) and uBAM output" | **The broader promise** — a blanket "works on" with no exception. The plan misses it and names only :145. |
| `guide/outputs.md:145` | "For paired-end runs, FastQC is invoked once per output file — so PE FASTQ produces two reports (one per mate)" | As the plan says. Note it is wrong in a **second** way: under `--retain_unpaired` there are four output files and two reports. |
| `modes/clump-only.md:77` | "`--fastqc` — runs on the reordered output" | True for `--fastqc`, false for `--fastqc_args` alone (§1.3). Fix the code or qualify the line. |
| `cli.rs:333` (`--fastqc_args` help) | "Implies --fastqc." | False on the `--clump_only` path (§1.3). |
| `cli.rs:327-331` (`--fastqc` help) | no restriction named | Step 5's clause. Agree with O3 (one clause, name the specialty modes). |

Checked and **clean** — no stale promise: `guide/flags.md` (no FastQC mention at all),
`modes/hardtrim.md:36` "## Compatibility" (covers multi-input, hardtrim5-vs-3, output naming;
says nothing about FastQC), `modes/clock.md`, `modes/implicon.md` (both state the mode exits after
trimming), `modes/passthrough.md:65` (capable arm, accurate), `README.md:113` and `:134` (SE trim
and uBAM output — both capable, both accurate), `guide/overview.md:23` ("Optional demultiplexing and
FastQC" — pipeline-order prose, not a compatibility claim).

`reference/changelog.md:47-52` describes the uBAM-out FastQC fix and mentions "`--fastqc` (or
`--fastqc_args`)". It is historical and the synced changelog page should not be edited — but note
that entry is the origin of A4's belief, and it was accurate for the sites it was about.

---

## 7. Action items

### Critical

1. **Resolve the sixth arm: `--clump_only --fastqc_args` without `--fastqc` is a live silent
   no-op.** The gate at `clump_only.rs:373/553/844/1094` is bare `fastqc: bool`, fed `cli.fastqc`
   at `main.rs:636/671/715/778/819`. Either add `Cli::fastqc_requested()` and pass it at those five
   call sites (recommended — 5 one-word edits), or state the deferral explicitly and file it.
   Correct A4 ("false for 5 of 12 sites") and A5 either way, and drop `main.rs` from §7's
   "Not touched" if fixing.

### Important

2. **Site the guard block late in `validate`** (after the `check_restartable_input` loop,
   `cli.rs:~1020`), not after `cli.rs:634`. At 634 it pre-empts the `--hardtrim5` range check, the
   hardtrim5+hardtrim3 conflict, and `validate_paired_input`'s odd-count / R1==R2 / duplicate-pair
   messages for all three paired modes (§1.4). Zero cost.
3. **Reword arm 5's message so it does not assert the input is a uBAM** (§1.5). Name the shape
   (`--paired` with one input file) and add "or pass R1 and R2 as two files" to the remedy, so the
   FASTQ-single-file user gets advice that works. Keeps O1's decision, removes its false statement.
4. **Add the two missing capable-arm `--fastqc_args` V-rows** (V-gap 1) — the
   `--clump_only --fastqc_args` one is red today and is the test that makes item 1 self-proving.
5. **Add the arm-5 refusal to the integration case** (V-gap 3): `--paired --fastqc <ip.bam>` →
   non-zero. It is the only refusal whose condition is non-trivial, and the only one whose accept
   twin has a binary-level test but whose refuse case would not.

### Optional

6. **Restate V10 and Remaining-risk #1.** `cargo test` already covers `ci.yml:584`, `:1146` and
   `:1173` (tripwire `se_trim_fastqc*`, `integration_ubam_out.rs:546`, `:572`); only `:1191` is
   CI-only. Also flag that `ci.yml:583`'s `! command -v fastqc` makes step `:584` unrunnable on a
   machine with FastQC installed.
7. **Complete step 3's deletion list** (§1.6): `compare()`'s `fastqc_capable` parameter (line 348),
   the argument at line 311, and nine literal `true`s in the `comparison_*` callers.
8. **Two small factual fixes.** §2.1's table lists 11 rows while claiming 12 — `clump_only.rs:555`
   is missing. And add `guide/outputs.md:143` to the docs task (step 6), plus a note that :145 is
   also wrong under `--retain_unpaired`.

---

*Reviewer B, independent. No file under the repo was modified other than this report.*
