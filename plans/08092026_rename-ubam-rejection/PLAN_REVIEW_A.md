# Plan Review A — reject `--rename` with `--output-format ubam` (#408)

**Reviewer:** A (independent)
**Plan:** `plans/08092026_rename-ubam-rejection/PLAN.md`
**Tree reviewed:** `dev` @ `cbb5ecd` (plan states base `d0bb76b`, which is `HEAD~1`)
**Method:** every factual claim re-derived from source; binary built from a `git archive HEAD` copy under `$TMPDIR` with a private `CARGO_TARGET_DIR`; Perl v0.6.11 read from the `0.6.11` tag.

**Verdict:** the decision (reject, not repair) survives review — the Perl-parity evidence it rests on is **correct**. But three things must change before implementation: an assumption the plan deferred (**A2**) is **false**, the blanket scope destroys a categorically-lossless path that already has a dedicated regression test, and the remediation text tells users to do something that reproduces the exact loss being refused.

---

## Summary of verification

| # | Claim | Verdict |
|---|---|---|
| 1 | Reproduction (FASTQ keeps annotation, uBAM QNAME drops it; whitespace-free survives) | **Confirmed exactly** |
| 2 | Perl v0.6.11 appends `:clip5:` to the end of the whole ID line | **Confirmed** (with one branch nuance + 6 under-cited call sites) |
| 3 | `append_to_id` is shared by both formats via `trim_read` | **Confirmed** — plus a second locus the plan misses (`specialty.rs`) |
| 4 | §3.4a is the right home; ordering within the block immaterial | Home ✓, style ✓, in-block ordering ✓ — **but the block as a whole masks the `--clump_only` message** |
| 5 | A2: nothing asserts `--rename` + uBAM succeeds | **FALSE** — `tests/integration_ubam_out.rs:394` asserts exactly that |
| 6 | Docs errors at `outputs.md:105` and the compatibility list | **Confirmed**, and one more gap found (`--basename` documented nowhere in that guide) |
| 7 | Over-rejection is an acceptable cost | **Mis-framed** — the uBAM→uBAM case is lossless *by spec*, not by luck |
| — | The loss is "silent" | **FALSE at this base** — #406 already emits a `NOTE` naming the dropped text |

---

## 1. Logic review

### 1.1 The reproduction is exact (claim 1 ✓)

Built at `cbb5ecd`, fixture with one space-bearing and one whitespace-free header:

```
--rename --clip_R1 3 -a AGATCGGAAGAGC --dont_gzip     → @withspace 1:N:0:ACGTAC:clip5:ACG
                                                        @nospacehere:clip5:ACG
--rename --clip_R1 3 -a AGATCGGAAGAGC --output-format ubam
                                                      → QNAME [withspace]                    ← annotation gone
                                                        QNAME [nospacehere:clip5:ACG]        ← annotation kept
```

Both halves of the plan's table reproduce, and the whitespace-free contrast that the "data-dependent" claim rests on holds. The mechanism is as described: `append_to_id` (`src/fastq.rs:173-181`) splices before the first TAB but has no SPACE equivalent, so the suffix lands after the description; `parse_name_and_data` (`src/bam.rs:709-723`) then splits on the first space and discards the remainder.

### 1.2 The Perl-parity argument is correct (claim 2 ✓ — the decision's premise holds)

This was the claim the whole decision rests on, so I traced each variable to its assignment rather than reading the append in isolation. All three call-site families append to a **raw, unsplit line**:

| Perl site | Variable | Assignment | Description retained at append? |
|---|---|---|---|
| `1098`, `1115` (RRBS+quality branch) | `$l1` | `my $l1 = <TRIM>;` (`:993`) | **Yes** — raw cutadapt line |
| `1282`, `1298` (poly-A branch) | `$l1` | `my $l1 = <TRIM>;` (`:1219`) | **No** — see nuance below |
| `1399`, `1414` (default branch) | `$l1` | `my $l1 = <TRIM>;` (`:1372`) | **Yes** — raw cutadapt line |
| `1910`, `2000` (`--hardtrim5/3`) | `$identifier` | `my $identifier = <$in>;` (`:1893`) | **Yes** — raw input line |
| `2209`, `2226`, `2243`, `2260` (paired clip) | `$id_1`/`$id_2` | `my $id_1 = <IN1>;` (`:2161`) | **Yes** — raw input line |

Nothing splits on whitespace before the append. The only preprocessing is `$l1 =~ s/\r|\n//g` (line-ending strip), which is exactly what our `id.trim_end()` does. So `@read 1:N:0:IDX` + `--clip_R1 3` under Perl yields `@read 1:N:0:IDX:clip5:ACG` — identical to ours. **Splicing before the space would diverge from shipped Perl. The plan's deciding evidence is sound and the maintainer's choice rests on a true premise.**

**Nuance worth recording (does not change the conclusion).** The poly-A branch is a genuine counterexample: `trim_galore:1255` does `$l1 =~ s/\s+/_/g;` *before* the clip append, so on that path the header is already whitespace-free and `:clip5:` does land on what is effectively a read name. It gets there by mangling the description into underscores rather than by splitting it off, so it is not a precedent for splicing — but it is the one place where Perl's behaviour matches the `readname:clip5:GATC` intent its own warning at `:3456` advertises. If the plan's Context keeps the "irony" aside, this branch is the honest completion of it.

**Under-citation.** The plan cites 5 lines; there are **12**. The missing ones matter for a reason beyond completeness: `1910`/`2000` are the `--hardtrim5`/`--hardtrim3` sites, and that whole mode family is absent from the plan (see §1.4).

### 1.3 `append_to_id` is format-blind (claim 3 ✓), but has a second caller the plan misses

`trim_read` (`src/trimmer.rs:89`) takes `(&mut FastqRecord, &TrimConfig, bool)` — no output-format parameter — and calls `append_to_id` at `:264` / `:273`. It is reached from both the FASTQ writers and the uBAM writers (`src/parallel.rs:463`, `:464`, `:1027`; `src/trimmer.rs:479`, `:480`, `:730`, `:731`). So no fix inside it could be format-scoped without threading the format down. **Confirmed.**

However, `append_to_id` has **six further callers outside `trim_read`**, all in `src/specialty.rs`: `:62`, `:106`, `:166`, `:219` (the `--hardtrim5`/`--hardtrim3` clip annotations) and `:312`/`:313`, `:413`/`:414` (clock/implicon UMI). The plan's Context reasons only about `trim_read`, which understates the surface.

### 1.4 `--hardtrim5` + `--rename` + uBAM is a second, live instance of the same bug — unmentioned

`--rename`'s own help text (`src/cli.rs:300`) says it applies to `--hardtrim5/3`, and those modes are **not** in §3.4a. Verified live:

```
--hardtrim5 10 --rename --output-format ubam   → exit 0, writes in.10bp_5prime.bam
  NOTE: ... First dropped: "1:N:0:ACGTAC:clip5:TTTGGGGGGGGGGCCCCCCCCCC"
  QNAME [withspace]                              ← annotation gone
  QNAME [nospacehere:clip5:TTTG…]                ← annotation kept
```

Good news: the plan's guard is a bare `if self.rename` inside the `matches!(self.output_format, UBam)` block, so it covers the hardtrim path **correctly by construction**. But the plan never says so, its reproduction does not exercise it, and its Validation table has no row for it. That is luck rather than design, and a reader auditing coverage later cannot tell which it was.

### 1.5 §3.4a is the right home and the style matches (claim 4, mostly ✓)

`Cli::validate` opens with the §3.4a block at `src/cli.rs:587-632`: seven arms (`--dont_gzip` `:592`, `--clumpify` `:598`, `--passthrough` `:603`, `--clock` `:606`, `--implicon` `:612`, `--demux` `:617`, `--retain_unpaired` `:625`), each a bare `if self.<field> { anyhow::bail!(…) }`. The proposed code is a byte-for-byte stylistic match, including the multi-line string continuation and the parenthesised reason-plus-remedy shape of the `--retain_unpaired` arm it sits beside. **A1 confirmed:** `--rename` appears in `cli.rs` only at `:199` (a doc comment), `:302-303` (the field), `:851-855` (the `--clump_only` rejection) and `:974` (unrelated prose). No existing §3.4a arm.

**Ordering within §3.4a is immaterial as claimed** — every arm bails, none is a precondition for another. But that is the wrong scope for the question (see §1.6).

### 1.6 The new arm silently changes the `--clump_only --rename` message — Behavior #3 is wrong

The plan states (Behavior 3) that "the `--clump_only` rejection stays as-is and keeps its own wording." That holds only for FASTQ output. §3.4a is at `:587`; the `--clump_only` `--rename` bail is at `:852`. So with `--output-format ubam` the new arm fires **first**.

Current behaviour, verified:

```
--clump_only --rename --output-format ubam
  → Error: --clump_only preserves record contents byte-identically; --rename would mutate read IDs
```

After the change this becomes the new §3.4a message, whose stated reason (*"the annotation is appended after any header description, and BAM read names cannot contain whitespace"*) is **not the reason `--clump_only` refuses `--rename`** — `--clump_only` refuses it because the mode is contractually byte-identity-preserving and does no clipping at all, so there is no `:clip5:` to place anywhere. The user is handed a mechanical explanation for a contractual refusal.

The existing test `tests/integration_clump_only.rs:330 rejects_rename_flag` does **not** pass `--output-format ubam`, so it stays green and will not catch this. Fix options, cheapest first: gate the new arm with `&& !self.clump_only`, or hoist the `--clump_only` `--rename` check above §3.4a. Either way Behavior #3 needs rewording.

### 1.7 The "silent" framing is false at this base — #406 already discloses the loss

The plan says "silently discards" (Goal), "silent, data-dependent loss" (Goal), "converts silent loss into a refusal" (Context), "it closes silent data loss" (step 5), and validation 5 predicts the unpatched build "*succeeds* there and **silently** drops the annotation."

At `cbb5ecd` the run prints, on stderr, once per file:

```
NOTE: BAM read names cannot contain whitespace, so FASTQ header text after the first
space is not carried into uBAM output. First dropped: "1:N:0:ACGTAC:clip5:ACG"
```

The dropped text **includes the `:clip5:ACG`**. This is #406's own fix — and `CHANGELOG.md:32-42`, in the same `### Unreleased / #### Bug fixes` section this plan proposes to append to, documents it as "*A one-time `NOTE:` now reports it and echoes the text that was dropped, so what is lost is visible rather than inferred.*" The plan's base commit `d0bb76b` **is** that fix.

This does not overturn the decision — a stderr NOTE in a pipeline log is easy to miss, it fires once per file rather than per record, and it does not say "the flag you explicitly passed has been neutralised." That is still a good reason to refuse. But a CHANGELOG bullet claiming to close *silent* data loss would contradict the bullet ten lines above it. The honest framing is narrower and still compelling: *the flag the user explicitly requested is silently neutralised; only the generic description-dropping NOTE hints at it, and only for the first record.*

### 1.8 The remediation advice does not work

The proposed message ends "*use FASTQ output, or convert afterwards*". The second half is wrong. `samtools import` splits QNAME on whitespace exactly as we do:

```
FASTQ from --rename:  @withspace 1:N:0:ACGTAC:clip5:ACG
samtools import       → QNAME withspace                    ← same loss
samtools import -T '*' → QNAME withspace                   ← same loss
```

A user who follows the advice reproduces the refused failure with an extra step. "Use FASTQ output" is sound; "convert afterwards" must either be dropped or qualified (the annotation survives conversion only if the description is stripped or folded into the name first).

---

## 2. Assumptions

**A1 — `--rename` absent from §3.4a.** ✓ Verified independently (§1.5).

**A2 — "no test or CI grep asserts that `--rename` + uBAM **succeeds**", deferred to implementation.** ✗ **FALSE, and the deferral was not safe.** `tests/integration_ubam_out.rs:393-440`:

```rust
fn ubam_out_rename_with_preserve_tags_keeps_tags_intact() {
    // Code-review C1 regression guard: --rename + --preserve-tags must
    // NOT corrupt the last preserved tag value. …
    …args(["--clip_R1","5","--rename","--output-format","ubam",
           "--preserve-tags","CB,UB"])
      .arg("test_files/ubam_test_with_tags.bam") …
    assert!(status.success(), "trim_galore exited non-zero");   // ← line 420
```

Line 420 asserts precisely what A2 says nothing asserts. The fixture guard at `:400` is not a get-out: `test_files/ubam_test_with_tags.bam` exists (587 bytes, committed 27 Jun), so the test runs rather than skipping. **`cargo test` will fail on the first implementation attempt.** The plan's step 6 runs the suite, so this would be caught — but as a mystery failure during implementation rather than a planned edit, which is the cost the #400 lesson was supposed to prevent. It also means the implementation outline is missing a step.

This matters for more than bookkeeping: the test is the **only** regression guard for the C1 tab-splice in `append_to_id`, and it is the thread that leads to §2's real problem.

**A3 — `--rename` not in the byte-identity matrix.** ✓ Verified. `.github/` has exactly two `rename` hits, both unrelated (`ci.yml:8` a branch-rename comment, `ci.yml:1002` a file rename inside a diff step).

**A4 — FASTQ-path behaviour unchanged, Perl parity untouched by construction.** ✓ Sound. A CLI-validate bail cannot reach `append_to_id`.

**Unstated assumption, and it is the load-bearing one:** that every `--rename` + uBAM-output run is at risk. It is not — see §3.1.

**Base-commit drift:** the plan says base `d0bb76b`; the tree is `cbb5ecd` (`#397`, collision refusals). `d0bb76b` is `HEAD~1`. No conflict with this change, but the header should say `cbb5ecd`.

---

## 3. The scope problem

### 3.1 uBAM-in → uBAM-out is lossless *by specification*, not by luck

The plan's Self-Review accepts over-rejection on these grounds: "whitespace-free headers lose nothing today, so this refusal is stricter than strictly necessary for them. Accepted deliberately: the alternative is a refusal whose firing depends on the shape of the first record, which is worse than a predictable one."

That dichotomy is false, because it treats "whitespace-free" as a property of the *data*. For uBAM input it is a property of the *format*: **SAM forbids whitespace in QNAME**, as our own code says at `src/bam.rs:703` ("The BAM spec rejects whitespace in QNAME"). A `FastqRecord` synthesised from a BAM record (`src/bam.rs:898`) is `@QNAME` plus an optional TAB-delimited tag tail — a space is unreachable. And the TAB case is exactly what the C1 splice already handles.

So uBAM→uBAM with `--rename` cannot lose the annotation. Verified end-to-end:

```
--clip_R1 5 --rename --output-format ubam --preserve-tags CB,UB  test_files/ubam_test_with_tags.bam

in :  SRR24827378.1                 CB:Z:ATCGATCG-1  UB:Z:GCTAGCTA
out:  SRR24827378.1:clip5:AATTA     CB:Z:ATCGATCG-1  UB:Z:GCTAGCTA     ← annotation AND tags intact
      SRR24827378.2:clip5:ACGTA     CB:Z:ATCGATCG-1  UB:Z:GCTAGCTA
(zero "not carried into uBAM" NOTEs emitted)
```

The plan would refuse this. It is a working, tested, spec-guaranteed-lossless path — single-cell uBAM carrying `CB`/`UB` through a UMI-aware clip is precisely the workflow `--rename` exists for (`cli.rs:301`: "Useful for UMI handling"), and it is the workflow the C1 fix was written to protect.

**To be explicit: this is not a request to revisit reject-versus-repair.** The maintainer's decision stands and no code inside `append_to_id` need change. The question is only *which inputs* the rejection covers.

### 3.2 The precise guard is cheap and is deterministic

Reject `--rename` + uBAM output when **at least one input is FASTQ**. That is decided by input file format — known before any trimming, independent of record shape — so it is every bit as predictable as a blanket refusal.

It cannot live in `Cli::validate` (no format detection there), but `main.rs` already has the pattern *and* the data. `src/main.rs:247-254` computes `input_formats: Vec<InputFormat>` once, and §3.4b at `:298-315` is an existing format-detection-time rule of exactly this shape (`--preserve-tags` + all-FASTQ → hard error under uBAM output, warning otherwise). A sibling arm right below it needs `input_formats.iter().all(…)` and a `bail!`.

Trade-offs, stated plainly:

| | Blanket, in `cli.rs` §3.4a (as planned) | FASTQ-input-scoped, in `main.rs` §3.4b |
|---|---|---|
| Refuses the lossy FASTQ→uBAM direction | yes | yes |
| Refuses the lossless uBAM→uBAM direction | yes (regression) | no |
| Fires before any I/O | yes | yes (§3.4b runs before `ensure_output_dir`) |
| Predictable | yes | yes — keyed on file format, not record shape |
| Existing C1 test | **must be deleted or rewritten** | **stays green unchanged** |
| Cost | 6 lines, one file | ~8 lines, one file, one extra `all()` |
| Mixed FASTQ+uBAM inputs (legal per `CHANGELOG:25`) | refused | refused (`!all_bam`) |

The blanket version is simpler to explain in one sentence, which has real value. But it buys that simplicity by deleting a regression guard for a fix the project made deliberately three commits ago, and by refusing a workflow that provably works. My recommendation is the scoped guard.

**In fairness, there is an in-repo precedent that argues the other way, and it should be weighed.** The `--phred64` + uBAM guard (#358) faced the same choice and chose uniformity, in a comment at `src/main.rs:264-269`:

> Rejected uniformly rather than per-mode. `--hardtrim5/3` and `--clump_only` currently accept the flag as an inert no-op, so this does remove working invocations — but a mode-dependent rule ("rejected unless you're in hardtrim, or clump-only-to-FASTQ, …") is worse to document, worse to test, and one refactor away from being wrong.

That is a direct answer to my recommendation, from this codebase, about this feature area — including the same willingness to "remove working invocations." If Felix applies the same reasoning here, the blanket guard is the consistent choice and §3.3 becomes the required work rather than a fallback.

Two things distinguish the cases, and I think they are enough to tip it the other way, but not decisively:

- **#358's rejected alternative was *mode*-dependent** (hardtrim vs clump-only vs trim — an open-ended list that grows with every new mode). Mine is *input-format*-dependent: a closed, two-valued distinction that cannot grow. `input_formats` is already computed for exactly this kind of decision, and §3.4b at `:298-315` is a format-dependent rule the same codebase accepted without reservation.
- **#358's uniformity cost nothing that worked correctly** — the flag was an "inert no-op" on the paths it newly refused, so nobody lost output they wanted. Here the refused path produces *correct, tag-preserving output today* and has a test asserting it. That is a materially larger cost than removing a no-op.

The same comment also confirms the mechanics of my suggestion: `:271-273` says the #358 guard is "sited here, not in `Cli::validate()` §3.4a: validate() runs before input format detection and cannot see `any_bam`. Same reason §3.4b below lives in main.rs." So `main.rs` is the established home for precisely this shape of rule.

### 3.3 If the blanket rejection is kept anyway

Then two consequences must be handled explicitly rather than absorbed:

1. `ubam_out_rename_with_preserve_tags_keeps_tags_intact` must become a **rejection** test, and the C1 tab-splice guard must be **re-homed to FASTQ output**, not lost. That route still exercises the splice — verified:
   ```
   --clip_R1 5 --rename --preserve-tags CB,UB --dont_gzip   test_files/ubam_test_with_tags.bam
   → @SRR24827378.1:clip5:AATTA<TAB>CB:Z:ATCGATCG-1<TAB>UB:Z:GCTAGCTA
   ```
   So the splice stays reachable and testable after the rejection; without this edit it becomes untested behaviour that a future refactor can quietly break.
2. The CHANGELOG bullet should say the lossless uBAM→uBAM case is being withdrawn too, and why (one rule beats two). Users of that path get an error where they previously got correct output, and it should not surprise them.

---

## 4. Efficiency / integration

Nil either way, as the plan says: one branch on the argument-parsing path (or one `all()` over an already-materialised `Vec<InputFormat>`). No I/O, no allocation of consequence. Nothing to add.

---

## 5. Validation sufficiency

The table's shape is good — validation 5 as an expected-fail control is the right instinct, and validation 2 guarding against over-reach is the right worry. Gaps:

1. **No row for the existing conflicting test** (A2's real answer). Needs an explicit row: the C1 guard is retargeted, and the tab-splice remains covered somewhere.
2. **No row for `--hardtrim5`/`--hardtrim3`** — a live instance of the identical loss (§1.4), covered only incidentally.
3. **No row for uBAM input.** With the scoped guard this is the central positive case ("uBAM→uBAM still succeeds and still carries `:clip5:` plus tags"). With the blanket guard it is the case whose withdrawal needs asserting. Either way it must be tested; today the table cannot distinguish the two designs.
4. **No paired row.** `--rename` annotates R1 and R2 independently (`trimmer.rs:479-480`, Perl `:2209`/`:2226`), and paired uBAM output is one interleaved BAM. `Cli::validate` runs before the paired split so the guard is format-agnostic, but nothing asserts it.
5. **Validation 5's prediction is wrong as written** — the unpatched build does not fail *silently*; it prints the #406 NOTE (§1.7). Correct the expectation or the control will look like it found a discrepancy.
6. **Validation 6** should also confirm `--basename` is documented in `outputs.md` after the edit (§6.2), since the plan's own step 3 leaves that conditional.

---

## 6. Docs (claim 6 ✓, plus one more gap)

### 6.1 Both reported errors are real

`docs/src/content/docs/guide/outputs.md:105` — "`--rename PREFIX` replaces the input filename stem in the output names. Useful for pipelines that thread sample IDs through trimming separately from input filenames." Wrong in every clause: `--rename` is `pub rename: bool` (`cli.rs:302-303`), takes no value, and never touches filenames. The described behaviour is `--basename`'s (`changelog.md:962`). **Confirmed.**

The §Feature compatibility list is missing `--rename`. **Confirmed.** Minor correction: the plan cites `:88-95`; `:88` is the heading and the bullets run `:90-97` (seven, not six — the plan's own count is off by the `--clumpify` bullet).

### 6.2 Nothing else in that list is stale — but `--basename` is missing entirely

I checked each of the seven bullets against §3.4a: `--dont_gzip` (`:592`), `--clock`/`--implicon` (`:606`/`:612`), `--demux` (`:617`), `--passthrough` (`:603`), `--retain_unpaired` (`:625`), `--clumpify` (`:598`). All seven present arms are documented, and the `--clumpify` bullet's `--clump_only` parenthetical is accurate. **The plan is right that `--rename` is the only gap.**

The plan's step 3 hedges: "Check whether the surrounding section already documents `--basename`; if not, one clause." Resolved: **`--basename` does not appear anywhere in `outputs.md`.** It is documented only in `modes/clump-only.md:79`, `modes/passthrough.md:28`, and the changelog. So the guide's `## Renaming outputs` section — the one place a reader looks for output naming — currently describes a flag that does not exist and omits the flag that does the job.

That makes the section heading itself part of the problem. `## Renaming outputs` is about output *filenames*; `--rename` is a read-ID annotation and does not belong under it at any length. Cleanest resolution: make that section about `--basename` (which genuinely renames outputs), and document `--rename` where read-ID mutation belongs — or as a clearly separate subsection with a disambiguating note, since `--rename` / `--basename` is an easy pair to confuse and the current text proves it. The plan's step 3 as written ("replace the sentence with what `--rename` actually does") would leave a read-ID feature sitting under a filename heading.

---

## 7. Alternatives

1. **FASTQ-input-scoped rejection in `main.rs` beside §3.4b — recommended.** §3.2. Same refusal for every lossy case, keeps the spec-lossless one, keeps the C1 test green unchanged.
2. **Blanket rejection in §3.4a — as planned.** Simplest rule, and backed by the #358 precedent at `main.rs:264-269` (§3.2). Costs the uBAM→uBAM path and forces the C1 guard to be re-homed (§3.3). A legitimate choice if made knowingly; the plan currently makes it without knowing the path is spec-safe.
3. **Warn instead of reject.** Rejected — already effectively the status quo via #406's NOTE, and §1.7 shows why that is not enough: the NOTE describes description-dropping generically and never says the requested flag was neutralised.
4. **Splice before space in `append_to_id`.** Correctly ruled out; §1.2 confirms it would diverge from shipped Perl on ten of twelve call sites.
5. **Move the annotation into a BAM aux tag** (e.g. `XC:Z:ACG`) instead of the QNAME. Out of scope for #408 and a real design question (tag-namespace choice, `--preserve-tags` interaction), but it is the only option that would let `--rename` work for FASTQ→uBAM at all. Worth an issue rather than silence, since the rejection message will prompt users to ask for it.

---

## Action items

### Critical

- **A2 is false — plan an edit to `tests/integration_ubam_out.rs:393-440`.** It asserts `status.success()` on `--clip_R1 5 --rename --output-format ubam --preserve-tags CB,UB` with a committed fixture, so it runs and will fail. Add an explicit implementation step; do not leave this to be discovered by `cargo test`. Correct A2's text — the deferral was not safe.
- **Decide the guard's scope deliberately, and record the decision.** uBAM→uBAM with `--rename` is lossless by SAM spec (`bam.rs:703`) and verified lossless end-to-end including aux tags (§3.1) — the blanket form refuses a working, tested, single-cell-relevant path, which the plan's Self-Review does not currently account for (it treats the safe case as data-dependent luck). I recommend narrowing to "at least one FASTQ input" in `main.rs` beside §3.4b. **But note the counter-precedent at `src/main.rs:264-269`**, where #358 chose uniform rejection over a conditional rule in this same feature area and explicitly accepted removing working invocations; if Felix follows it, that is defensible and the plan should cite it rather than leave the trade unexamined. Either way: if the blanket form is kept, re-home the C1 tab-splice guard to the FASTQ-output route (§3.3, verified to still exercise it) and disclose the withdrawal of the lossless path in the CHANGELOG.
- **Fix the remediation text.** "convert afterwards" reproduces the identical loss — `samtools import`, with and without `-T '*'`, drops the description exactly as we do (§1.8). Drop that clause or qualify it.

### Important

- **Stop calling the loss "silent."** #406 already prints a `NOTE` naming the dropped text, and `CHANGELOG.md:32-42` documents that in the very section this plan appends to (§1.7). Reframe to what is actually silent: the explicitly-requested flag is neutralised with no flag-specific diagnostic. Update Goal, Context, step 5, and validation 5's prediction.
- **Correct Behavior #3 and decide the `--clump_only` overlap.** §3.4a precedes `cli.rs:852`, so `--clump_only --rename --output-format ubam` will switch from the byte-identity message to the new mechanical one, which is the wrong explanation for that mode. Gate with `&& !self.clump_only` or hoist the clump_only check (§1.6). The existing test will not catch this.
- **Acknowledge the `--hardtrim5`/`--hardtrim3` path.** Same loss, live today, covered by the guard only incidentally (§1.4). Say it is intended and add a validation row.
- **`--basename` is absent from `outputs.md` entirely**, and `## Renaming outputs` is a filename heading. Resolve the section's structure rather than swapping one sentence into it (§6.2).

### Optional

- Add the six under-cited Perl lines (`1115`, `1298`, `1414`, `1910`, `2000`, `2243`, `2260`) — `1910`/`2000` are the hardtrim sites and connect to the gap above.
- Record the poly-A branch nuance (`trim_galore:1255` folds whitespace to `_` before appending), which completes the "irony" aside honestly (§1.2).
- Note that `append_to_id` has six callers in `specialty.rs` beyond `trim_read` (§1.3).
- Fix line references: base `cbb5ecd` not `d0bb76b`; compatibility bullets `:90-97` not `:88-95` (seven items).
- Add validation rows for paired-end and for uBAM input (§5).
- Consider filing a follow-up for carrying the annotation as a BAM aux tag (alternative 5) — the rejection message will invite the request.
