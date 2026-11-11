# Plan review — `phase3/PLAN.md` v2 (defend and correct the `--clumpify` memory promise, #439)

**Reviewer A**, fresh context, independent. Target: `plans/08162026_clumpify-memory-accounting/phase3/PLAN.md` (v2).
Repo `/Users/fkrueger/Github/TrimGalore`, branch `dev`, HEAD **`fdade7e`** (verified `git rev-parse`, equal to the
plan's stated baseline). Read-only: nothing in `src/`, `tests/`, `docs/` or `.github/` edited, no branch switched.

> ### Scope and provenance — read first
>
> **I am not the author of `phase3/PLAN_REVIEW_A.md`.** That file (42,525 B, 17 Aug 18:34) is a different
> Reviewer A lineage. I was briefed on the **parent** `../PLAN.md` at tip `2ab0918` and delivered
> `../PLAN_REVIEW_A_v2.md` (50,656 B, 20:59). This document is my first reading of `phase3/PLAN.md`,
> written at the team lead's explicit request.
>
> **I ran zero `trim_galore` invocations.** There is therefore no fixture path, rep interleaving, free
> disk, or exit code to report for "my numbers" — I have none. Everything below is one of:
> - **[V-code]** a `git`/`grep`/`sed`/`Read` observation of source at `fdade7e`, command quoted;
> - **[V-arith]** a Python reimplementation of `clumpify_min_memory_bytes` / `resolve_layout` from the
>   constants as read at `fdade7e`, compared against the plan's published figures;
> - **[NV]** not verified by me — every RSS measurement in this plan falls here (the ramp, the
>   perturbation matrix, the Linux extrapolation, the `threading.md` tables, all of §0).
>
> **§0 is withdrawn per the team lead** (their harness ran with 1.6 GiB free while writing 2.9 GB/rep;
> a controlled re-run gives 756–789 MiB at `1G/c8` plain output with 9 GiB free, against the 1081.8 MiB
> §0 claims). I spent no effort re-measuring it. §1 below is entirely about what that withdrawal does
> to the rest of the document, which is the largest single finding in this review.

---

> ## ⚠ STATUS — this review targets v2; **v3 landed at 21:08:38 and resolves both Criticals**
>
> I read `phase3/PLAN.md` at 19:06 (41,678 B, titled v2) and wrote the review below. It was rewritten
> to **v3** (36,637 B) while I was writing. v3 independently reaches the same two conclusions:
>
> | my finding | v3's disposition |
> |---|---|
> | **Critical 1** — sever §0's twelve dependents | **Done.** "Retracts §0"; the queued issue is narrowed by **D6** to the false deprecation message plus the `clump-only.md:78` asymmetry, explicitly "**No memory claim**", citing #245 and `io.rs:14-32`. |
> | **Critical 2** — don't gate on an unmeasured Linux baseline | **Done, and better than my proposal.** **D5**: a baseline-relative band on a *committed* ratio plus the absolute promise assertion (§5.3). My staging suggestion is a cruder version of the same idea. |
> | **§1 process finding** — assert free disk per rep | **Done.** PR 3 step 3 "Assert free disk before the first cell (§5.4)", plus Validation 14: a cell run under artificially low free disk must exit non-zero "and the message names disk, not memory". |
> | multi-pair disposition | **Partly done** — v3 adds the double-banner parse fix ("pinning the first match") and gives the cell its own Validation 11 row. |
>
> **One claim in my §1 was wrong, and v3 had already caught it.** I listed §0's three CLI behaviours as
> surviving the withdrawal intact. Two of the three do not. Specifically **"there is no way to obtain
> gzipped output from plain input" is false**: `main.rs:424` derives output encoding from
> `is_gzipped(&cli.input[0])`, i.e. the input *filename*, while the reader sniffs bytes — so naming a
> plain file `*.gz` (or symlinking it) yields gzipped output. v3 ships exactly that as the documented
> workaround under D6. The claim that "always gzip" was never tried is also wrong: it shipped in
> v2.1.0-beta.5 and was reverted under **#245**. Only the third behaviour survives, and it is the one I
> actually verified: `--clumpify` refuses `--dont_gzip` as pointless (`cli.rs:747-751`) and then clumps
> plain text anyway on plain input. **§1's "Survives" list is corrected in place below.**
>
> **Findings below that v3 does *not* address — the live remainder of this review:**
> Important **3** (the 25 bp cell still cannot fail — v3 lists short reads as Remaining risk (e) and an
> open question, but adds no band assertion), **4** (`2G/c4`'s known residual now has no guard and no
> tracking issue), **6** (Validation 12 still mutates only the helper form; the CI job is the *only*
> guard that catches a re-inlined literal, which is PR 3's strongest justification and is still absent
> from §5 — v3 repeats the concession at line 255 without drawing the consequence), and **7** (pin the
> sum `work + 1 + result_per_core == 4`, not two independent literals). Everything in §6 — the ~20
> verified `file:line` and arithmetic claims — applies unchanged to v3.

---

## Verdict

**Implement, after severing §0 and staging the CI job.** This is the strongest document in the arc. I
checked every `file:line` reference and every published number I could reach without running the binary
— **~20 claims, zero errors** (§6). The guard rework is correct and the algebra behind it is right. Two
things block a clean implementation:

1. **§0's withdrawal has ~12 load-bearing dependents in this document**, including the Bar on the first
   issue in §7 and the stated justification for the gzipped CI fixture. Left unsevered, v3 files an
   issue asserting a data-integrity breach that the lead has retracted. **Critical.**
2. **PR 3 lands a gating assertion on an unmeasured Linux baseline.** Validation 14 forbids tuning the
   threshold if it goes red, and offers no other branch — so the realistic outcome is a disabled job.
   **Critical**, and cheap to fix with a staging step.

Everything else is Important or below.

---

## 1. §0's withdrawal — the cascade, and what survives  **[Critical]**

§0 is the spine of v2: the revision history leads with it, and it is cited by name in D3, §5, §6, §9,
Assumptions 7 and 8, Self-Review ¶2, and Validation 16. Withdrawing the measurement does **not**
withdraw all of it, and the distinction matters, so here is each dependent with a disposition.

**Falls — must be struck from v3:**

| Where | The claim that dies |
|---|---|
| Revision history, v2 bullet | "records a defect found while resolving the one place the reviewers contradicted each other" as the headline |
| §0 lines 25–45 | the whole 2×2 table's plain-output column; "`1G/c8` breached reproducibly"; "4 of 5 over the prediction and 3 of 5 over the budget" |
| §0 line 42 | "**Every breach lives on the plain-text output branch**" |
| §0 "coupled part" 3 (line 74) | "there peak RSS reaches **1.057 × `--memory`** at `--cores 8` … (measured, 3 of 5 reps)" |
| §7, issue 1, Bar column | "**data-integrity**: breaks the `--memory` promise at 1.057× budget" |
| §5, first bullet | "Gzip is not an optimisation here — §0 shows the plain-output branch is a different memory regime" |
| §6, third bullet | "**Lead with the gzip-output figure** … §0 is the reason" |
| §9, second dated note | "the plain-output branch is a different regime" — do not write this into `CALIBRATION.md` |
| Self-Review ¶2 | "The measurement lesson is bigger than either guard" — the entire paragraph |
| Assumption 1's parenthetical | "(Narrowed from v1's broader 'the model is correct as shipped', which §0 shows is true only on the gz path)" |
| Remaining risk (e) | "`--clumpify` on plain input keeps breaching the budget until that issue is picked up" |

**Survives — and v3 should say explicitly that it survives, or a later reader will bin it with §0:**

- **Corrected — only ONE of §0's three CLI behaviours survives, not all three.** (See the status block
  above; v3 had already established this and I had it wrong here.) Behaviour 2, "there is no way to
  obtain gzipped output from plain input", is **false** — `main.rs:424` takes output encoding from
  `is_gzipped(&cli.input[0])`, the input *filename*, while the reader sniffs bytes, so naming a plain
  file `*.gz` yields gzipped output; v3 ships that as D6's documented workaround. The "always gzip"
  alternative was also already tried and reverted (**#245**). What does survive is behaviour 3, which is
  the one I verified directly **[V-code]** — `src/cli.rs:747-751`:
  ```rust
  if self.dont_gzip {
      anyhow::bail!(
          "--clumpify and --dont_gzip are mutually exclusive (clumping plain text is pointless)"
      );
  ```
  So `--clumpify` genuinely refuses to clump plain text as pointless, and genuinely does it anyway when
  the input is plain. That is a coherent **UX/consistency defect** and still worth an issue.
- **The `clump-only.md:78` asymmetry** (`--dont_gzip` documented as supported under `--clump_only`)
  stands — it is a docs-vs-docs contradiction, no measurement involved.
- **Assumption 8's factual half** — "`--clumpify` output is gzipped. False for plain input" — stands.
- **The `.fq.gz` assertion in PR 3 should be KEPT**, with its rationale rewritten. It no longer guards
  "a pathological memory regime"; it guards **measurement validity** — a CI job that silently measured
  the plain-output branch would not be comparing against the gz-path calibration. That is a sufficient
  reason on its own and it is the reason §0's process lesson was right even though its numbers were not.
- **The gzipped fixture is still right**, on the disk argument alone (2.96 GB → ~490 MB), which §5
  already states independently. `timeout-minutes: 30` and the disk sizing are unaffected.

**Consequent re-scoping v3 owes:**

- **D3's issue drops from data-integrity to user-visible-wrong-behaviour.** Against the repo's stated
  filing bar that is still fileable (the `--gzip` deprecation message is *false* for plain input, and
  there is no way to request gzip), but the plan must not file it with the 1.057× claim attached. Note
  this also weakens D3's rationale for keeping Phase 3 free of the behaviour change — the change is now
  smaller and less urgent, which if anything strengthens D3, not weakens it.
- **Validation 16 must split**: the three commands still reproduce (assert that); the breach does not
  (delete that half).
- **§6's threading figures need their branch stated.** The plan says "all three harnesses measured
  plain output". With §0 withdrawn, we no longer know that plain output is a different regime — so the
  §6 tables may be fine as they stand. But they are **[NV]** for me and the plan's own instruction to
  "lead with the gzip-output figure" now has no support. **Re-derive or drop that instruction**; do not
  file §6 with an unexplained branch caveat.

**Process finding, and I think the most valuable thing in this review.** §0 is the **fourth**
measurement retraction in this arc — after v1's coefficient, v1's two guards, and the 4.86%
MiB-vs-`1e6` harness error. **All four were instrumentation defects, not model defects.** The plan's
Validation preamble ("Four false negatives in one recent session came from harnesses that could not
fail — prove each of these can") is exactly the right response, but it is applied to three CI rows only.
It should be a standing requirement on the *measurement* harness, and the specific precondition that
would have caught §0 is trivial: **assert free disk ≥ 3× the expected write volume before each rep, and
record it in the results row.** Add that to the generator/harness spec in §Signature. Every retraction
in this arc would have been caught by a precondition assertion of that shape.

---

## 2. PR 3 gates on an unmeasured baseline  **[Critical]**

§5 expects ~0.85 (`1G/c4`) and ~0.72 (`--fastqc`) on Linux, extrapolated from macOS by −9 points
(Assumption 6: "all ratios here are macOS/arm64"; Remaining risk (a): "the Linux baseline is
unmeasured"). PR 3 step 6 then asserts `peak ≤ prediction` as **gating on the first run**.

Validation 14 handles the bad case correctly in spirit — "if any cell exceeds 1.0 on unmodified code,
stop and treat it as a finding, not a threshold to tune" — but it gives no operational branch. A new
required job that goes red on `dev` for unmodified code creates immediate pressure, and the only
permitted response is "stop", which nobody does at 5pm. The realistic outcome is that the assertion gets
loosened or the job gets `continue-on-error`'d permanently — i.e. the guard this whole phase exists to
build is the thing most likely to be discarded.

Note the withdrawn §0 is itself the evidence for this risk: an *intermittent* breach at `1G/c8` (4 of 5
reps) is precisely the failure shape that would land here, and it took three rounds and a controlled
re-run to attribute. Whatever its cause, the arc has already demonstrated that this measurement is
capable of intermittent reds that are not model drift.

**Recommendation, concrete and one line each:**

1. Land PR 3 with `continue-on-error: true` on the job, all assertions live, both ratios echoed.
2. Collect ≥10 runs spanning a cold and a warm cargo cache (the plan already flags cold-cache cost as
   the reason for `timeout-minutes: 30`).
3. Flip to gating in a follow-up PR that records the observed Linux min/max in the job comment.

This costs one extra PR and converts Validation 14 from an instruction nobody can follow into a
scheduled decision with data. It also gives Assumption 6 a closing date.

---

## 3. The 25 bp cell is a harness that cannot fail  **[Important]**

§8 adds a 25 bp cell because σ ≈ 1.55 measured there against the shipped σ = 23/16 = **1.4375**
**[V-code]** (`src/clump.rs:245-247`) — i.e. a *known* read length at which the model under-reserves,
in the promise-breaking direction. The cell is specified as "logged, non-gating".

So CI will faithfully print a ratio above 1.0 and pass, every run, forever. That is the exact pattern
the plan's own Validation preamble condemns, and §8's justification for adding the cell ("the guard is
structurally unable to see it") is undercut by making the new cell unable to fail either.

**Recommendation:** assert the 25 bp ratio against a **band** recorded from the first runs (e.g.
`≤ 1.15`), so today's known breach is tolerated but a *worsening* is caught. And file the short-read
under-reservation as its own issue now — a known promise breach living only as a logged CI line is, in
reviewer B's phrase from §7, "the shape of thing that rots".

---

## 4. Dropping `2G/c4` removes the only cell with a known residual  **[Important]**

§5 drops `2G/c4` on B's argument (pool 672 MiB, needs ~2.3× the data to saturate, most RSS, sole reason
Assumption 4's RAM half is load-bearing). The reasoning is sound and I would keep the decision.

But §4c states that the parent plan "names it as the next lever someone will pull for the `2G/c4`
residual" — so the one cell exhibiting a known unexplained residual now has **no guard at all**, and the
plan does not record that. **Add it to Remaining risks and give the residual a tracking issue**, or it
gets rediscovered as a novel finding in six months. This is bookkeeping, not a design change.

---

## 5. Two refinements to the guards  **[Important]**

### 5a. The CI job is the *only* guard that catches a re-inlined literal — say so, it is PR 3's best argument

Remaining risk (c) concedes that 4c "catches a helper edit, not a re-inlined literal", and §4's honesty
note says re-inlining "is the path of least resistance" because the literal is what stands there today.
Both true **[V-code]** — `src/parallel.rs:117` and `:943` are today
`let queue_depth = if clump_layout.is_some() { 1 } else { 2 };`, i.e. inline conditionals, not a helper
call.

The plan never connects this to PR 3. It should: **the CI job measures real RSS, so it catches a
re-inlined literal that 4c structurally cannot.** That is the strongest available justification for
spending a runner on this, and it is currently missing from §5 while §5 instead argues from fixture
saturation. Also update **Validation 12** to run the mutation *both* ways — as a helper edit (proves 4c)
and as a re-inlined literal at `:117` (proves the CI job earns its keep). As written, Validation 12
tests only the form 4c already covers.

### 5b. Pin the sum, not two independent literals

I verified `k`'s decomposition **[V-code]** (`src/clump.rs:246-248`):

> *"Per worker: `k × B` for the queued batch, the batch in hand, and the compressed output waiting on
> the result channel, k = 4."*

Against the shipped depths: work channel **1** (`:117`, `:943`), plus 1 in hand, plus result channel
`cores * 2` → **2 per core** (`:130`, `:959`). So `1 + 1 + 2 = 4 = k`. The identity is exact, and the
plan's §Signature (`work: 1, result_per_core: 2` for clumpy; `work: 2` for plain) **matches shipped
behaviour precisely** — so 4c is a genuine pure refactor with no behaviour change, and the plan's
Integration section ("five changed stderr lines … no output-file bytes change") is correct. I had
expected to find a silent regime change here and did not; credit where due.

**Refinement:** have the pin assert `work + 1 + result_per_core == 4` with a failure message naming
`k` and `dyn_denominator`, in addition to pinning the two values. Two independent literal pins do not
communicate that their *sum* is what `k` encodes, which is the invariant a future editor needs to see.
The plan's own comment rule wants the invariant on the thing it constrains.

---

## 6. Claims I verified — a clean sweep

Recorded compactly because the density of correct references here is unusual and v3 should not re-derive
them.

**[V-arith]** Reimplemented `clumpify_min_memory_bytes` from `STATIC_TRIM_BYTES = 224 MiB`,
`STATIC_FASTQC_PER_THREAD_BYTES = 24 MiB`, `FASTQC_CHARGED_THREAD_CAP = 16`, `MIN_BIN_BYTES = 1 MiB`,
`MARGIN 10/11`, `dyn_denominator = 23·n_bins + 64·workers`, `n_bins = max(16, 4·bin_cores)`:

| plan figure | mine | |
|---|---|---|
| floor `c2` plain 281 | 280.50 → ⌈281⌉ | ✔ |
| floor `c6` plain 311 | 310.75 → ⌈311⌉ | ✔ |
| floor `c32` plain 590 | 589.60 → ⌈590⌉ | ✔ |
| floor `c2` fastqc 334 | 333.30 → ⌈334⌉ | ✔ |
| floor `c32` fastqc 1012 | 1012.00 | ✔ |
| floor `c64` 933 / 1356 | 932.80 / 1355.20 | ✔ |
| `clump_only c1` 277 | 276.10 → ⌈277⌉ | ✔ |
| `clump_only` table 329 / 513 / 579 / 611 / 775 / 876 | 328.90 / 512.x / 578.05 / 610.78 / 774.40 / 875.60 | ✔ all six |
| `clump_only` plain 352 (`c16`) / 454 (`c32`) | 352.00 / 453.20 | ✔ |
| crossover exactly `--cores ≥ 11` | `c10` 578 < 580 (no); `c11` 610 > 586 (yes) | ✔ |
| widest window 154.4 MiB at `c16` | 774.40 − 620 = 154.40 | ✔ |
| `1G/c8` → 32 bins × 9 MiB | 9.06 → 9 | ✔ |
| `4G/c8` → 32 bins × 44 MiB | 44.87 → 44 | ✔ |
| banner `≈ 930 of 1024`; `≈ 1363 of 1500` | 930.9 → 930; 1363.6 → 1363 | ✔ (Validation 1 correct) |

The `div_ceil`-to-MiB convention the plan states is the one that reproduces its numbers; truncation does
not. v1's two corrected cells (277-vs-329 wrong series; 578-vs-579 truncation) are both real errors and
both correctly fixed.

**[V-code]** All of these check out at `fdade7e`:

- `src/main.rs:207-213` — the banner, currently `"clumpify: {} bins × {} MiB; predicted peak ≈ {} MiB (gzip level {})"`. §1's guard-consequence analysis is **correct**: after D2 the peak loses its own unit token, so the tripwire's `"MiB (gzip level"` clause pins the *budget's* divisor. Changing `predicted_peak_bytes()/(1024*1024)` → `/(1000*1000)` gives `≈ 975 of 1024 MiB budget (gzip level 1)` and both existing clauses still match.
- `tests/integration_preflight_tripwire.rs:597-603` — two `contains` clauses exactly as described.
- `src/clump_only.rs:275 / 421 / 764 / 929` — the four `resolve_layout` sites.
- `src/clump_only.rs:289 / 436 / 784 / 1031` — Q4's terminology collision is real: `:289`, `:784`, `:1031` say "bins × N MiB **budget**", `:436` omits "budget". Five stderr lines total with D2. ✔
- `src/parallel.rs:117, :130, :943, :959` — all four depth sites, as §5's file table lists them.
- `src/parallel.rs:33` — `const BATCH_SIZE: usize = 4096`, and its docstring does carry the read-length-blind "~300 bytes" assumption §6 cites.
- `ClumpLayout::new` (`src/clump.rs:150`) is `pub` with **exactly one caller**, `src/parallel.rs:1615`, in test code. `pub(crate)` is therefore safe. One nuance the plan may want: the crate is a library as well as a binary, so narrowing a `pub` fn is semver-breaking in principle — irrelevant here given no production caller, but do not let the CHANGELOG claim "no API change" if it says that anywhere.
- §7's EOF-send asymmetry is real and the line numbers are right: `let _ = …send(…)` at `:714`, `:754`, `:1202` against `.is_err()` at `:794`, `:825`, `:1193`, `:1228`. **Two details worth adding to the issue:** `:1193` (checked) and `:1202` (unchecked) are nine lines apart in the same function, which reads as oversight rather than intent; and the neighbouring comment cites **#434**, so this is a gap in #434's own fix. Both strengthen the filing.
- Assumption 3 — `--clumpify` requires `--cores ≥ 2` — confirmed at `src/cli.rs:743-745`.

**§4's algebra is right and I reached the same conclusion independently** on the parent plan:
`b = ⌊usable·16/d⌋ ⇒ b·d ≤ usable·16 ⇒ predicted ≤ budgeted ≤ budget` for *any* σ, k, static. So v1's
sweep was vacuous, and `resolve_layout_predicted_peak_fits_budget` (`clump.rs:577`) is a tautology that
the margin has made *more* trivially true. Splitting into 4a/4b/4c is the correct response, and
Assumption 2's observation — that exactness is load-bearing in two places with opposite signs — is the
sharpest sentence in the document.

**One correction to my own parent-plan review, which v3 may have inherited:** I recommended adding a
`cores = 1` calibration cell and called it "the single most informative cell in the design". **That cell
is unreachable on the clumpify path** (`cli.rs:743-745`). The identifiability lever is
`cores ∈ {2,3,4}` — three points, not four, with `cores = 3` the only unused one. This plan's
Assumption 3 already has it right.

---

## 7. Smaller notes

**PR ordering is safe, and the plan should say why.** D4 ships docs first while D2 changes stderr in
PR 2, which invites "do the docs then describe a banner that does not exist yet?". Checked: the only
verbatim stderr block PR 1 touches is `clumpy.md:113`, the *fallback warning* (`need ≥ 552` → `311`),
which D2 does not change. `clumpy.md:102` says only "the binary prints those resolved values at
startup". So there is no window where PR 1's docs contradict the shipped binary. One clause recording
that closes the question. **[Optional]**

**§8's multi-pair cell has no failure disposition.** The reasoning is good — FastQC is an addend
*because* mimalloc keeps the pool's pages resident across a phase boundary, so pair 2's pool may sit on
pair 1's retained pages — and nothing tests it. But the cell is "logged, non-gating" with no statement
of what happens if it breaches. Given the mechanism is the same one that forced FastQC to be an addend,
a breach here is likely and would mean `--memory` is a per-pair rather than per-process bound, which is
a documentation change at minimum. **Say what a red means before running it.** **[Important]**

**Validation 9's fallback is the right instinct.** Reading the floor's top end off
`--cores 32 --fastqc --memory 1011M`'s warning rather than trusting a reimplementation is exactly the
discipline this arc has been missing; my own table above is a reimplementation and therefore weaker
evidence than that invocation would be. **[Optional]**

**D1 is right and its narrowness is well argued.** A fallback that byte-copies the input at full I/O
cost and exits 0 would be worse than a hard failure, and the corrected regression table shows the
surface is three core counts wide with a 154 MiB worst-case window. The `--cores 1` edge case (floor
277 MiB, "or decrease --cores" unactionable) is caught. No objection. **[Optional]**

---

## Action items

### Critical
1. **Sever §0 throughout** — strike the twelve dependents in §1's first table; keep the five survivors
   and say explicitly that they survive. Re-scope D3's issue from data-integrity to user-visible wrong
   behaviour, rewrite the `.fq.gz` assertion's rationale as measurement validity, and rest the gzipped
   fixture on the disk argument alone.
2. **Add a free-disk precondition to the measurement harness spec** (≥ 3× expected write volume,
   recorded per row). All four retractions in this arc were instrumentation defects; this is the one
   that would have caught §0.
3. **Stage PR 3**: land with `continue-on-error: true`, collect ≥10 runs across cold and warm cache,
   flip to gating in a follow-up recording the observed Linux range. Validation 14's "stop" is not an
   operational branch on a required job.

### Important
4. **Give the 25 bp cell a band assertion** so a worsening fails, and file the short-read
   under-reservation (σ ≈ 1.55 vs shipped 1.4375) as its own issue.
5. **Record `2G/c4` as a deliberate blind spot** in Remaining risks and give its residual a tracking
   issue. Keep the decision to drop it.
6. **Extend Validation 12** to run the depth mutation as a re-inlined literal at `parallel.rs:117` as
   well as a helper edit, and add to §5 that the CI job is the only guard catching the re-inlined form —
   it is PR 3's strongest justification and is currently absent.
7. **Pin the sum** `work + 1 + result_per_core == 4` in 4c, message naming `k`, alongside the two value
   pins.
8. **State what a multi-pair red means** before running that cell.

### Optional
9. Record that PR 1 → PR 2 ordering is safe because `clumpy.md:113` is the fallback warning, not the
   banner.
10. Add the `:1193`/`:1202` proximity and the #434 linkage to the EOF-send issue.
11. Note that narrowing `ClumpLayout::new` is semver-breaking in principle; harmless here.
12. Correct any inherited "add a `cores = 1` cell" recommendation from my parent-plan review — that
    cell is unreachable under `--clumpify`.

---

*Reviewer A. Code claims verified against `fdade7e` with commands quoted; floor/bin/banner arithmetic
reimplemented in Python from the constants as read. No `trim_galore` run was executed by me, so no RSS
figure in this document is mine, and every measurement in the plan is marked **[NV]** above.*
