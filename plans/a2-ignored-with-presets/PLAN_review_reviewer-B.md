# PLAN review — Reviewer B

**Plan:** `plans/a2-ignored-with-presets/PLAN.md`
**Target:** FelixKrueger/TrimGalore#369 — `-a2` ignored when `--illumina` is used
**Repo state:** `dev` @ `4276810`, `cargo build --release` run fresh before every smoke test
**Method:** every load-bearing claim below was checked by reading source, grepping, or running the release binary. Scratch output went to the session scratchpad, outside the repo. No repo file was modified except this report.

---

## Verdict

The diagnosis is correct and unusually well evidenced — I reproduced §2.2 exactly, all seven rows. The fix shape (one override point reached by all branches) is the right design and I would not propose an alternative.

What the plan is missing is **branch coverage of the resolution space it is changing**. It enumerates "five presets + explicit `-a` + auto-detect" and treats that as exhaustive. It is not: `--consider_already_trimmed` reaches the auto-detect branch through a *synthetic* preset, and the plan's unconditional-override rule produces an asymmetric, self-contradictory result there. Separately, the plan introduces a class of change it explicitly promises not to introduce (turning working invocations into hard failures), and its headline regression guard (V7) is structurally blind to the one tuple field most likely to be mis-wired by the reshape.

Two Critical, nine Important, four Optional.

---

## 1. Claims I verified as sound

Stated briefly so the plan gets credit and the caller can skip re-checking.

| Claim | Verdict | How |
|---|---|---|
| §2.2 seven-row discard table | **Confirmed exactly** | Ran all seven invocations with `-a2 AAATCAAAAAAAC` on `test_files/BS-seq_10K_R{1,2}.fastq.gz`, read `Command line parameters` from each R2 report. Every sequence matches the plan's table, including the auto-detect row. |
| §2.1 root cause | Confirmed | `src/main.rs:985-1048`; `cli.adapter2` is read at **exactly one** place in the whole crate (`src/main.rs:988`) plus the `--clump_only` rejection at `src/cli.rs:716`. Fallback at `src/trimmer.rs:112` reads as quoted. |
| §2.3 not a v2.3.0 regression | Confirmed | `git show v2.1.0:src/main.rs`, `v2.2.0`, `v2.3.0` — the `if cli.illumina` branch is byte-identical at all three and on `dev`. |
| §2.4 `unless (defined $a2)` is the contract | Confirmed | `git show 0.6.11:trim_galore` lines 562, 572 (smallRNA / BGISEQ R2 defaults) and 584 (`$a2 = ''` fallback). |
| §2.4 `-a` + preset → `die` | Confirmed | 0.6.11 lines 3131-3147, five separate `die`s (illumina / nextera / small_rna / stranded / bgiseq). |
| §2.4 `-a2` without `--paired` → `die` | Confirmed | 0.6.11 lines 3151-3153. And confirmed v2.x silently accepts: `--illumina -a2 AAATCAAAAAAAC single.fq.gz` → exit 0, no warning, no `Adapter 2` line. |
| §2.6 display code already prints the R2 line | Confirmed | `src/main.rs:833-843`; the `eprintln!("Adapter 2 (Read 2): {}", …)` is at **`main.rs:835`**, gated only on `!adapters_r2.is_empty()`. The plan's cited range contains it. |
| **A4** — `--small_rna` cutoff is R1-driven | **Confirmed** | `src/main.rs:848-850`: `let first_adapter_seq = adapters_r1.first()…; if first_adapter_seq == "TGGAATTCTCGG" { … 18 }`. Keyed on `adapters_r1` only; `-a2` cannot reach it. Matches Perl, which also keys on the resolved R1 *sequence* (0.6.11:552). A4 is safe. |
| A5 conclusion — no CI case passes `-a2` | Confirmed (but see I6) | `grep -nE '(^\|[[:space:]])(-a2\|--adapter2)([[:space:]=]\|$)' .github/workflows/ci.yml` → no matches. |
| §7 downstream consumers need no change | Confirmed | `trimmer.rs:112` and `TrimConfig::r2_adapter_count` (`trimmer.rs:49-55`) both key on `is_empty()`; `r2_adapter_count` is called only from paired paths (`trimmer.rs:400,711`, `parallel.rs:221,317`). |
| §4 step 5 — `main.rs` has no `#[cfg(test)]` | Confirmed | `grep -n "cfg(test)" src/main.rs` → no matches. |
| §4 step 7 — `#### Fixes` exists under `### Unreleased` | Confirmed | `CHANGELOG.md:84`. |
| §4 comment-style paragraph | Consistent with `CLAUDE.md` | The proposed `// user -a2 always wins over preset/auto-detected R2 defaults` is one line stating the fact. The NOTE/WARNING strings are user-facing output, not comments — no conflict, and their register is plain. |

Two further things in the plan's favour that it does not claim:

- **No existing test asserts the buggy behaviour.** The only `adapter2` reference in any test is `src/cli.rs:1345`, a pure parse assertion (`cli.adapter2 == vec!["GCAT","AAAA"]`). The fix requires no test rewrites.
- **The plan's "silent empty override" hazard is provably unreachable.** I went looking for a path where `cli.adapter2` is non-empty but `parse_adapter_specs` returns `Ok(vec![])` — which would silently wipe a preset's R2 default and fall back to R1's adapter. There is none: `adapter.rs:328-330` bails on an empty spec, `:320-322` on an all-empty multi-adapter string, and `read_fasta_adapters` bails on both an empty sequence (`:443-448`) and a record-free file (`:472-477`). So keying the override on `cli.adapter2.is_empty()` (§3.1 step 2) rather than on the parsed result is safe. Worth stating in the plan so an implementer does not "defensively" add a second check.

---

## 2. Logic review — findings

### CRITICAL

#### C1. `--consider_already_trimmed` is an unenumerated resolution branch, and the unconditional override breaks it

**Confirmed wrong.** The plan's §3.1/§3.2 treat auto-detect as one branch with one behaviour. It is two.

`adapter.rs:219-227` — when `--consider_already_trimmed <N>` suppresses trimming, auto-detection returns a **synthetic preset**:

```rust
let final_preset = if suppressed {
    AdapterPreset { name: "already trimmed (adapter trimming suppressed)", seq: "", seq_r2: None }
} else { preset };
```

So `adapters_r1 = [(label, "")]` — one entry with an empty sequence, which `trim_read` skips (`trimmer.rs:~127`, `if adapter_seq.is_empty() { continue; }`) — and `adapters_r2 = []`.

Current behaviour, verified by running `--paired --consider_already_trimmed 10000 -a2 AAATCAAAAAAAC`:

```
No auto-detected adapter sequence exceeded the user-specified 'already adapter-trimmed' limit of 10000 counts. Only quality trimming will be carried out.
Adapter: already trimmed (adapter trimming suppressed) ()
R1 report: Command line parameters: … -O 1 -a  BS-seq_10K_R1.fastq.gz
R2 report: Command line parameters: … -O 1 -a  BS-seq_10K_R2.fastq.gz
```

Apply §3.1 step 2 verbatim ("If `cli.adapter2` is non-empty, it replaces the candidate, **unconditionally, on every path**") and you get:

- R1: adapter trimming suppressed (empty sequence).
- R2: **full adapter trimming with `-a2`**.
- Display: `Adapter: already trimmed (adapter trimming suppressed) ()` immediately followed by `Adapter 2 (Read 2): AAATCAAAAAAAC`.
- The §3.3 NOTE does **not** fire, because `seq_r2` was `None` — nothing was "displaced". So nothing announces the conflict.

I could not execute the post-fix behaviour (the fix does not exist), so the asymmetry is *derived* from the plan's stated rule plus the code I read — but the derivation is mechanical and I am confident in it.

There is a defensible argument for the override: `-a SEQ --consider_already_trimmed N` **already** bypasses suppression today, because the `-a` branch returns at `main.rs:986-991` before auto-detection ever runs. So "an explicit adapter beats suppression" is already established for R1, and extending it to R2 is consistent. But the *asymmetry* (R1 untrimmed, R2 trimmed, label contradicting the R2 line) is the part no user would predict, and it arises only from passing `-a2` alone.

**Action:** decide explicitly and put it in §3.2 as a fifth row. Either (a) allow the override and fix the display so the label/NOTE reflect it, or (b) skip the override when `detection.suppressed` and warn. `DetectionResult.suppressed` (`adapter.rs:78`) is already available at the decision point, so (b) costs one condition. Add a V-row either way — this is the one case where the plan's own semantics produce output no test currently covers.

#### C2. The fix converts currently-succeeding invocations into hard errors — the exact thing §3.4 promises not to do

**Confirmed by running the binary.** Today, `-a2` is only validated on the `-a` branch, so garbage in `-a2` is silently discarded everywhere else:

```
$ trim_galore --illumina -a2 'ZZZQQQ' illumina_10K.fastq.gz     → exit 0, output written
$ trim_galore -a AGATCGGAAGAGC -a2 'ZZZQQQ' illumina_10K.fastq.gz
  Error: Adapter sequence must contain only DNA characters (A, C, G, T, N, X), got: 'ZZZQQQ'
```

§3.1 requires `-a2` to go through `parse_adapter_specs` **on every path** (correctly — that is A3). The consequence is that `--illumina -a2 ZZZQQQ`, `--nextera -a2 ""`, and `-a2 file:missing.fa` all move from exit 0 to a hard error. Real command lines with a typo'd or stale `-a2` that currently run to completion will start failing.

This directly contradicts the plan's own stated principle. §3.4: *"this change should not convert working invocations into failures"* — used there to justify a warning over Perl's `die`. §7 lists four user-visible changes; all four are output/stderr changes. The exit-code change is a fifth class and is absent.

I think erroring is the **right** call — silently ignoring a malformed adapter is worse. But it must be an owned decision with a CHANGELOG line, not a side effect discovered by a user. It also slightly weakens the §3.4 argument, which should be re-stated on its own merits rather than on "we never break working commands".

**Action:** add to §7 as change 5; add a CHANGELOG sentence; add a V-row asserting `--illumina -a2 ZZZQQQ` exits non-zero with the DNA-characters message.

### IMPORTANT

#### I3. Specialty modes still silently discard `-a2`, and the §3.4 warning misses them

**Confirmed by running.** §1's secondary goal is "stop silently discarding `-a2` when it cannot apply". §3.4 implements that only for `!cli.paired`. But:

```
$ trim_galore --hardtrim5 30 --paired -a2 AAATCAAAAAAAC R1.fq.gz R2.fq.gz
  → exit 0, hard-trims both files, no mention of -a2 anywhere
```

`--hardtrim5/3`, `--clock`, `--implicon` bypass the trimming pipeline entirely and never call `resolve_adapter`. With `--paired` set, the §3.4 warning does not fire, so this remains a silent-discard path after the fix — the same defect class the plan says it is closing.

The mirror case is worse. **Without** `--paired`, `--hardtrim5 30 -a2 SEQ R1 R2` *will* fire the new warning, with a wrong explanation: "applies to Read 2 of a pair and is ignored without `--paired`". It is ignored because hardtrim does no adapter trimming at all. (Verified: `--clock -a2 SEQ R1 R2` and `--hardtrim5 30 -a2 SEQ R1 R2` both run to completion today.)

**Action:** either extend the warning to "given but not used" for specialty modes with a mode-appropriate reason, or make the message reason-neutral so it is true in both cases. Note `--clump_only` already gets this right with a hard error (`cli.rs:716`) — that is the precedent to follow.

#### I4. V7 is structurally blind to the tuple field most likely to be mis-wired

**Confirmed by running.** §11 names the reshape touching Read 1 as the Medium risk and V7 (md5 R1 before/after) as the guard. V7 covers all seven branches' R1 *sequences*, which is good. But the reshape moves **four** fields per branch, and field 3 (`Option<(usize, usize)>`, the poly-G piggyback) is invisible to R1 md5 on the plan's own fixture.

Field 3 is `Some(…)` **only** on the auto-detect branch — an optimisation to avoid a second file scan (`main.rs:900-905`). Verified empirically:

```
--illumina  run: "Scanning for poly-G content..." appears 1×
auto-detect run: "Scanning for poly-G content..." appears 0×   (piggybacked)
```

Two mis-wirings, neither caught by V7:

1. **Field 3 dropped on the auto-detect branch** → falls through to `adapter::detect_poly_g(input_file)` at `main.rs:904`, a second full pass over up to 1M reads per input file. Trimmed output is byte-identical, so V7 passes, `cargo test` passes, and §6's claim that the restructure "adds no I/O" is silently false.
2. **Field 3 set to `Some((0, 0))` on a preset branch** → `threshold = (0/10_000).max(10) = 10`, `0 > 10` is false → poly-G trimming silently **disabled** where auto-detection would have enabled it. On 2-colour data that changes trimmed output.

Case 2 is invisible on the plan's fixture: `BS-seq_10K_R1` reports `Poly-G trimming: not enabled (auto-detection found 2 of 10000 reads (0.02%) …)` — below the threshold of 10, so poly-G is off either way and R1 md5 is identical whether field 3 is right or wrong.

**Action:** add a V-row asserting the `Scanning for poly-G content` line count per branch (0 for auto-detect, 1 for each preset) and that the `Poly-G trimming:` verdict line is unchanged before/after. This is cheap and it is the only check that sees field 3.

#### I5. V13 under-scopes the validation job by roughly 6×

**Confirmed.** V13 says "Re-run the five existing CI matrix invocations" and A5 describes the job as SE / PE / hardtrim5 / clock / demux. That summary comes from `CLAUDE.md`, which is out of date. The actual job has ~30 steps named `Validate …`, and `ci.yml` contains 24 byte-identity assertions (14 `md5sum` lines), including `--clump_only` byte-identity SE and PE, `--clump_only` cross-run determinism, `--clump_only` uBAM record parity, uBAM `@PG` chain preservation, and the collision pre-flight cases.

The fix should not move any of them, but "five invocations" is the wrong instruction to hand an implementer.

**Action:** restate V13 as "run the whole `validation` job" (locally or on a branch push). Simpler to write and strictly stronger.

#### I6. §9 does not distinguish one-off manual checks from permanent regression gates

The root cause of #369 is a *missing permanent gate* — §2.5 says so explicitly. Yet §9 reads as a manual checklist, and §4 step 5's test prescription is "unit coverage for the resolution matrix in `main.rs`'s or `adapter.rs`'s test module — whichever is reachable".

A unit test of a pure override helper **cannot catch the #369 bug class**. #369 was not "the helper computes the wrong answer"; it was "six branches never called the helper". Only a per-branch end-to-end assertion catches that. V1-V3 are exactly the right assertions and they need to live in `tests/`, not in a plan table.

The pattern already exists: `tests/integration_clump_only.rs` and `tests/integration_passthrough.rs` both run the binary and grep trimming reports.

**Action:** promote V1, V2, V3, V9 and the C1/C2 rows to a new `tests/integration_adapter2.rs`, asserting the R2 report's `Command line parameters` line for all seven branches. Keep the unit test of the override helper as well, but do not let it substitute.

#### I7. Docs (step 8) is too vague, and the reality is worse than the plan suspects

Step 8 asks whether any page "documents `-a2` in a way that implies it works with presets". Verified answer: the docs document `-a2`-without-`-a` as the **recommended** form, and that form is a complete no-op today.

- **`docs/src/content/docs/guide/adapters.md:31`** — presented as "The cleanest (v2.x)": `trim_galore --paired -a2 AGCTAGCG -a2 TCTCTTATAT -a2 TTTCGGATTTAT -n 3 R1.fq.gz R2.fq.gz`. No `-a` → auto-detect branch → all three `-a2` values discarded. The documented example does nothing.
- **`adapters.md:38`** — same problem for the embedded-string form: `-a2 " AGCTAGCG -a … " --paired R1 R2`, no `-a`.
- **`adapters.md:50`, `flags.md:27`, `flags.md:39`** — claim `A{N}` expansion works for `-a2`. Verified false unless `-a` is also given: `--illumina -a2 'A{10}'` prints no expansion line at all.
- **`flags.md:12`** — "`--small_rna` … auto-sets `--adapter2` to the Illumina small RNA 5' adapter (`GATCGTCGGACT`) for paired-end data." After the fix this needs the "unless you supply your own" hedge. The adjacent `flags.md:13` already uses exactly that phrasing for `--clip_R2` — copy it.
- **`reference/migration.md:26`** — sells "Repeatable `-a` / `-a2`" as a v2 improvement over Perl. True only with `-a`.

That two documented examples were no-ops for three releases is worth a CHANGELOG sentence in its own right.

#### I8. §2.5's reproduction command does not do what the plan says

**Confirmed wrong, minor.** The plan states: *"`grep -n 'a2\|adapter2' .github/workflows/ci.yml` returns **nothing**."* It returns **8 matches** — spurious hits inside `sha256sum` and the pinned `dtolnay/rust-toolchain@29eef336d9b2848a0b548edc03f92a220660cdb8` SHA.

The conclusion is right; the command is not. A reader re-running it gets output and may conclude the opposite. Replace with the anchored form: `grep -nE '(^|[[:space:]])(-a2|--adapter2)([[:space:]=]|$)' .github/workflows/ci.yml` (no matches).

#### I9. SE `-a2` produces self-contradictory output after the fix, and V10 cannot see it

**Confirmed by reading.** `main.rs:833-843` gates the R2 display **only** on `!adapters_r2.is_empty()` — there is no `cli.paired` guard. After the fix, `resolve_adapter` populates `adapters_r2` from `-a2` on the SE path too (`setup_trimming` is called for SE at `main.rs:802` and `:1836`; it discards the R2 list but prints from it first). So a single-end run with `-a2` will emit both:

```
WARNING: -a2/--adapter2 applies to Read 2 of a pair and is ignored without --paired.
...
Adapter 2 (Read 2): AAATCAAAAAAAC
```

V10 asserts "Warning on stderr, exit 0, output byte-identical to omitting `-a2`". The FASTQ *is* byte-identical (`trim_read` gates on `is_r2`, and the SE report hardcodes `adapters_r2: Vec::new()` at `main.rs:1161`), so V10 passes while the stderr contradicts itself.

**Action:** guard the display on `cli.paired`, or do not populate `adapters_r2` when `!cli.paired`. Extend V10 to assert the absence of the `Adapter 2 (Read 2):` line.

#### I10. §3.5's "empty `-a2`" bullet is ambiguous and, on one reading, contradicts C2

§3.5 says *"Empty `-a2` on a preset with an R2 default: still gets the default."* Two readings:

- flag absent (`cli.adapter2` is an empty `Vec`) → default stands. Correct, and V5 tests it.
- `-a2 ""` (vec is `[""]`, non-empty) → `parse_adapter_specs` bails with "Empty adapter sequence" (`adapter.rs:328-330`). **Not** "still gets the default".

§11 lists "`-a2` given but empty after parsing" as an enumerated edge case but no V-row covers it and no expected behaviour is stated. Disambiguate the bullet; the second reading is the C2 decision.

#### I11. §7 does not mention the uBAM output path, which also inherits the change

**Confirmed by reading.** `setup_trimming` — and therefore `resolve_adapter` — is called from six sites: `main.rs:674, 771, 802` (FASTQ dispatch) and `main.rs:1754, 1806, 1836` (the `--output-format ubam` dispatch). So `--paired --illumina -a2 SEQ --output-format ubam` changes behaviour too, and §7's integration section does not say so.

**Action:** one sentence in §7; one V-row exercising `--output-format ubam --paired --illumina -a2 SEQ` (paired uBAM output is one interleaved BAM, so the check is on the report, not on a second file).

---

## 3. Assumptions

| ID | Status |
|---|---|
| A1 (R1 must not move) | Sound as a requirement; the guard (V7) is insufficient — see I4. |
| A2 (`unless (defined $a2)` is the contract) | **Verified.** Two independent sources as claimed. See O4 for a nuance the plan omits. |
| A3 (`-a2` parse should match `-a`) | Correct, and it is the right call — but it is the direct cause of C2, which the plan does not connect. |
| A4 (`--small_rna` cutoff is R1-driven) | **Verified sound.** `main.rs:848-850`. The plan was right to flag it and right about the answer. |
| A5 (no CI case passes `-a2`) | Conclusion verified; stated command wrong (I8); V13's scoping wrong (I5). |
| A6 (`-a` + preset stays permissive) | Sound, user-confirmed, verified as the current behaviour (`--illumina -a AAAACCCCGGGG` → exit 0, `-a` wins). |

**Unstated assumptions I found:**

- **U1.** That "auto-detect" is a single branch. False — `--consider_already_trimmed` makes it two (C1).
- **U2.** That `resolve_adapter`'s four returned fields are independent, so reshaping for field 2 cannot disturb fields 0/1/3. Field 3 is coupled to a downstream I/O decision (I4).
- **U3.** That `-a2` is only ever discarded on the trimming path. Specialty modes discard it too (I3).
- **U4.** That making `-a2` parse everywhere is behaviour-preserving for valid input. It is — but not for invalid input (C2).

---

## 4. Efficiency

§6 is correct as far as it goes: one extra `parse_adapter_specs` over a near-always-empty vector, once per input file, outside any per-read path. No concern.

The one claim that is *conditionally* false is "adds no I/O" — true only if field 3 stays correctly wired. If it does not, the auto-detect path silently gains a full second pass over up to 1M reads per input file (I4, case 1). Nothing in §9 checks this. That is an efficiency claim resting on an unverified invariant, which is precisely how the original defect shipped.

The 5-tuple has no runtime cost.

---

## 5. Validation sufficiency

Where the code could silently produce wrong results and no proposed validation would notice:

| Failure mode | Caught by plan as written? | Fix |
|---|---|---|
| Field 3 dropped on auto-detect → duplicate full scan | **No** — output byte-identical | Assert `Scanning for poly-G content` count per branch (I4) |
| Field 3 wrongly `Some` on a preset → poly-G silently off | **No** on the chosen fixture (2/10000, below threshold) | Same, plus assert the `Poly-G trimming:` verdict line |
| Adapter *label* swapped while sequence stays right | **No** — V7 md5s the FASTQ, not the report | md5 the `*_trimming_report.txt` too, not just the `.fq.gz` |
| `--consider_already_trimmed` + `-a2` asymmetry | **No** — no V-row, no truth-table row | C1 |
| `-a2` silently dropped in specialty modes | **No** | I3 |
| Invalid `-a2` newly fatal | **No** | C2 |
| SE stderr self-contradiction | **No** — V10 only checks the FASTQ | I9 |
| uBAM output path | **No** | I11 |

Two structural points:

- **V7 should md5 the reports, not only the FASTQ.** The label (tuple field 0) flows into `Adapter: {label}` and into the report config's `adapters`. A reshape that returns the right sequence with the wrong label is invisible to a FASTQ md5 and obvious in a report diff. One extra file per invocation, no extra runs.
- **V15 is the best row in the table** and the plan is right to single it out. Keep it exactly as written.

---

## 6. Alternatives

**On §5's 5-tuple.** The plan's reasoning ("already at the edge of readable… not worth the churn now") is undercut by a fact it does not state: `resolve_adapter` has **exactly one caller** (`main.rs:818`, verified by grep). Converting to a named struct is a ~6-line change, not churn. Two concrete improvements, in increasing order of ambition:

1. **Shrink the new field to `Option<String>`.** The plan proposes `Option<(String, String)>` = (preset name, displaced seq). The preset name is *already* returned as field 0 in every case where a displacement can occur, so the tuple duplicates it. Just the displaced sequence suffices. Smaller diff, no lost information.
2. **Introduce `struct ResolvedAdapters { label, r1, r2, poly_g_counts, displaced_r2 }`.** With one call site this is cheap, and it directly mitigates I4: named fields make a mis-wired `poly_g_counts` a visible error at the construction site rather than a positional slip in a 5-tuple. Given that I4 is the plan's own Medium risk, the struct is the cheapest mitigation available — it converts a class of runtime bug into something a reader can see. `SetupResult` (also a 4-tuple, `main.rs:77`) does not need to change either way, since §4 step 3 has `setup_trimming` consume the displacement locally.

I'd take (2). If the maintainer prefers minimal diff, (1) is strictly better than the plan as written.

**On where the NOTE is printed.** The plan puts it in `setup_trimming` after the display, and returns the displacement to get it there. An alternative is printing it inside `resolve_adapter`, which already prints (`main.rs:1035-1037`: `"Auto-detecting adapter type..."` and `detection.message`) — that would avoid the signature change entirely. But it would place the NOTE *before* `Adapter 2 (Read 2):`, which reads worse. **The plan's choice is correct**; I mention this only to record that the obvious shortcut was considered and is worse.

**On C1.** The narrower alternative to an unconditional override: skip it when `detection.suppressed`. `DetectionResult.suppressed` (`adapter.rs:78`) is already in scope, so this is one condition. It preserves `--consider_already_trimmed`'s meaning at the cost of one documented exception to "`-a2` always wins". Either answer is defensible; the plan must pick one.

---

## 7. Action items

### Critical

1. **C1** — Add `--consider_already_trimmed` as an explicit row in §3.2 and decide whether `-a2` overrides suppression. If yes, fix the label/NOTE so the output is not self-contradictory. Add a V-row.
2. **C2** — Own the exit-code change: `-a2` with invalid content goes from silently ignored to fatal on preset and auto-detect paths. Add to §7 as change 5, add a CHANGELOG line, add a V-row. Re-argue §3.4's warning-over-error choice on its own merits, since "we never break working commands" is no longer true.

### Important

3. **I3** — Close or correctly warn about `-a2` in specialty modes (`--hardtrim5/3`, `--clock`, `--implicon`), including with `--paired`. Fix the §3.4 message so its stated reason is true in every case it fires.
4. **I4** — Add a V-row for the poly-G piggyback: `Scanning for poly-G content` count per branch (0 auto-detect / 1 per preset) and an unchanged `Poly-G trimming:` verdict. Strengthen V7 to md5 the `*_trimming_report.txt` as well as the `.fq.gz`.
5. **I5** — Restate V13 as "run the whole `validation` job" (~30 steps, 24 byte-identity assertions), not "the five existing invocations".
6. **I6** — Promote V1/V2/V3/V9 and the C1/C2 cases into `tests/integration_adapter2.rs`, following `tests/integration_clump_only.rs`. State in §4 step 5 that a unit test of the override helper cannot catch the #369 bug class.
7. **I7** — Replace step 8's "check the docs" with the verified edit list: `guide/adapters.md:31,38,50`; `guide/flags.md:12,27,39`; `reference/migration.md:26`. Note in the CHANGELOG that two documented examples were no-ops.
8. **I8** — Correct §2.5's grep to the anchored form.
9. **I9** — Guard the `Adapter 2 (Read 2):` display on `cli.paired` (or don't populate `adapters_r2` in SE). Extend V10 to assert that line's absence.
10. **I10** — Disambiguate §3.5's "empty `-a2`" bullet; add the `-a2 ""` V-row.
11. **I11** — Add the uBAM output path to §7 and a V-row for `--output-format ubam --paired --illumina -a2 SEQ`.

### Optional

12. **O1** — Shrink §5's new field to `Option<String>`; preferably convert to a named struct (one call site; also mitigates I4).
13. **O2** — Align §3.4's wording with the house pattern at `cli.rs:939-963` (`WARNING: <flag> … Ignoring.`). Good news: placing a warning in `validate()` is well precedented — that function already ends with a deprecation-warning block, so this is not a new side effect and the plan's siting is right.
14. **O3** — Add a `--poly_a` + `-a2` V-row. Verified: `--paired --poly_a -a2 SEQ` currently gives R2 the auto-detected `AGATCGGAAGAGC`; after the fix it gets `SEQ`. Perl **dies** on this invocation (0.6.11:3088-3095 requires both `-a` and `-a2` or neither under `--polyA`). Rust's `--poly_a` is a standalone homopolymer step (`trimmer.rs:196-213`, `quality.rs:172-174`) independent of adapters, so a user porting Perl's `-a2 "T{150}"` idiom gets a 150-T adapter alignment **plus** the built-in poly-T 5' clip. Pre-existing on the `-a` branch; the fix extends it to the `-a2`-only branch.
15. **O4** — Record one more Perl divergence next to A6: 0.6.11 keys its R2 defaults on the resolved R1 *sequence* (`if ($adapter eq 'TGGAATTCTCGG')` at :552, `if ($adapter eq 'AAGTCGGAGG…')` at :569), not on the preset flag — so Perl sets the smallRNA R2 default even for an explicit `-a TGGAATTCTCGG`. Rust keys its *length cutoff* the same way (`main.rs:848-850`, which is why A4 holds) but not the R2 default. Out of scope; worth recording so it is not rediscovered as a bug.
