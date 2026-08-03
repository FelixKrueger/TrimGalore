# PLAN review — Reviewer A

**Plan:** `plans/a2-ignored-with-presets/PLAN.md` (#369, `-a2` ignored with presets)
**Repo state:** `dev` @ `4276810`, `cargo build --release` run fresh before every measurement below.
**Method:** every claim marked *verified* was checked by reading the cited source, by `git show`, or by running `./target/release/trim_galore` on `test_files/BS-seq_10K_R{1,2}.fastq.gz` with scratch output under `$TMPDIR`. No repo file was modified. Claims marked *inferred* are reasoning from code I read but did not execute.

## Verdict

The diagnosis is correct and unusually well evidenced — I reproduced all seven rows of §2.2 exactly, and every load-bearing claim in §2.3–§2.6 and A4/A5 holds. The fix direction is right.

But the plan's central safety property is **factually wrong**, and it is the one the plan itself nominates as "the one property whose breakage would be worst and least obvious". §3.5, A1 and V7 all assert that Read 1 output is byte-identical after the fix. On paired-end input it is not, and cannot be: the joint pair-level length filter drops different pairs once Read 2 is trimmed differently. I measured a 135-pair swing on the project's own fixture. An implementer running V7 as written will see a Read 1 md5 mismatch on V1, V2 and V3 and will have no way to tell "expected consequence of the fix" from "I broke Read 1".

There is also a cheaper and lower-risk shape for the fix than §4/§5 propose, which removes the plan's own Medium risk and the 5-tuple entirely.

---

## 1. Logic review

### C1 (Critical) — "Read 1 output is byte-identical" is false for paired-end. §3.5, A1, V7 must be re-specified.

*Verified by running the binary.* Two invocations that trim Read 1 with the identical adapter (`AGATCGGAAGAGC`):

| Invocation (current `dev`) | R1 `_val_1.fq.gz` md5 | R2 `_val_2.fq.gz` md5 | pairs dropped by length filter |
|---|---|---|---|
| `-a AGATCGGAAGAGC -a2 AAATCAAAAAAAC` | `a571f332bad6343d62f1ec37120e1a27` | `1504b7c1f2c34d9d605b4db702537213` | **139** (1.39%) |
| `--illumina -a2 AAATCAAAAAAAC` | `f8d879e31d9b797667796cec0f208dd5` | `5b94fb44b2ea4d415ba3e99f2b49bfea` | **4** (0.04%) |

Both R1 reports show the same cutadapt-stage figures — `Reads written (passing filters): 10,000 (100.0%)`, `Total written (filtered): 644,525 bp (99.2%)`. So R1's *per-read trimming* is genuinely identical; the `_val_1.fq.gz` files differ solely because the joint filter (`Number of sequence pairs removed because at least one read was shorter than the length cutoff`) discards 135 more pairs once R2 is trimmed with the user's adapter.

The second row is exactly what `--illumina -a2 SEQ` does today; the first row is what it must produce after the fix. **Post-fix, `_val_1.fq.gz` changes for every invocation in V1, V2 and V3 — by design.**

Consequences for the plan:

- §3.5 bullet 1 ("Read 1 adapter selection and precedence" unchanged) is fine as worded, but §2 line 20 ("Read 1 behaviour is correct today and **must be byte-identical after this change**") and A1 ("must not move… Testable via V7") are wrong as stated.
- V7's expectation column ("Identical in all cases") will fail on V1–V3. Re-specify as two separate properties:
  - **R1 trimming invariance** — for V1–V3, compare the R1 report's `Reads written (passing filters)` and `Total written (filtered)` bp, plus the `Reads with adapters` count. These must be unchanged.
  - **R1 file byte-identity** — only for invocations where `-a2` does not newly apply: V4 (`-a` + `-a2`, already correct), V5 and V6 (no `-a2`), and V13.
- §7 "Behaviour changes visible to users" should say the trimmed output changes for **both** reads of affected pairs, not just Read 2. A user seeing R1 change may otherwise reasonably file a follow-up bug.

This also has a real-world upside worth putting in the CHANGELOG: the reporter's ~20% mapping-rate drop is consistent with a *pair-retention* change of this kind, not just Read 2 mis-trimming.

### C2 (Critical, cheap) — a Perl-free differential oracle exists and the plan misses it; it also would have caught C1.

*Verified by running the binary* (table above). After the fix, `--illumina -a2 SEQ` and `-a AGATCGGAAGAGC -a2 SEQ` become the same invocation in every respect that reaches the FASTQ output: same R1 adapter sequence, same R2 adapter sequence, same length cutoff (keyed on the R1 *sequence* at `src/main.rs:848-856`), same poly-G decision (both return `None` for the auto-detect poly-G slot, so both fall through to `adapter::detect_poly_g`). Therefore:

> **Post-fix `--illumina -a2 AAATCAAAAAAAC` must produce R1 `a571f332bad6343d62f1ec37120e1a27` and R2 `1504b7c1f2c34d9d605b4db702537213`** — the md5s that today's `-a AGATCGGAAGAGC -a2 AAATCAAAAAAAC` already produces.

The same holds for the auto-detect row (`--paired -a2 SEQ`), because detection picks Illumina on this fixture (verified: `Command line parameters: … -a AGATCGGAAGAGC` in the auto-detect run's R1 report).

This is strictly stronger than V1/V2/V3 as written, requires no Perl install, and pins both reads at once. Compare the FASTQ outputs only — the reports legitimately differ in the adapter label (`user-specified` vs `Illumina`) and in the echoed command line.

It also settles V15 in advance: the chosen `-a2` sequence **does** discriminate on this fixture. R2 `Reads with adapters` moves 4,882 (48.8%) → 5,444 (54.4%) and the R2 md5 changes. A vacuous CI gate was a real risk — a rare `-a2` sequence would have produced identical md5s and a gate that guards nothing — and the plan does not currently record that it checked. Record the numbers.

### I3 (Important) — the override belongs in the caller; that deletes the plan's Medium risk *and* the §5 signature change.

*Inferred from reading `src/main.rs:815-849, 985-1049`.* §4 Step 1 reshapes all seven branches of `resolve_adapter` from early returns into a single value, and §11 books the resulting "restructure touches the Read 1 path incidentally" as the plan's only Medium risk. That risk is self-inflicted. `resolve_adapter` has exactly **one** caller (`setup_trimming`, `main.rs:817`), and the caller already holds the candidate R2 list immediately before the display block that has to print the NOTE:

```rust
let (adapter_label, adapters_r1, mut adapters_r2, autodetect_poly_g) =
    resolve_adapter(cli, input_file)?;

// user -a2 always wins over preset/auto-detected R2 defaults
let displaced = if !cli.adapter2.is_empty() {
    let user_r2 = adapter::parse_adapter_specs(&cli.adapter2)?;
    let d = adapters_r2.first().map(|(_, s)| (adapter_label.clone(), s.clone()));
    adapters_r2 = user_r2;
    d
} else {
    None
};
```

That is a ~10-line addition that touches **zero** existing branches, so Read 1 provably cannot move (nothing in the R1 path is edited), the 5-tuple in §5 is not needed at all (the caller computes `displaced` itself), and §3.1's "exactly one place where the Read 2 adapter is decided / none of them can forget" is satisfied just as well — the single place has simply moved one frame up. The only tidy-up is deleting the now-redundant `parse_adapter_specs(&cli.adapter2)` from the `-a` branch (`main.rs:988`), a one-line deletion.

Cost: `resolve_adapter`'s returned `adapters_r2` becomes a *candidate*, which needs one doc-comment line. That is a smaller readability debt than a 5-tuple.

I'd take this over §4/§5 unless there is a reason the plan hasn't stated.

### I4 (Important) — `-a2 ""` and malformed `-a2` become hard errors on the preset paths; §3.5 says the opposite.

*Verified by running the binary.* `--paired -a AGATCGGAAGAGC --adapter2 ""` → `Error: Empty adapter sequence` (from `adapter::parse_adapter_spec`, `src/adapter.rs:328-330`). `--paired --illumina --adapter2 ""` → **exit 0** today, because that branch never parses `-a2`. Post-fix, the second invocation starts failing.

Two problems:

- §3.5 bullet 4 — "Empty `-a2` on a preset with an R2 default: still gets the default" — is at best ambiguous ("no `-a2` given" vs `-a2 ''`) and at worst states behaviour the fix will not deliver. Disambiguate.
- §3.4's stated principle for choosing a warning over an error is "this change should not convert working invocations into failures". `--illumina -a2 ''` is precisely such a conversion. Either accept it deliberately (defensible — the `-a` path already errors, so this is consistency; but then say so and CHANGELOG it) or treat an empty/whitespace-only `-a2` as absent.

The same applies to invalid characters: `--illumina -a2 'ACGT!'` is silently swallowed today and will start erroring. Not enumerated anywhere in §9 or §11.

Related ordering point: if the override is applied "once before returning" (§4 Step 1), a malformed `-a2` is only diagnosed **after** the 1 M-read auto-detection scan has run. Parse `-a2` before the branch chain (or in `Cli::validate`) so the error is immediate. Cheap, and it matters on large inputs.

### I5 (Important) — `--consider_already_trimmed` + `-a2` yields a run whose own stderr message becomes false.

*Verified by running the binary.* `--paired --consider_already_trimmed 0 --adapter2 SEQ` today prints:

```
No auto-detected adapter sequence exceeded the user-specified 'already adapter-trimmed' limit of 0 counts. Only quality trimming will be carried out.
Adapter: already trimmed (adapter trimming suppressed) ()
```

and both reports show `Command line parameters: … -a  <file>` (empty adapter) with `Reads with adapters: 0 (0.0%)`. `src/adapter.rs:219-227` sets `seq: ""` and `seq_r2: None` on suppression, and `src/trimmer.rs:129-131` skips empty adapters.

Post-fix, `-a2` overrides unconditionally → **R2 is adapter-trimmed while R1 is not**, and "Only quality trimming will be carried out" is then untrue. *Verified against Perl:* 0.6.11 does the same thing (`trim_galore:2535` sets the R1 adapter to the literal `X`; `$a2` is untouched), so the behaviour is Perl-faithful and I would not change it. The **message** is the defect. Add the row to §3.2/§11 and either qualify the message when `-a2` is present or route it through the §3.3 NOTE.

### I6 (Important) — §3.4's warning condition is too narrow; specialty modes still discard `-a2` silently, including in paired mode.

*Verified by running the binary.* `--paired --clock --adapter2 SEQ` → exit 0, not one line of output mentions an adapter. `--hardtrim5 30 --adapter2 SEQ` → same. §1 names "stop silently discarding `-a2` when it cannot apply" as a secondary goal; a `!cli.paired` condition meets half of it. `--paired --clock -a2 SEQ` is paired, so it would emit no warning at all while still ignoring the flag.

Widen the condition to "not paired **or** a specialty mode is active" (`--hardtrim5/3`, `--clock`, `--implicon`; `--clump_only` is already a hard error at `src/cli.rs:716-718`, verified). Same site, same one-line change.

### I7 (Important) — the new SE warning will contradict output the binary already prints.

*Verified by running the binary.* Single-end `-a AGATCGGAAGAGC --adapter2 AAATCAAAAAAAC` today prints:

```
Adapter: user-specified (AGATCGGAAGAGC)
Adapter 2 (Read 2): AAATCAAAAAAAC
```

`setup_trimming`'s display block (`main.rs:833-842`) is not gated on `cli.paired`. Post-fix, `--illumina -a2 SEQ single.fq.gz` will print `Adapter 2 (Read 2): SEQ` *and* `WARNING: -a2/--adapter2 … is ignored without --paired`. Directly contradictory, on a tool whose stderr is read closely.

Gate the display on `cli.paired` (a one-line change that also fixes the pre-existing SE `-a`+`-a2` case), or emit the warning immediately after that line so the two read as one thought. §2.6 leans on this display block, so the plan needs to own its SE behaviour.

### I8 (Important) — better test siting than §4 Step 5: move `resolve_adapter` into the library.

*Verified:* `src/main.rs` contains no `#[cfg(test)]` (grep: no match), `Cli` is `trim_galore::cli::Cli` and public (`main.rs:9`), and `src/cli.rs:1329-1346` already builds a `Cli` in a unit test via `Cli::parse_from`. Step 5's fallback ("a pure helper belongs in `adapter.rs`") tests a fragment.

Moving `resolve_adapter` itself into the library (`adapter.rs`, or a small `resolve.rs`) makes the whole §3.2 matrix unit-testable via `Cli::parse_from(["trim_galore", "--paired", "--small_rna", "-a2", "SEQ", …])` — and six of the seven branches need **no** file I/O at all, since only the auto-detect branch touches the input. That is the difference between testing the override helper and testing the actual resolution matrix the bug lives in. It is also a pure move; the risk is low and clippy/`fmt` will catch any slip.

Aside, the existing `-a2` test at `cli.rs:1329` pairs `-a2` with `-a` — the one path that worked. Worth a sentence in the plan: the test suite encoded the bug's blind spot, exactly as CI did.

### Smaller logic points

- **§10.2 / §7: "v2.0.0" never existed.** *Verified* via `git tag`: the earliest v2 tag is `v2.1.0-beta.1`, and `CLAUDE.md` says master carries "v2.1.0 GA and onwards". Both §7 ("differ from v2.0.0–v2.3.0") and §10.2 ("affects v2.0.0–v2.3.0") name a release that does not exist, and this wording is destined for a user-facing CHANGELOG. Use "v2.1.0–v2.3.0" or "every v2 release".
- **§9 V13: "the five existing CI matrix invocations."** *Verified* — the `validation` job (`.github/workflows/ci.yml:244`) has more than five md5 steps plus several negative-path steps. Cosmetic, but V13 should name the steps it means so it is actually executable.
- **JSON report changes too.** *Verified:* `parameters.adapters_r2` exists in the JSON report (`src/report.rs:813`) and is populated today on the `--small_rna` path (`[{"name":"smallRNA_r2","sequence":"GATCGTCGGACT"}]`). It will change from `[]` to the user's sequence on the `--illumina`/`--nextera`/`--stranded_illumina`/auto-detect paths. Add to §7's list of four visible changes; it is a fifth.
- **uBAM output inherits the fix for free.** *Verified by running the binary:* `--paired --illumina -a2 SEQ --output-format ubam` completes and emits `*_val.bam`, and the paired uBAM path resolves adapters through the same `setup_trimming` (`main.rs:1806`). Worth one line so a reader does not wonder; a V-row is optional.
- **`--nextera` / `--nextseq` / `--rrbs` are not additional branches.** *Verified:* `cli.nextseq` reaches `TrimConfig` only as `nextseq: cli.nextseq.is_some()`; neither it nor `--rrbs` touches adapter resolution. The plan's enumeration of seven paths is complete.

---

## 2. Assumptions

| # | Verdict | Evidence |
|---|---|---|
| A1 | **Wrong as stated** | See C1. R1 *trimming* does not move; the R1 *file* does, via the joint pair filter. |
| A2 | Verified | `git show 0.6.11:trim_galore` — `unless (defined $a2)` at lines 562, 572, and the terminal default at 584. |
| A3 | Verified | `adapter::parse_adapter_specs` (`src/adapter.rs:356`) is reached for `-a2` only from `main.rs:988`. Brace expansion confirmed live: `-a2 'A{10}'` on the `-a` path prints `Adapter sequence A{10} expanded to AAAAAAAAAA` and the R2 report shows `-a AAAAAAAAAA`. |
| A4 | Verified twice | Code: `main.rs:848-856` keys the 18 bp cutoff on `adapters_r1.first()`, never on `adapters_r2`. Empirically: `--small_rna` with and without `-a2` both print the "Reducing length cutoff to 18bp" line and both report `length cutoff of 18 bp`. A4 holds. |
| A5 | Verified | `grep -n 'a2\|adapter2' .github/workflows/ci.yml` → no match. The only `illumina` hits are the fixture filename `illumina_10K.fastq.gz`. No matrix md5 can move. |
| A6 | Verified | `--paired --illumina -a AAAACCCCGGGG` → exit 0, `Adapter: user-specified (AAAACCCCGGGG)`. Permissive today, as the plan says. |

Unstated assumptions worth surfacing:

- **U1 — that `-a2` semantics are adapter-shaped in every mode.** *Verified in the Perl source:* `--polyA` sets `$a2 = extend_adapter_sequence("T",150)` at `trim_galore:518-522` with **no** `unless (defined $a2)` guard (so 0.6.11 *overwrites* a user `-a2` there — a genuine exception to §1's "always"), and it applies it as a **5'** adapter, `-g $a2`, at line 1197. In this codebase `--poly_a` is a boolean that drives a separate trimming stage (`quality::poly_a_trim_index`, `revcomp` for R2) and substitutes no adapter — *verified:* `resolve_adapter` has no `poly_a` branch. So there is **no missing branch**, and I confirm the plan's enumeration is complete on this point. But post-fix, a user carrying over the Perl-era `--polyA -a2 T{150}` recipe gets a 150 bp sequence fed into **3'** adapter matching (today it is silently dropped). One line in the plan as a documented divergence.
- **U2 — that Perl keys its R2 defaults off the preset flag.** It does not: `if ($adapter eq 'TGGAATTCTCGG')` (line 557) and `if ($adapter eq 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA')` (line 570). So 0.6.11 also supplies the smallRNA/BGI R2 default when the user writes `-a TGGAATTCTCGG` explicitly; this port does not (the `-a` branch supplies no R2 default). A separate pre-existing divergence in the same function, out of scope for #369 — but the plan claims to "restore the 0.6.11 contract", so it should say which part of the contract it is *not* restoring, or a future reader will assume the function is now Perl-faithful.
- **U3 — that the fixture discriminates.** Now verified (C2), but the plan asserted V15 without recording the measurement.

---

## 3. Efficiency

§6 is right about the resolution code: one `parse_adapter_specs` over a near-always-empty slice, once per input file, outside any per-read path. `parse_adapter_specs(&[])` returns an empty vec with no allocation of consequence. Nothing to object to.

Two things §6's "no measurable change" does not cover:

- **R2 adapter length is now user-controlled on the preset/auto-detect paths.** The DP is O(read_len × adapter_len) per read, and `src/alignment.rs:19` / `:272` limit the Myers bit-parallel prefilter to adapters ≤ 64 bp — longer patterns skip the prefilter and run the full scalar DP on every read. A `-a2 T{150}` (the Perl poly-A idiom, U1) therefore costs roughly an order of magnitude more per R2 read than the 13 bp Illumina adapter it replaces, with the prefilter disabled. Reachable today only via the `-a` path; the fix makes it reachable from every path. Not a reason to change the design — just don't claim "no measurable change" unqualified.
- **`-a2 file:…` is re-read once per input file / per pair,** since `resolve_adapter` runs per input. Already true for `-a`, so no regression, but a 20-pair run now opens the R2 FASTA 20 times.

---

## 4. Validation sufficiency

The V-table is unusually thorough for a bug fix of this size, and V15's instinct (a gate that passes pre-fix guards nothing) is exactly right. Gaps, in order of how quietly they would let a wrong implementation through:

1. **V7 is mis-specified and will fail on V1–V3** (C1). This is the plan's own nominated highest-value check.
2. **V1/V2/V3 assert metadata, not behaviour.** They read the report's `Command line parameters` line. That line is written from `ReportConfig.adapters_r2` (`main.rs:2290-2298`), which is a *sibling copy* of the `TrimConfig.adapters_r2` the trimmer uses (`setup_trimming` builds both from the same local at `main.rs:944-951`). Close enough that a mis-thread is unlikely — but "the report echoes the right sequence" is not "the reads were trimmed with it". Add the C2 md5 oracle, or at minimum assert the R2 `Reads with adapters` count moves (4,882 → 5,444 on this fixture).
3. **No row for the new hard-error surface** (I4): `-a2 ''` and `-a2 'ACGT!'` on a preset path.
4. **No row for `--consider_already_trimmed` + `-a2`** (I5).
5. **No row for specialty-mode `-a2`** (I6) or for the SE display/warning contradiction (I7).
6. **V12 has no expected `-n` interaction.** Repeated `-a2 S1 -a2 S2` sizes the R2 per-adapter stats vector via `TrimConfig::r2_adapter_count()` (`trimmer.rs:49`), which switches from `adapters.len()` to `adapters_r2.len()` once R2 is non-empty. Post-fix, `--illumina -a2 S1 -a2 S2` gives R1 one adapter and R2 two — a shape reachable today only via `-a`. Assert the R2 report's per-adapter breakdown has two entries; an off-by-one there would panic or mis-attribute silently.
7. **CI Step 6 covers the least common path.** `--illumina -a2 SEQ` is the *reported* case, but §2.2's own analysis says the auto-detect row is "the most likely to affect others". Perl 0.6.11 accepts `--paired -a2 SEQ` with auto-detection just as happily (the `die` at `trim_galore:3131` fires on `-a` + preset, not `-a2` + preset — *verified*), and both implementations detect Illumina on this fixture. Add the auto-detect invocation, or prefer it.

One thing the plan should also state: V4's "identical output to `dev` before the fix" is the *only* pre/post byte-identity check that survives C1 on a `-a2` invocation, so it is load-bearing. Keep it.

---

## 5. Alternatives

1. **Apply the override in `setup_trimming` instead of restructuring `resolve_adapter`** — I3. Lower risk, smaller diff, no signature change. My recommendation.
2. **Named struct instead of the 5-tuple.** If §4/§5 are kept as written: I checked whether `clippy::type_complexity` would fire on the proposed alias, since CI is `-D warnings`. It does **not** — I built a scratch crate outside the repo with the exact alias and `cargo clippy --all-targets -- -D warnings` was clean. So §5 is CI-safe and the choice is pure style. A `struct AdapterResolution { label, r1, r2, autodetect_poly_g, displaced_r2_default }` is still worth the six lines given two of the five fields are `String`-shaped and one call site destructures positionally — but it is not a blocker. Adopting I3 makes the question moot.
3. **Move `resolve_adapter` into the library** — I8. Testability, not behaviour.
4. **Hard error for SE `-a2`** (Perl parity). Already booked as open item §10.1; I agree with the warning, but note that I4's empty-`-a2` case means the plan is not actually holding the "no new failures" line consistently. Pick one principle and apply it.
5. **Treat empty/whitespace `-a2` as absent** rather than as an error — the narrower fix for I4, and it keeps `--illumina -a2 ''` working.

---

## 6. Action items

### Critical

1. **Rewrite §3.5/A1/§2-line-20 and re-specify V7.** R1 *trimming* is invariant; the R1 *val file* changes on V1–V3 through the joint pair-length filter (139 vs 4 pairs dropped, measured). Split into "R1 cutadapt-stage stats unchanged" (all rows) and "R1 file byte-identical" (V4/V5/V6/V13 only). Update §7 to say both reads' output changes for affected invocations.
2. **Add the differential oracle as the primary validation.** Post-fix `--illumina -a2 AAATCAAAAAAAC` must equal today's `-a AGATCGGAAGAGC -a2 AAATCAAAAAAAC` — R1 `a571f332bad6343d62f1ec37120e1a27`, R2 `1504b7c1f2c34d9d605b4db702537213`. Same for the auto-detect path. Record in the plan that the sequence discriminates (R2 48.8% → 54.4% adapter hits) so V15 is known to be satisfiable.

### Important

3. **Adopt the caller-side override** (I3): kills the Medium risk in §11 and the §5 signature change; ~10 lines, zero existing branches touched.
4. **Decide and document the empty/malformed `-a2` behaviour** (I4), fix the ambiguous §3.5 bullet, and parse `-a2` before the 1 M-read auto-detect scan so bad input fails fast.
5. **Enumerate `--consider_already_trimmed` + `-a2`** (I5) and stop "Only quality trimming will be carried out" from being printed when it is false.
6. **Widen §3.4's warning to specialty modes** (I6) — `--paired --clock -a2 SEQ` is currently a silent discard with no warning under the plan as written.
7. **Resolve the SE display/warning contradiction** (I7) by gating `Adapter 2 (Read 2):` on `cli.paired`.
8. **Move `resolve_adapter` into the library** (I8) so the §3.2 matrix is unit-testable across all seven branches; six of them need no fixture I/O.
9. **Add V-rows** for I4, I5, I6, I7 and for the two-`-a2` per-adapter stats shape (§4 item 6 above).
10. **Prefer or add the auto-detect invocation to CI Step 6** — Perl accepts `--paired -a2 SEQ`, and §2.2 identifies that row as the widest-impact one.

### Optional

11. Fix "v2.0.0" → "v2.1.0" in §7 and §10.2 before it reaches the CHANGELOG; there is no v2.0.0 tag.
12. Add `parameters.adapters_r2` in the JSON report to §7's list of visible changes (fifth item).
13. One line in §6 qualifying "no measurable change": a long `-a2` (e.g. the Perl `--polyA` idiom `T{150}`) exceeds the 64 bp Myers prefilter limit and runs the full DP per R2 read.
14. One line noting U1 (Perl's `--polyA` overwrites `-a2` and applies it 5') and U2 (Perl keys its R2 defaults off the resolved sequence, not the flag) as parts of the 0.6.11 contract this fix deliberately does not restore.
15. Name the docs targets for Step 8: `docs/src/content/docs/guide/adapters.md:31` and `:38` show `--paired -a2 …` examples with **no** `-a`, which are no-ops today; the flag table at `:78` has no caveat.
16. Note that `--output-format ubam` inherits the fix via `setup_trimming` (`main.rs:1806`) — verified working, no code needed.
17. V13 should name the specific `validation` job steps it re-runs; "the five" undercounts them.
