# Code Review B — `fix/polya-behind-polyg` (#504)

**Reviewer:** B (independent, fresh context)
**Date:** 2026-10-20
**Target:** `git diff dev..HEAD` on `fix/polya-behind-polyg` — `b3ad41a`, `23cb2a6`, `869912c`
**Spec:** `plans/10202026_polya-behind-polyg/PLAN.md` rev 2
**Mode:** report only — no fixes applied, nothing under `src/`, `tests/`, `docs/`, `.github/` or `CHANGELOG.md` was touched.

## Verdict

**The code change is correct and I found nothing wrong with it.** The reorder is a one-statement
move, provably inert unless `--poly_a` is passed, the stats stay attributed per stage, and
`bp_after_cutadapt` is untouched. I verified the fix end-to-end against the pre-change binary and
swept three real fixtures for over-trim.

The defects I found are in the **documentation commit** and in the **tests**. One is High: a
user-facing factual claim about Perl v0.6.x in `migration.md` is wrong, and wrong in the opposite
direction from what rev 2 of the plan "corrected". Everything else is Medium or below.

---

## What I verified by running something

| # | Check | Result |
|---|---|---|
| V1 | Full suite: `cargo test --release` | **516 lib + 236 integration, 0 failed** (plan expected 509 lib; +7 is the new unit tests) |
| V2 | `cargo fmt --all -- --check` | clean |
| V3 | `cargo clippy --all-targets --release -- -D warnings` | no output — **but a warm-cache run (0.28 s, nothing recompiled)**, so treat as indicative, not verified |
| V4 | The reported bug, end to end, dev binary vs branch binary on the plan's synthetic set | **dev:** `flush_A40` 60, `A40_G15` **100**, `A20_G20` **80**, `A60_G20` **120**, `no_tail` 60. **branch:** all five → 60, sequence byte-equal to the insert. Fix confirmed; the integration test genuinely fails on `dev`. |
| V5 | Over-trim differential sweep (Validation 3), `--poly_a --poly_g --length 0`, dev vs branch | `illumina10K_with_polyA` 10 shortened / 2 longer / worst cut 16 bp; `BS-seq_10K_R1` 7 / 9 / 15 bp; `SRR24766921_RRBS_R1` 15 / 5 / 8 bp; `PolyA.fastq.gz` 0 / 0. **Every shortened output is a prefix of the old output** — no mid-read divergence anywhere. Matches the plan's 12-read footprint. |
| V6 | `bp_after_cutadapt` / #232 invariant | `Total written (filtered)` **identical** in all three fixtures (1,435,389 / 644,525 / 490,527 bp). |
| V7 | Stat attribution shift | poly-A up (2129→2139, 402→409, 113→128), poly-G **down** (217→215, 230→221, 486→481; bases 888→882, 970→938, 2614→2587). See L6. |
| V8 | R2 prefix-stability claim, `polyAT_R1/R2`, 13,000 pairs | Holds for every record. But see M2 for what the test actually measures. |
| V9 | Perl v0.6.10 `--polyA`, read from `git show 0.6.10:trim_galore` | See H1 — the adapter was **`A{20}`**, not 10 A's. |
| V10 | Poly-G auto-detection in the `--2colour` test | `Poly-G trimming: not enabled (auto-detection found 2 of 2 reads (100.00%) — below threshold)`. See L3. |
| V11 | Does clippy flag `&PathBuf` params? | Reproduced the two signatures in a standalone crate: **no lint**. L5 is style only, not a CI risk. |
| V12 | `--poly_a` can only come from the CLI | `grep` for `poly_a: true` across `src/ tests/ examples/` → nothing; the only sources are four `poly_a: cli.poly_a` sites in `main.rs`; no `impl Default for TrimConfig`. The reorder is a strict no-op without the flag. |

---

## Critical

None.

---

## High

### H1 — `migration.md:58` states the wrong Perl adapter, and calls the correct `--help` text a fallback

`docs/src/content/docs/reference/migration.md:58`:

> The adapter was auto-detected per sample as `AAAAAAAAAA` (10 A's) or the poly-T equivalent. The
> `A{20}` form named in the v0.6.x `--help` text was only the fallback used when neither probe
> matched anything.

Both halves are false. From `git show 0.6.10:trim_galore`, lines 501–530:

```perl
if ($polyA){ # specialised PolyA trimming
    ($adapter,$adapter_name) = autodetect_polyA_type();   # returns 'AAAAAAAAAA' or 'TTTTTTTTTT'
    if ($validate){                                       # paired-end
        if ($adapter =~ /A+/){
            $adapter = extend_adapter_sequence("A",20);   # <-- unconditionally overwritten
            $a2      = extend_adapter_sequence("T",150);
        } elsif($adapter =~ /T+/){ ... }
    } else {                                              # single-end
        if ($adapter =~ /A+/){ $adapter = extend_adapter_sequence("A",20); }
        elsif($adapter =~ /T+/){ $adapter = extend_adapter_sequence("T",20); }
    }
}
```

`extend_adapter_sequence` is `return (${letter}x$number)` (`:3113`), so the string Cutadapt received
was **20 A's** (or 20 T's), with `T{150}` / `A{150}` as `-a2` for the 5' end of Read 2. The 10-mer
returned by `autodetect_polyA_type()` is used **only to decide A vs T** — it is a probe, never the
adapter. `extend_adapter_sequence('A',20)` inside `autodetect_polyA_type()` (`:2376`) is the
inconclusive-detection branch and is **also 20**, not 10.

So the `--help` text ("a default of A{20} for Read1 and T{150} for Read2", `:3357`) and the
validation message (`:2754`) are both **accurate**, and the code does not disagree with them. The
claim that they disagree originated in `PLAN_REVIEW_B.md` §3a, was adopted as a correction in PLAN
rev 2 ("rev 1 got all three wrong"), and has now been published in the docs. Rev 1 was right on this
point.

This matters beyond the sentence: PLAN rev 2's **Decision** section argues that reverting to a
homopolymer adapter would be worse than rev 1 claimed because "`A{10}` (what Perl actually used) cuts
6,401 reads / 460,805 bp against `A{20}`'s 3,818 / 128,130". The premise is wrong — `A{20}` *is*
what Perl used. The decision to keep the Cutadapt scan stands on its own merits (it matches modern
Cutadapt, and the mode isn't byte-comparable anyway), so **no code change follows from this**. Only
the docs and the plan need correcting.

**Exact change I would make** — replace line 58 with:

```markdown
- The adapter was `A{20}` — twenty A's — or `T{20}`, with `T{150}` / `A{150}` as `-a2` for the 5' end of Read 2. Auto-detection only chose between A and T, by counting whether `AAAAAAAAAA` or `TTTTTTTTTT` occurred in more of the first million reads.
```

And, since this now becomes genuinely useful migration advice, consider opening the `-a 'A{N}'`
subsection (line 65) with it:

```markdown
`-a 'A{20}'` reproduces exactly what v0.6.x `--polyA` passed to Cutadapt for single-end data.
```

The other two Perl claims in that list are **correct** and I verified them:

- "carried no `-q` and no `--length` of its own" — the `elsif($polyA)` dispatch (`:1130`) invokes
  `cutadapt $cores -e $error_rate -O $stringency -a $adapter` with nothing else (R2 uses `-g $a2`).
- the `32:A:` prefix and `_PolyA:32` suffix — `$l1 =~ s/\@/\@${diff}:A:/` and `$l1 .= "_PolyA:$diff"`
  at `:1211-1212`, guarded by `if ($diff > 0)`.

---

## Medium

### M1 — `poly.md`'s recommendation table still prescribes the workaround this change removes

`docs/src/content/docs/modes/poly.md:62` adds:

```markdown
| 3'seq on a 2-colour instrument | Add `--poly_a` and use `--2colour 20` rather than `--quality 20`. |
```

and line 54 says "If poly-A tails survive on 2-colour data, use `--2colour 20` in place of
`--quality 20`."

That was the right advice yesterday. After this fix, `--poly_a` plus auto-detected `--poly_g` handles
`[insert][poly-A][poly-G]` on its own — which is the whole claim of the CHANGELOG entry two files
over. A 3'seq user who reads the release note ("now removes poly-A tails that sit behind a 3'
poly-G run") and then lands on the table is told they still need `--2colour`. The two pages
contradict each other.

`--2colour` is still the right fallback for the two residual cases (`--no_poly_g`, and the mirror
arrangement), and the prose a few lines above already frames those. The table row should match.

**Exact change:**

```markdown
| 3'seq on a 2-colour instrument | Add `--poly_a`. Poly-G is auto-detected and removed first, so a tail behind the G-run goes too. |
```

and reword line 54 to attach the `--2colour` advice to the residuals rather than to 3'seq generally:

```markdown
In either of those two cases, `--2colour 20` in place of `--quality 20` removes the G-run during quality trimming, ahead of both poly stages.
```

### M2 — the paired R2 test's anti-triviality guard does not guard the thing it claims, and it silently drops 23% of the records

`tests/integration_poly.rs::paired_r2_5prime_clip_is_independent_of_3prime_trimming`. Two separate
problems; both measured.

**(a) The `shortened > 0` guard is satisfied almost entirely by records the assertion is vacuous
for.** I instrumented the test's own two output files plus the input fixture:

```
R2 records compared:                         10,012
3'-shortened by poly-G (`shortened`):           215
5'-clipped by poly-T at all:                    129
BOTH shortened and 5'-clipped:                    5
```

The assertion under test is "the 5' clip did not move". For a record with **no** 5' clip in either
run, `seq.starts_with(with)` is true because the output is a prefix of the input in both runs — it
holds whatever the 5' stage does, because the 5' stage did nothing. Only the **5 records** in the
intersection exercise the claim. The guard counts 215, so it reports "this proves something" on the
strength of 210 records that prove nothing. It is non-zero by luck of the fixture, and nothing pins
that it stays non-zero.

**Exact change** — count the intersection and assert on that instead:

```rust
    let mut both = 0;
    for (id, seq) in &untrimmed {
        let with = trimmed.get(id).expect("every record in both runs");
        assert!(
            seq.starts_with(with.as_str()),
            "{id}: R2 5' clip moved — {with} is not a prefix of {seq}"
        );
        if with.len() < seq.len() && with.len() + seq_offset(id) < seq.len() {
            both += 1;
        }
    }
    assert!(
        both > 0,
        "no R2 record is both 3'-trimmed and 5'-clipped, so the prefix check proves nothing"
    );
```

(`seq_offset` needs the input fixture to tell a 5' clip from a 3' trim; simpler is to read
`polyAT_R2.fastq.gz` once and count records where `!original.starts_with(&untrimmed[id])` **and**
`trimmed[id].len() < untrimmed[id].len()`. Either way the constant to assert is 5 on today's
fixture — I would assert `>= 1` rather than `== 5` so a fixture refresh doesn't break it.)

**(b) `records()` returns a `HashMap` keyed by the header line, and the fixture has duplicate read
IDs.** The run processes **13,000** pairs; the map holds **10,012** entries. 2,988 records (23%) are
silently collapsed. `assert_eq!(trimmed.len(), untrimmed.len())` passes because both sides collapse
identically, and the surviving value is the same input record in both runs, so the comparison is not
*wrong* — it is just a quarter smaller than it looks, with nothing saying so. Keying on
`(line_index, header)` or returning a `Vec` and zipping (the records are in input order in both runs,
and `Pairs removed: 0`) would compare all 13,000.

### M3 — comment policy: the new test comments are the exact shape `CLAUDE.md` rules out

`CLAUDE.md` → *Code comments*: one line by default, two maximum with a reason; state the current
fact; never justify the change in the source. It says explicitly that this "covers code, tests and
config that gets committed". Four places in the diff exceed it:

1. **`tests/integration_poly.rs:1-10`** — a ten-line `//!` header that restates the bug, the
   algorithm, the reason for choosing integration over unit tests, and the testing philosophy. This
   is the "pre-empting *but doesn't this break X?*" pattern the rule names; the argument belongs in
   `b3ad41a`'s commit message, which already carries it. Suggested replacement, the whole header:

   ```rust
   //! #504 — `--poly_a` removes a poly-A tail that sits behind a 3' poly-G run.
   ```

2. **`paired_r2_5prime_clip_is_independent_of_3prime_trimming`, 5 lines** — carries the gate formula
   and a mechanism derivation. Suggested: `// Poly-T comes off R2's 5' end and poly-G off its 3', so the two stages never compete.`

3. **`test_poly_g_behind_poly_a_survives`, 3 lines** — ends with "pinned so the trade stays visible",
   which is justification. Suggested: `// The mirror arrangement: the tail blocks poly-G, so the G-run stays.`

4. **`reads_without_poly_a_are_untouched`, 3 lines** — "so a regression that over-trims is visible
   here rather than hiding behind a tail-free assertion" is justification, and it is also inaccurate
   (M5). Suggested: `// Without --poly_a both poly stages are skipped and the tails survive.`

The **source** comments are fine. The new two-line insert in `src/trimmer.rs:195-198` states a
current fact with no history and no justification, which is exactly right. It could be one line —
`// Runs ahead of poly-A: a tail behind a terminal G-run is invisible to a 3'-anchored scan.` — but
two is within the rule. The pre-existing #321 block below it is untouched and out of scope.

### M4 — the integration module header makes two claims that are false

Independent of M3's length problem, as long as that header exists it should be true:

- *"only an end-to-end run pins that order"* — not so. I walked the four unit tests in
  `src/trimmer.rs` and `test_poly_a_tail_behind_poly_g_is_trimmed` and
  `test_poly_g_behind_poly_a_survives` both fail under a reverted order (the first would see
  `INSERT_60 + A*a` instead of `INSERT_60`; the second would see 59 bp with `poly_g_trimmed == 16`
  instead of 75 bp with 0). The `trim_read`-level tests pin the order perfectly well; the integration
  tests add end-to-end confidence, not the only signal.
- *"Every case asserts an exact output length"* — true of three cases, false of
  `paired_r2_5prime_clip_is_independent_of_3prime_trimming`, which asserts a prefix relation and a
  count.

### M5 — nothing tests `--poly_g` without `--poly_a`, which is what nearly every user runs

`reads_without_poly_a_are_untouched` passes `--no_poly_g` and no `--poly_a`, so **both** poly blocks
are skipped. As a control it pins "the reorder is inert when neither flag is set", which V12 already
establishes statically. What it does not cover — and nothing else does either — is **poly-G enabled,
poly-A off**: the configuration of every 2-colour user who never heard of `--poly_a`, reached through
auto-detection. That is the path with real users on it, and it has no regression guard.

It is one case in the file you already have. **Exact change** — add to `poly_a_tail_behind_poly_g_is_trimmed`'s
file:

```rust
#[test]
fn poly_g_alone_leaves_the_tail() {
    // --poly_g without --poly_a: the G-run goes, the tail stays.
    let dir = fresh_tmpdir("tg_poly_g_only");
    let input = dir.join("reads.fastq");
    let expected = format!("{INSERT}{}", "A".repeat(40));
    write_fastq(
        &input,
        &[("A40_G15", format!("{expected}{}", "G".repeat(15)))],
    );

    run(&[
        "--quality", "20", "--phred33", "--poly_g", "--dont_gzip",
        "-o", dir.to_str().unwrap(),
        ADAPTER[0], ADAPTER[1], input.to_str().unwrap(),
    ]);

    assert_eq!(sequences(&dir.join("reads_trimmed.fq"))[0], expected);
}
```

I confirmed that expectation against both binaries (V4: `A40_G15` → 100 bp on `dev` with both flags,
which is `INSERT + A40`).

---

## Low

- **L1 — two 60 bp insert constants that differ only in the last base.** `src/trimmer.rs:1182`
  `INSERT_60` ends in `…GAATCG`; `tests/integration_poly.rs:28` `INSERT` ends in `…GAATCT`. I diffed
  them: identical at positions 0–58, differ only at index 59 (`G` vs `T`). Both carry a doc comment
  saying "60 bp of insert", so they read as the same constant. That terminal `G` is load-bearing for
  `test_poly_a_r2_head_is_unaffected_by_the_stage_order`'s `59` and `16`: if anyone edits the last
  base those two numbers silently become wrong in a way that still compiles. Make the coupling
  checked rather than implied, e.g. `assert_eq!(result.poly_g_trimmed, 15 + 1, "the insert's trailing G goes with the run");`
  and `assert_eq!(record.seq, &INSERT_60[..INSERT_60.len() - 1]);`. Giving the integration constant a
  one-clause note that its last base is deliberately not `G` would also help.

- **L2 — the magic `75`.** `test_poly_g_behind_poly_a_survives` asserts `record.seq.len() == 75`
  with no arithmetic shown. `INSERT_60.len() + 15` says the same thing and survives an edit to the
  constant. (`59` and `16` in the R2 test are at least gestured at by the comment above them; `75`
  is not.)

- **L3 — `two_colour_quality_trimming_also_exposes_the_tail` depends on an auto-detection threshold
  it does not state.** The test passes `--2colour 20 --poly_a` and neither `--poly_g` nor
  `--no_poly_g`, so poly-G auto-detection runs. Its claim ("`--2colour` scores G as no-signal, so the
  G-run goes in the quality step") only holds if poly-G is off. I checked: it is off, but only
  because the threshold is `(reads_scanned / 10_000).max(10)` (`src/main.rs:1231`) and the fixture
  has 2 reads — detection reports *100.00% of reads have poly-G tails, below threshold*. Add a third
  read with a G-run and the test stops testing `--2colour`. Pass `--no_poly_g` explicitly; the test
  still passes (verified: `--2colour` alone takes both reads to 60 bp on the pre-change binary too,
  which is also why this case is a documentation test rather than a regression guard — correctly
  noted in its comment).

- **L4 — the paired test's mechanism argument is incomplete.** "the 5' gate `errors * 5 <= i + 1`
  depends on the position only, never on the read length. So however much poly-G takes off the 3'
  end, the 5' clip is the same" — the gate is position-only, but the 5' loop runs
  `.enumerate().take(n)`, so truncating the 3' *does* remove candidate positions `i >= n'`. The
  conclusion survives for a stronger reason: a new 5' maximum beyond the cut point would need the
  trailing T-run to exceed twice the G-run, and that is precisely the condition under which the 3'
  G-scan would not have cut there (its score at that G is `#G - 2·#T_after < 0`). So the two cannot
  both happen. If the comment stays at all (M3), it should say that, or just cite the measurement.

- **L5 — `&PathBuf` parameters.** `write_fastq`, `sequences` and `records` take `&PathBuf`; it is the
  only occurrence in `src/` or `tests/`, and `&Path` is the idiom. I checked it is not a CI risk
  (V11: clippy does not fire on these signatures).

- **L6 — PLAN §Integration's stat claim is wrong; the CHANGELOG line it called for was dropped.** The
  plan says "Bases that are currently counted as poly-G will be counted as poly-G still, but reads
  that currently report poly-G only will now also report poly-A", and "worth one line in the CHANGELOG
  entry". Measured (V7): poly-G **falls** in all three fixtures — 888→882, 970→938, 2614→2587 bp;
  217→215, 230→221, 486→481 reads — because the mirror arrangement moves the blocked stage onto
  poly-G. The shift is 0.02–0.2% of reads, so omitting it from the CHANGELOG is defensible; the plan
  text is what needs the correction, not the entry. No user-facing artifact repeats the claim.

- **L7 — the docs correction stops at the two spots the plan named.** Poly-A is still presented as a
  v2-only capability in `docs/src/content/docs/guide/overview.md:13` ("newer capabilities (e.g.
  `--poly_g` auto-detection, a generic `--poly_a` trimmer) extend beyond the Perl version"),
  `docs/src/content/docs/reference/credits.md:10` ("extra capabilities: … a generic poly-A trimmer"),
  and `README.md:19`/`:30` ("Adds poly-G auto-detection … a generic poly-A trimmer"). The sharpest is
  `overview.md:13`, since "extend beyond the Perl version" is the exact claim `migration.md:28` just
  stopped making. `performance/threading.md:85` ("both built in without external tools") is literally
  true — v0.6.x needed Cutadapt — so I would leave that one. Not blocking; worth a follow-up line in
  the plan so it isn't rediscovered.

- **L8 — @marenmsk is not credited.** The plan's commit-3 item 5 says "one sentence, repo format,
  crediting @marenmsk". The entry has no credit. Having read the file, **no** entry in it credits a
  reporter, so the entry is right and the plan item was the outlier — noting it only so the
  plan-vs-implementation diff is complete, not as something to change.

- **L9 — no uBAM coverage.** `--output-format ubam --poly_a` is reachable (only `--clump_only`
  rejects the poly flags, `src/cli.rs:863-871`) and inherits the reorder through the shared
  `trim_read` — I confirmed `trim_read` is the sole implementation and that
  `trimmer.rs:650`/`:735-736` (the `_to_bam` entry points) and `parallel.rs:495-496,1133` all route
  through it. Nothing tests it. Low, given the single implementation.

- **L10 — `--poly_a` help text, `src/cli.rs:366-367`.** "Runs after adapter and poly-G trimming, so
  tails behind either are also removed." Nothing survives behind an adapter — a 3' adapter match
  already removes everything downstream of itself, which is precisely the point `migration.md:61`
  makes. "so a tail behind a poly-G run is removed too" would be exact and no longer. Pre-existing
  looseness carried forward, not introduced.

- **L11 — temp-directory hygiene.** `fresh_tmpdir` takes a fixed slug, `remove_dir_all`s it, and
  never cleans up. Two concurrent `cargo test` runs on this machine collide (which is live right now,
  with two reviewers in one tree). Eight of the seventeen integration test files already suffix
  `std::process::id()`; the other nine, including `integration_partial_output.rs` which this file was
  modelled on, do not. Repo convention is split, so this is a note on which half to follow rather
  than a deviation.

---

## Answers to the specific questions asked

1. **Is the stage order correct for every combination?** Yes. The reorder is a statement move inside
   `trim_read`, which V12 confirms is the single implementation — specialty modes (`--hardtrim5/3`,
   `--clock`, `--implicon`) bypass it entirely and do no poly trimming (`grep poly src/specialty.rs`
   is empty), and the uBAM and parallel paths call into it. R1/SE is the only place the two stages
   compete; on R2 poly-G takes the 3' end and poly-T the 5', and V8 confirms 13,000 R2 records come
   out prefix-identical. `--rrbs` is step 2.5, ahead of both, and its relative position is unchanged.
   Stats are computed locally inside each block from that block's own before/after length, so
   attribution per stage stays exact; the aggregate movement in V7 is real behaviour, not
   mis-attribution. `bp_after_cutadapt` is captured at `trimmer.rs:161`, ahead of RRBS and both poly
   stages, untouched by the diff and confirmed byte-stable in V6.

2. **Comment quality.** Source comments pass. Test comments do not — M3, four places, with
   replacements.

3. **Do the tests pin what they claim?** Mostly yes, with one real weakness.
   Order-sensitive (would fail on a revert): `test_poly_a_tail_behind_poly_g_is_trimmed`,
   `test_poly_g_behind_poly_a_survives`, and `poly_a_tail_behind_poly_g_is_trimmed` (empirically
   confirmed — V4 shows `dev` produces 100/80/120 where the test demands 60).
   Order-**in**sensitive by design, correctly so: `test_poly_a_tail_flush_with_3prime_end_is_trimmed`,
   `test_stage_order_does_not_touch_reads_without_tails`,
   `test_no_poly_g_leaves_tails_behind_a_g_run`, `test_poly_a_r2_head_is_unaffected_by_the_stage_order`,
   `two_colour_quality_trimming_also_exposes_the_tail`, `reads_without_poly_a_are_untouched`.
   No assertion is self-satisfying in the sense of holding regardless of the code under test — except
   the prefix assertion in the paired test for 210 of its 215 guarded records, which is M2.
   The new `src/quality.rs::test_poly_a_short_reads` is sound: I checked all four assertions by hand
   across `n = 3..=10` and each is a real consequence of the scan (`0` for all-A, `n` for all-T at the
   5', `n`/`0` for the tail-free control), with the control's first two bases non-T so the 5' score
   never goes positive.

4. **Test brittleness and magic numbers.** L1 (the terminal base of the shared constant, measured),
   L2 (`75`), L3 (the auto-detection threshold the `--2colour` case silently relies on).

5. **Docs claims about Perl v0.6.x.** H1 — the `AAAAAAAAAA` claim is wrong and the "`A{20}` is only a
   fallback" claim is wrong. The `-q`/`--length` claim and the `32:A:` / `_PolyA:32` claim are both
   correct, verified against `0.6.10:trim_galore` lines 1130, 1211-1212 and 3352-3375. Internal links
   resolve: `/modes/poly/` → `docs/src/content/docs/modes/poly.md`, `/reference/migration/` →
   `reference/migration.md`, both matching the repo's existing root-relative-with-trailing-slash
   style (`guide/flags.md:43`, `credits.md:10`). The `#stage-order` anchor matches `### Stage order`.
   Note that nothing enforces this — `astro.config.mjs` loads no link validator, so `docs-build`
   would not catch a broken internal link.

6. **CHANGELOG.** In line with the file: the two existing `Unreleased` neighbours are bold headline
   plus one to three sentences, and the new entry is shorter than both. The `#504` issue link and the
   `https://trimgalone.com`-style doc link both match the file's conventions (no `www`, as every other
   entry uses, even though `astro.config.mjs` sets `site: 'https://www.trimgalore.com'`). The anchor
   is real. Against the global one-sentence rule it is two sentences and carries the mechanism
   ("because the poly-A scan is anchored at the 3' end and a G-run in front of the tail blocked it") —
   trimmable to *"Poly-G is trimmed first; see [the poly-X docs](…) for the order and its two edge
   cases."* if you want it to meet the rule rather than the file.

---

## Fixes applied

None, per instruction. Every change above is given as the exact edit.

## Recommended order

1. **H1** — `migration.md:58`. It is published, specific, checkable, and currently wrong.
2. **M1** — the `poly.md` table row, which contradicts the release note.
3. **M5, M2** — the two test gaps, in that order: M5 is five lines and covers the path with users on
   it; M2 is the subtler one.
4. **M3, M4** — the comment pass; M4 disappears if M3's header replacement is taken.
5. **L1–L3** — cheap, and L1 is the one most likely to bite a future editor.
6. L4–L11 at leisure; L6 and L7 want a line in the plan rather than a code change.
