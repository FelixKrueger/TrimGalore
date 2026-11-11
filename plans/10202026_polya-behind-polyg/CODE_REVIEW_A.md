# Code review A — `fix/polya-behind-polyg` (issue #504)

**Reviewer:** A (independent, fresh context)
**Target:** `git diff dev..HEAD` — `b3ad41a`, `23cb2a6`, `869912c`
**Date:** 2026-10-20
**Mode:** report only, no fixes applied (tree shared with Reviewer B). No file outside this report was touched.

---

## Summary

**The code change is correct and I verified it end to end.** The reorder is the minimal fix, it reaches
every entry point through the single `trim_read`, the stats stay attributed to the stage that did the
work, `bp_after_cutadapt` is captured upstream of both poly stages and is unaffected, and the real-data
footprint is 12 of 10,900 reads — matching the plan's measurement exactly, reproduced independently.
The Perl byte-identity invariant is untouched: no CI validation case passes `--poly_a`.

**The problem is in the docs commit.** `migration.md` now asserts that Perl v0.6.x `--polyA` trimmed
with `AAAAAAAAAA` (10 A's) and that the `A{20}` in the upstream `--help` was only a fallback. The Perl
source says the opposite: the 10-mer is a direction probe that is **overwritten two lines later**, and
`A{20}` / `T{150}` is what Cutadapt actually received on every run. The PR publishes a correction that
is itself wrong, and tells readers the upstream help text was misleading when it was accurate. Same
bullet also states Perl's polyA mode carried "no `--length`" — it did (`0.6.11:trim_galore:1307`).

Secondary: one integration test can pass vacuously; the CHANGELOG drops the reporter credit the plan
called for; the added source comment is a 3-line paragraph where the plan said one line.

**Verdict: do not merge the docs commit as written.** `b3ad41a` and `23cb2a6` are sound modulo the
test and comment nits below.

---

## What I verified by running something (vs. reasoning about)

| # | Check | Result |
|---|---|---|
| V1 | `cargo test --release`, full suite | **516 lib + 4 `integration_poly` + all others pass, 0 failed** |
| V2 | `cargo clippy --all-targets --release -- -D warnings`, forced recompile (`touch src/lib.rs tests/integration_poly.rs`) | **clean** — the `&PathBuf` params do **not** trip `ptr_arg` (they are passed to generic `AsRef<Path>` sinks, which the lint skips) |
| V3 | Pre-change (`/tmp/claude-501/tg_dev_b1/…`) vs post-change binary on the synthetic 5-read case | dev leaves 100/80/120 bp; HEAD gives 60 bp for all five. **The headline test would fail on a revert.** |
| V4 | Both binaries on `test_files/illumina10K_with_polyA.fastq.gz`, `--poly_a --poly_g --length 0` | 10,900 reads, **12 changed**: 10 shorter by 3–16 bp (tail now removed), 2 longer by 3 bp (G-run now retained). Plan's "12-read footprint" reproduced. |
| V5 | Both binaries, `--paired --poly_a --poly_g` on `polyAT_R{1,2}` | R1: 12 of 13,000 changed. **R2: 0 of 13,000 changed** — the R2-independence claim holds against the real pre-change build, not just across flag settings. |
| V6 | Does poly-G auto-detection fire in `two_colour_quality_trimming_also_exposes_the_tail`? | **No.** `threshold = (reads_scanned / 10_000).max(10)` → 2 reads never clears a floor of 10; binary prints "not enabled … below threshold". The test genuinely isolates the `--2colour` path. |
| V7 | Perl claims against `0.6.11:trim_galore` (tag fetched) and `0.6.9` | see F1/F2 — both versions identical on the relevant lines |
| V8 | `#504` author via `gh issue view` | `marenmsk` |
| V9 | Doc anchors and link targets | `id="stage-order"` present in `docs/dist/modes/poly/index.html`; `/reference/migration/` and `/modes/poly/` match the convention at `guide/flags.md:43`, `reference/credits.md:10`. No stale references to the replaced `## Note on -a 'A{N}'` heading. |
| V10 | CI exposure | `.github/workflows/ci.yml` uses `--poly_a` nowhere; the only poly flag is `--no_poly_g` at `:968`. Validation matrix unaffected. |

Everything else below is reasoning about code I read.

---

## Findings

### F1 — HIGH — `migration.md` inverts the Perl adapter claim and contradicts the upstream `--help`

`docs/src/content/docs/reference/migration.md:58`:

> The adapter was auto-detected per sample as `AAAAAAAAAA` (10 A's) or the poly-T equivalent. The
> `A{20}` form named in the v0.6.x `--help` text was only the fallback used when neither probe matched
> anything.

Both halves are wrong. `0.6.11:trim_galore:512-541` (identical in `0.6.9`):

```perl
512:  if ($polyA){
513:      ($adapter,$adapter_name) = autodetect_polyA_type();
514:
515:      if ($validate){                                     # paired-end
516:          if ($adapter =~ /A+/){
517:              $adapter =  extend_adapter_sequence("A",20);   # <-- overwritten
518:              $a2 = extend_adapter_sequence("T",150);
...
530:      else{                                               # single-end
532:          if ($adapter =~ /A+/){
533:              $adapter =  extend_adapter_sequence("A",20);   # <-- overwritten
```

`autodetect_polyA_type()` returns the 10-mer, but the caller immediately **discards it as a value and
keeps only its direction**: the `=~ /A+/` test picks A-vs-T, then `$adapter` is reassigned to
`extend_adapter_sequence("A",20)` — 20 A's — in *both* branches, unconditionally. The 10-mer never
reaches Cutadapt. The trimming adapter on every `--polyA` run was `A{20}` for R1/SE and `T{150}`
5'-anchored (`-g`) for R2, exactly as the `--help` at `:3355` says. The "fallback when neither probe
matched" path returns `extend_adapter_sequence('A',20)` too, so even that case is `A{20}`.

Why this matters beyond the sentence:

1. The doc now tells readers the upstream help text was inaccurate when it was correct. That is worse
   than saying nothing.
2. `PLAN.md` §Decision rests on the same premise — "`A{10}` (what Perl actually used) cuts 6,401 reads
   / 460,805 bp against `A{20}`'s 3,818 reads / 128,130 bp … 3.6× the bases". The Perl-parity baseline
   is the `A{20}` column, not the `A{10}` one. The conclusion probably survives (3,818 is still 2.5×
   today's 1,537), but the stated margin is 3.6× too generous, and anyone who later revisits "restore
   Perl semantics" will measure against the wrong sequence.

Nothing in the shipped code depends on this — `b3ad41a` is a pure reorder — so it is a docs-only fix.

**Change I would make** (replace the bullet):

```markdown
- The adapter was `A{20}` — 20 A's — for Read 1 and single-end, and `T{150}` anchored at the 5' end
  (`-g`) for Read 2. Auto-detection per sample used a 10-mer probe (`AAAAAAAAAA` vs `TTTTTTTTTT`) only
  to choose the direction; the sequence handed to Cutadapt was always the 20-mer.
```

### F2 — MEDIUM — `migration.md` says Perl's polyA mode carried no `--length`. It did.

Same paragraph, `migration.md:57`:

> It carried no `-q` and no `--length` of its own, so it expected input that had already been adapter-
> and quality-trimmed.

The `-q` half is right — `0.6.11:trim_galore:1196` invokes
`cutadapt $cores -e $error_rate -O $stringency -a $adapter`, with no `-q` and no `-m`. But Trim Galore
has never passed `-m` to Cutadapt on *any* path; the wrapper does length filtering itself, and the
polyA branch does it too:

```perl
1307:  if (length $seq < $length_cutoff){
1308:      ++$too_short;
1309:      next;
```

`$length_cutoff` defaults to 20 (`:580`, or 18 for smallRNA at `:556`), and paired-end polyA runs get
the same cutoff in the validation step at `:2002`. So `--length` behaved identically in polyA mode and
in normal mode. The clipping options (`clip_r1`, `three_prime_clip_r1`) and `--max_length` are applied
in that branch as well.

**Change I would make:**

```markdown
- The Cutadapt pass carried no `-q`, so it expected input that had already been adapter- and
  quality-trimmed. The documented workflow was two runs: trim adapters and qualities first, then
  `--polyA`. (`--length`, `--max_length` and the `--clip_r1`/`--three_prime_clip_r1` options still
  applied, as on every other path.)
```

### F3 — MEDIUM — the read-ID bullet omits the change that hits every read

`migration.md:59` records the `32:A:` prefix and `_PolyA:32` suffix, which only appear on reads that
were actually trimmed (`$diff > 0`). The line above it in the Perl source runs unconditionally:

```perl
$l1 =~ s/\s+/_/g;   # removing white spaces from readID
```

Every read in polyA mode came out with whitespace collapsed to `_`, trimmed or not — the upstream help
calls this out first ("white spaces are replaced with `_`"). For the bullet's own audience ("a pipeline
that greps read IDs"), that is the more consequential difference, because it changes the ID of 100% of
records rather than of the trimmed subset, and v2 does neither.

**Change I would make** — append one clause: `…to record how many bases went, and replaced whitespace
in every read ID with `_` whether or not anything was trimmed.`

Minor, same bullet: `32` reads as a literal. It is the per-read base count (`$diff`), carried over from
the help text's worked example. `prefixing `<n>:A:` and appending `_PolyA:<n>`` removes the ambiguity.

### F4 — MEDIUM — `two_colour_quality_trimming_also_exposes_the_tail` can pass on an empty file

`tests/integration_poly.rs:134`:

```rust
for seq in sequences(&dir.join("reads_trimmed.fq")) {
    assert_eq!(seq, INSERT);
}
```

No count assertion. If the output were empty — both reads dropped by the length filter, a writer that
never flushed, an output-naming change — the loop body never executes and the test reports `ok`. The
file's own module doc sets the standard this one misses: *"Every case asserts an exact output length,
not just the absence of a tail."* `poly_a_tail_behind_poly_g_is_trimmed` does it right at `:95`
(`assert_eq!(out.len(), 5, …)`).

I confirmed by hand that the test is doing real work today (V6: poly-G auto-detection is below
threshold, so the `--2colour` quality step is what exposes the tail, and the output is 2×60 bp). The
issue is that nothing holds it to that.

**Change I would make:**

```rust
    let out = sequences(&dir.join("reads_trimmed.fq"));
    assert_eq!(out.len(), 2, "both reads must be written");
    for seq in out {
        assert_eq!(seq, INSERT);
    }
```

While there: adding `"--no_poly_g"` to that invocation would make the isolation structural rather than
a consequence of the 10-read auto-detection floor, which a future fixture of ≥110 reads would silently
break. (I checked the result is unchanged: quality-trim removes the G-run, poly-A removes the tail,
output is `INSERT`.) That is also the idiom the sibling file already uses — `SKIP_PRESCANS` at
`tests/integration_partial_output.rs:34`.

### F5 — MEDIUM — the added comment is a 3-line paragraph where the plan specified one line

`src/trimmer.rs:198-200`:

```rust
    // Runs ahead of poly-A: the G-run is terminal, and a poly-A tail behind
    // one is invisible to a 3'-anchored scan.
```

`PLAN.md` commit-1 item 2: *"Add **one line** above the poly-G block giving the invariant (current
state, no history, per CLAUDE.md 'Code comments')."* CLAUDE.md: *"Default to a single line. Two is the
maximum and needs a reason."* It is also closer to justifying the change than to stating a fact — the
"because a 3'-anchored scan can't see past it" half is the commit-message argument.

It does not reference the previous state, so it is not the historical-narrative failure; it is purely
length.

**Change I would make** — one line, 94 characters with indent:

```rust
    // Runs ahead of poly-A: a tail behind a terminal G-run is invisible to a 3'-anchored scan.
```

The 10-line `#321` block below it is pre-existing and was moved verbatim as the plan required; I am not
flagging it, though the reorder is the natural moment to cut it if you ever want to.

Two test comments are over the same limit and lean the same way — `src/trimmer.rs:1222-1224`
("…pinned so the trade stays visible" — 3 lines, the last clause justifies the test's existence) and
`tests/integration_poly.rs:144-146` ("…so a regression that over-trims is visible here rather than
hiding behind a tail-free assertion" — 3 lines, second half is test-design rationale). Both would read
the same at two lines. The 5-line comment at `:189-193` is discussed in F7.

### F6 — LOW — the CHANGELOG entry drops the reporter credit the plan asked for, and is two sentences

`PLAN.md` commit-3 item 5: *"`CHANGELOG.md` — one sentence, repo format, crediting @marenmsk."* The
shipped entry credits no one. `marenmsk` is confirmed as #504's author (V8), and the file has an
established shape for this — `CHANGELOG.md:621` reads
`([#369](…), reported by\n  @MathieuUm)`, bare handle, inside the issue parenthetical.

The entry is also two sentences against CLAUDE.md's one (the neighbouring `#457` entry runs to four, so
the file's practice has drifted, but the rule is the rule). And the URL is `https://trimgalore.com/…`
while `astro.config.mjs` sets `site: 'https://www.trimgalore.com'` and `docs/dist/CNAME` is
`www.trimgalore.com`; six older entries use the `www.` form. I could not resolve either host from the
sandbox to confirm whether the apex redirects, so treat this as a consistency point, not a broken link.
The `#stage-order` fragment itself is real (V9).

**Change I would make:**

```markdown
- **`--poly_a` now removes poly-A tails that sit behind a 3' poly-G run**
  ([#504](https://github.com/FelixKrueger/TrimGalore/issues/504), reported by @marenmsk). Poly-G is
  trimmed first; see [the poly-X docs](https://www.trimgalore.com/modes/poly/#stage-order) for the
  stage order and its two edge cases.
```

### F7 — LOW — `paired_r2_5prime_clip_is_independent_of_3prime_trimming`: the guard works, the comment overstates, and it is not a revert detector

You flagged this one; here is what I found.

**The prefix assertion does establish 5'-clip independence**, with one theoretical hole. If the two
runs clipped different amounts off the 5' end, the sequences would start at different offsets in the
original read and `seq.starts_with(with)` would fail — *unless* the extra clipped bases happen to be a
repeat of what follows, which inside a poly-T head is not impossible. In practice a divergence requires
the poly-T head's optimum to lie beyond the point poly-G truncated, i.e. a read that is almost entirely
T with a terminal G-run; I tried to construct one and every candidate I traced
(`T*20+G*4`, `T*10+G*20+T*2+G*3`, `T*5+G*3`) ends identically under both orders, because the 3' G-run
and the 5' T-head are disjoint whenever both are shorter than the read.

**The anti-vacuity guard is sufficient.** `shortened > 0` proves poly-G actually trimmed R2 3' ends in
the fixture, so the two runs differ and the prefix check is comparing two genuinely different outputs.
`!trimmed.is_empty()` and the length equality cover the empty/mismatched-set cases, and
`.expect("every record in both runs")` catches an ID that appears in only one.

**But it is not a revert detector.** Under the old order the 5' poly-T clip was computed on the
full-length read and poly-G then truncated the 3' end, so the prefix relation and `shortened > 0` both
hold there too. The test pins a property of the new order; it would not fail if `b3ad41a` were
reverted. Worth knowing when you read the suite as a regression net — of the ten new tests, the
revert detectors are `test_poly_a_tail_behind_poly_g_is_trimmed`,
`test_poly_g_behind_poly_a_survives` (in the opposite direction) and
`poly_a_tail_behind_poly_g_is_trimmed`; the other seven are controls, each deliberately so.

**The comment at `:189-193` is stronger than the mechanism:** *"So however much poly-G takes off the 3'
end, the 5' clip is the same."* The 5' scan iterates over the sequence it is given, so truncating the
3' end removes candidate indices; `best_index` over the shorter prefix set is ≤ `best_index` over the
full one, equal only when the head does not extend into the removed region. "Never competes, because
the head and the run are disjoint" is the honest form, and it fits in two lines. (The plan asked for
"`--paired --poly_a` R2 byte-identity as a committed test, not a manual before/after" — what shipped
compares `--poly_g` against `--no_poly_g` inside one build, which is a weaker property. I measured the
stronger one by hand at V5: 0 of 13,000 R2 records differ from the pre-change binary. Pinning that
would need a golden file, so I would not ask for it; just note the plan item is partially met.)

### F8 — LOW — magic numbers, and an invariant recorded in the tests instead of on the constant

`src/trimmer.rs:1181` — `INSERT_60` ends in `G`. Three assertions depend on that:
`[..59]` and `poly_g_trimmed == 16` at `:1270`/`:1272` (15 run + 1 insert base), and
`test_stage_order_does_not_touch_reads_without_tails` passes only because a lone trailing G is below
the 3-bp floor in `homopolymer_trim_index`. Change that one base and three tests fail with messages
that point at the wrong thing. CLAUDE.md: *"Put an invariant on the thing it constrains, not in the
function that relies on it."* Note it on the constant, not in the test that trips over it:

```rust
    /// 60 bp of insert, used by the stage-order tests below. Its 3' base is a G, which the R2 test
    /// counts as part of the poly-G run.
    const INSERT_60: &str = "…";
```

Also note that `tests/integration_poly.rs:27`'s `INSERT` is the **same 60 bp except the last base**
(`…AATCT` vs `…AATCG`) — a difference a reader will not spot, and a load-bearing one: the integration
assertions are exact equality against `INSERT`, which requires a terminal base that is neither A nor G.
Its doc comment ("no poly-A or poly-G content of its own") gestures at this without saying it.

The literals themselves would be better derived:

| Site | Now | Suggested |
|---|---|---|
| `src/trimmer.rs:1229` | `assert_eq!(record.seq.len(), 75);` | `assert_eq!(record.seq.len(), INSERT_60.len() + 15);` |
| `src/trimmer.rs:1245` | `assert_eq!(record.seq.len(), 100, "nothing is trimmed");` | `assert_eq!(record.seq.len(), read_seq.len(), "nothing is trimmed");` |
| `src/trimmer.rs:1270` | `assert_eq!(record.seq, INSERT_60[..59]);` | `assert_eq!(record.seq, INSERT_60[..INSERT_60.len() - 1]);` |
| `src/trimmer.rs:1272` | `assert_eq!(result.poly_g_trimmed, 16);` | `assert_eq!(result.poly_g_trimmed, 15 + 1, "the run plus the insert's own trailing G");` |

The `:1245` one is a strict improvement rather than a restatement — `read_seq.len()` cannot drift out of
step with the read it was built from.

### F9 — LOW — `test_poly_g_tail` is now a duplicate of `test_poly_g_clear_tail`

`src/quality.rs:517` is `ACGTACGTGGGGG` → 8; `src/quality.rs:433` is `ACGTACGTGGGGGGGGG` → 8. Same
function, same expectation, same boundary, differing only in run length. The rename the plan asked for
is correct — the old name encoded the old order — but after it the test has no distinct purpose; its
coverage is strictly inside `test_poly_g_clear_tail`'s. I would delete `test_poly_g_tail`. Flagging
rather than acting because the plan explicitly specified the rename, so the call is yours.

Related placement nit: `test_poly_a_short_reads` (`:501`) is a **poly-A** test sitting below the
`// --- Poly-G / Poly-C trimming tests` banner at `:430`, between two poly-C tests. It belongs with the
poly-A block at `:349-428`. Its comment ("length filtering runs after trimming") states something about
`trim_read` that this test does not exercise — it tests `poly_a_trim_index` indices only. I traced all
four assertions for `n = 3..=10` and they are correct and non-vacuous (they pin the `n < 3` boundary and
the no-tail direction across the whole window), so this is placement and wording only.

### F10 — LOW — `&PathBuf` parameters

`tests/integration_poly.rs:32`, `:41`, `:177` take `&PathBuf`. Confirmed not a CI problem (V2:
`clippy::ptr_arg` skips arguments forwarded to generic `AsRef<Path>` sinks, and a forced full
`--all-targets -D warnings` run is clean), but they are the only `&PathBuf` function parameters in
`tests/` or `src/`, and the sibling file this one is modelled on uses `&Path`
(`tests/integration_partial_output.rs:37`). `-> fn write_fastq(path: &Path, …)` and the same for
`sequences`/`records`; no call site changes, deref coercion handles it.

### F11 — LOW — the `--help` text is unconditional about poly-G

`src/cli.rs:365-367`: *"Runs after adapter and poly-G trimming, so tails behind either are also
removed."* True when poly-G is active, which on 2-colour data it usually is by auto-detection — but
under `--no_poly_g` the second clause is exactly what does **not** happen, and that is the case the
reporter will hit. `docs/modes/poly.md` covers it; `--help` reads as a guarantee.
`…so tails behind either are also removed (poly-G trimming permitting)` or similar closes it in one
clause.

### F12 — LOW — the canonical pipeline-order list still omits the poly stages

`docs/src/content/docs/guide/overview.md:16-24` enumerates "The pipeline runs in this order" as
1 quality / 2 adapter / 3 RRBS / 4 length filter / 5 PE validation / 6 demux+FastQC — no poly-G,
no poly-A, no N-trim, no fixed clipping. Pre-existing, but this commit's stated job is "record the
stage order", and `trim_read`'s own doc comment was updated for exactly that reason (`src/trimmer.rs:85-86`,
a good change). A reader who wants the order will look at the overview page first and find the poly
stages absent. One inserted step between 3 and 4, linking `/modes/poly/#stage-order`, finishes the
commit's own goal.

### F13 — LOW — fixed temp-directory names

`fresh_tmpdir` uses a fixed `std::env::temp_dir().join(slug)`. Two concurrent `cargo test` runs, or two
CI jobs on one runner, clobber each other's `/tmp/tg_poly_*`, and the directories survive the run.
Half the suite guards against this with `std::process::id()` (`integration_adapter2.rs:33`,
`integration_library.rs:19`, `integration_output_collision.rs:32`, five others) and half does not
(`integration_ubam.rs:23`, `integration_passthrough.rs:27`, `integration_partial_output.rs:24`).
The new file copies the second group — including the helper, byte for byte, from
`integration_partial_output.rs:24-30`.

Worth noting for the record: `PLAN.md` commit-2 item 3 asked for `tempfile::TempDir` and cited
`tests/integration_partial_output.rs:525` as the pattern. That citation is wrong — that file uses the
hand-rolled helper, and `tempfile::` appears nowhere under `tests/` despite being a dev-dependency.
**The implementation followed the real convention over the plan's misdescription of it, which was the
right call.** If you want the `process::id()` suffix, that is the cheaper half of the fix:
`.join(format!("{slug}_{}", std::process::id()))`.

---

## Verified clean — things I checked and found no problem with

- **Stage order across every entry point.** All eight production `trim_read` call sites —
  `run_single_end` (`:324`), `run_paired_end` (`:484-485`), `run_single_end_to_bam` (`:650`),
  `run_paired_end_to_bam` (`:735-736`), `parallel.rs:495-496`, `parallel.rs:1133` — funnel through the
  one function, so SE, PE, `--cores N` and both uBAM-output paths inherit the fix uniformly. Specialty
  modes bypass trimming entirely. `--clump_only` rejects both poly flags at `src/cli.rs:863-867`.
- **Stat attribution.** `poly_a_trimmed` / `poly_g_trimmed` are each computed inside their own stage
  from that stage's own `seq_len_before`, so swapping the blocks cannot cross-attribute. The
  `> 0 → count a read` aggregation in `parallel.rs:517-531` / `:1142-1148` and
  `report.rs:97-100` is order-independent.
- **`bp_after_cutadapt`.** Captured at `src/trimmer.rs:160`, before the RRBS block and therefore before
  both poly stages — unchanged by the reorder, so the `#232` "Total written (filtered)" parity for RRBS
  is safe.
- **Empty-record handling.** Both stages guard with `!record.is_empty()`, and the poly-A guard is now
  downstream of a stage that can empty the record (an all-G read truncates to 0). Checked: skipping
  poly-A on an emptied record gives the same output as the old order, where poly-A ran first and found
  no A's.
- **`--rrbs`.** RRBS is step 2.5, upstream of both poly stages in both orders; relative order unchanged.
- **R2 semantics.** Poly-T off the 5' end, poly-G off the 3'; opposite ends, no interaction. Measured:
  0 of 13,000 R2 records change vs the pre-change binary (V5).
- **Perl byte-identity invariant.** No CI validation case uses `--poly_a` (V10).
- **The `--2colour` test isolates the quality path** (V6) — my first read of it was wrong; the
  auto-detection floor of 10 reads keeps poly-G off.
- **Doc links and the `#stage-order` anchor** (V9).
- **`docs/.../reference/changelog.md` not updated.** Not a finding: that synced copy is already many
  entries behind (`#457` and everything since are absent), so the convention is clearly to sync at
  release time.
- **No stale copies of the old `--poly_a` help string** anywhere in `src/`, `tests/` or `docs/src/`.
- **`b"CGTCGTCGTC"[..n].to_vec()`** in `test_poly_a_short_reads` — compiles and is correct; I hand-traced
  all four assertions for every `n` in `3..=10`.

## Out of scope, noted in passing

`src/cli.rs:864` and `:867` word the `--clump_only` rejections as "`--polyA` is not compatible" /
"`--polyG` is not compatible". Those are the Perl-era aliases (still accepted via `alias = "polyA"`), so
the messages are not wrong — but the whole point of the migration section this PR rewrites is that
`--polyA` and `--poly_a` are different things. Pre-existing and a separate change.

---

## Priority order

1. **F1** (HIGH) — fix or drop the `AAAAAAAAAA` bullet before the docs commit lands; it is wrong in a
   way that also invalidates the plan's `A{10}` baseline figure.
2. **F2, F3, F4, F5** (MEDIUM) — the `--length` claim, the whitespace omission, the vacuous-pass hole,
   and the 3-line comment. All one-line or one-sentence fixes.
3. **F6–F13** (LOW) — credit, comment precision, magic numbers, duplicate test, `&PathBuf`, `--help`
   caveat, overview page, temp dirs.

`b3ad41a` (the reorder) needs no change. `23cb2a6` needs F4 and the comment trims. `869912c` should not
ship as written.
