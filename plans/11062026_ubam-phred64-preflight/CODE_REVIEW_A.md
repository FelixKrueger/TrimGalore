# Code Review A — `--phred64` pre-flight widened to uBAM output

**Reviewer:** A (independent, fresh context)
**Target:** uncommitted changes on `dev`, base `49152c1`
**Files:** `src/main.rs`, `tests/integration_fastqc_phred64.rs`, `CHANGELOG.md`, `docs/src/content/docs/guide/{quality,outputs}.md`
**Mode:** recommend-only (fix clause suspended by the caller). No `src/`, `tests/`, `docs/` or `CHANGELOG.md` file was modified.
**Build:** nothing was built or mutated. `cargo fmt --all -- --check` was run (read-only, exit 0). All behavioural checks used the existing `target/release/trim_galore`, whose mtime (2026-11-08T10:09:50) post-dates `src/main.rs` (10:04:15), so it carries the change under review. Clippy was **not** run (it would write to `target/` and race the concurrent agents).

---

## Summary

**The change is correct.** The condition, its siting, and its coverage of every dispatch branch all hold, and I verified each empirically rather than by reading alone. No currently-working invocation is broken: genuine Phred+64 → uBAM still works in all four modes that exercise it, and I confirmed *why* it still works (see F-3 — it rests on a fixture property that nothing records).

**No Critical and no High findings.** The defects I found are a user-facing overclaim in the startup NOTE (F-1), a docs mis-categorisation (F-2), and an undocumented test-fixture coupling the change newly introduces (F-3). The rest are Low.

### What I verified by running the binary

| Scenario | Expected | Observed |
|---|---|---|
| Phred+33 → uBAM, no `--fastqc` | refuse, no output dir | exit 1, pre-flight message, `out/` absent |
| Phred+33 → `--hardtrim5 20` + uBAM | refuse, no output dir | exit 1, pre-flight message, `out/` absent |
| Phred+33 → `--clump_only` + uBAM | refuse, no output dir | exit 1, pre-flight message, `out/` absent |
| Phred+33 → FASTQ, no `--fastqc` | unchanged, succeeds | exit 0, `_trimmed.fq` written |
| Genuine Phred+64 → uBAM, SE trim | succeeds, true Phred | exit 0, `QUAL` = Q40 run then Q2 run |
| Genuine Phred+64 → uBAM, `--hardtrim5` | succeeds | exit 0, `.20bp_5prime.bam` |
| Genuine Phred+64 → uBAM, `--clump_only` | succeeds | exit 0, `_clumped.bam` |
| Genuine Phred+64 → uBAM, paired | succeeds | exit 0, `_R1_val.bam` |
| **gzipped** Phred+33 → uBAM | refuse | exit 1 — `open_sync_reader` is format-aware, so `.gz` is scanned |
| Paired, genuine R1 + mislabelled R2 | refuse, naming R2 | exit 1, message names `mix_R2.fastq` |
| Phred+33 with every byte ≥ ASCII 64 | accepted (blind spot) | exit 0; Q40 `I` stored as Phred 9 (`*`) |

The last row is the limit `PLAN.md` §13 records, reproduced: a Q40-throughout Phred+33 file is silently written 31 points low. See F-1 and F-2 for whether the prose describes that honestly.

### The condition, siting and coverage

```rust
if cli.phred64
    && (cli.fastqc_requested()
        || matches!(cli.output_format, trim_galore::cli::OutputFormat::UBam))
```

- **Boolean is right.** `fastqc_requested()` is `fastqc || fastqc_args.is_some()`, so an args-only run is covered; the uBAM disjunct is the only output format that consumes `cli.phred_offset()` on the write side. The fully-qualified `trim_galore::cli::OutputFormat::UBam` matches the other nine inline matches in `main.rs`, so it is house style, not a one-off.
- **Siting is right.** Line ~386, which is ahead of `naming::ensure_output_dir` (~552), ahead of all three output-collision pre-flights, and ahead of every `return`ing dispatch branch. The `!out_dir.exists()` assertions are not merely passing — I confirmed the directory is genuinely never created in all three refusal scenarios.
- **Early-returning branches are covered.** `--hardtrim5/3` and `--clump_only` all return well past line 386; verified by running each. `--clock`/`--implicon` need no coverage because `Cli::validate()` §3.4a rejects them with `--output-format ubam` outright.
- **uBAM input cannot reach the loop.** `cli.phred64 && any_bam` bails at ~369, so within this loop no `cli.input` entry can be `UnalignedBam`. The `continue` guard is therefore live only for the `--passthrough` carrier, which is what its comment claims. I confirmed `--passthrough` + `--output-format ubam` is itself rejected at startup, so on the uBAM disjunct the passthrough chain is inert — the comment stays accurate.

### The four new tests

All three refusal tests assert `assert_preflight_refused`, which greps for `"--phred64 declares Phred+64 input"` — a string unique to this `bail!` and distinct from the BAM-input refusal (`"--phred64 cannot be used with unaligned BAM input"`). **They therefore cannot pass on an unrelated non-zero exit**, including the `--clump_only` + `--length` trap the caller flagged. I looked for others of that shape and found none: `--hardtrim5 20 --output-format ubam` and `--clump_only --output-format ubam` are both otherwise-valid invocations (confirmed by running them on genuine Phred+64 input, where they succeed).

`run_raw`'s first argument is appended as `--output_dir`, not used as a cwd, so passing a not-yet-existent `out_dir` is correct and makes the `!out_dir.exists()` assertion load-bearing.

---

## Findings

### F-1 — Medium — the startup NOTE now overclaims, and the new test pins the vague half

`src/main.rs:567`, closing sentence of the NOTE:

> A mislabelled input is refused before anything is written when the output is uBAM, or when `--fastqc` is given.

This states an unconditional guarantee. It is not one: the check refuses on *any quality character below ASCII 64 within the first 10,000 reads*, and I reproduced the case it cannot settle (table row 11 — a Phred+33 file with every byte ≥ 64 is accepted and written 31 points low). `docs/guide/quality.md` is scrupulous about this in the same change; the NOTE is not.

The asymmetry matters because of who reads what. The NOTE fires on stderr at the moment the user is taking the risk; `quality.md` is read later, if at all. The sentence it replaced was precise about the condition —

> With `--fastqc` the FastQC step refuses a file holding any quality character below ASCII 64.

— and the rewrite trades that precision away. (The rewrite does fix a separate staleness in that old sentence: since #532 it is the pre-flight, not the FastQC step, that refuses. So the edit is a net gain in attribution and a net loss in precision; only the second half needs fixing.)

Compounding it, `note_qualifies_pass_through_by_output_format` asserts on `"writing uBAM they do not"` — the pass-through clause — and asserts nothing about the refusal sentence, so the overclaim is unpinned while the prose around it is pinned hard.

**Recommend** restoring the condition in the same breath, which costs no length:

> An input holding any quality character below ASCII 64 is refused before anything is written when the output is uBAM, or when `--fastqc` is given.

### F-2 — Medium — `outputs.md` files `--phred64` as a v1 scope gap, which it is not

`docs/src/content/docs/guide/outputs.md:101`. The new bullet joins a list headed:

> Most FASTQ-shaped features work with uBAM output. Rejected at startup, before anything is written (v1 scope):

Every other entry is a feature that genuinely does not work with uBAM output — `--demux`, `--clock`, `--retain_unpaired`, `--clumpify`. `--phred64` is the opposite: it is **supported** with uBAM output, deliberately so since #362, and what the new code rejects is *invalid input*, not the flag. Input validation is being catalogued as a scope limitation.

The `--rename` entry above it is not a precedent — `--rename` really is refused for a broad, legitimate class of runs. A reader skimming the bullet list (which is how a compatibility table gets read) will take away "`--phred64` is out of scope for uBAM in v1", which is wrong and contradicts `quality.md`.

It is also three sentences and ~60 words in a list of one-clause entries, so it does not read as a list item.

**Recommend** either lifting it out of the list into the Phred cross-reference, or at minimum leading with the supported status so the skim reads correctly:

> - `--phred64` — supported; only the run whose input contradicts the declaration is refused, since the declared offset is subtracted on write (see [Phred encoding](/guide/quality/#phred-encoding))

### F-3 — Medium — the change creates an undocumented fixture coupling across three test files

Four pre-existing guards now pass through the widened gate on every run:

- `tests/integration_ubam_out.rs` — `phred64_ubam_out_se_stores_true_phred`, `phred64_ubam_out_specialty_stores_true_phred`, `phred64_ubam_out_pe_two_file_stores_true_phred`
- `tests/integration_clump_only_ubam.rs` — `phred64_clump_only_ubam_out_stores_true_phred`

None passes `--fastqc`, so before this change none of them touched the pre-flight at all. They survive it for one reason, which I measured: the minimum quality byte in `test_files/phred64_test.fastq`, `phred64_test_R1.fastq` and `phred64_test_R2.fastq` is **66** (`B`, Q2 under Phred+64) — comfortably above the ASCII-64 floor.

That is now a load-bearing invariant on three fixtures, and nothing records it. The next person who extends `phred64_test.fastq` with a lower tail — a natural thing to do, since `-q` behaviour is exactly what that fixture exists to exercise — gets four failures in two *other* files, reporting a `--phred64` encoding refusal rather than anything resembling what they changed.

`CLAUDE.md` is explicit about where this belongs: *"Put an invariant on the thing it constrains, not in the function that relies on it."* The fixtures are plain FASTQ with no comment surface, so the nearest honest home is the fixture inventory.

**Recommend** one line in `CLAUDE.md`'s "Test fixtures" paragraph, or failing that in the new test block:

```
// test_files/phred64_test*.fastq must hold no quality byte below ASCII 64:
// the pre-flight refuses any --phred64 run whose input does.
```

### F-4 — Low — a refusal test carries flags its sibling's comment argues are pointless

`preflight_fires_on_ubam_output_without_fastqc` passes `-q 0 --length 0`. `preflight_fires_on_clump_only_ubam_output`, eight lines later, explains at length that it omits both because *"neither flag is needed, because the refusal precedes every dispatch branch"* — which is true, and equally true of the first test.

Harmless (the message assertion protects it either way), but the two tests now make opposite statements about the same fact, and the inherited `-q 0 --length 0` is a copy from `preflight_refuses_mislabelled_input_before_writing_anything`, where it was needed for a different reason. **Recommend** dropping both flags from the uBAM test.

### F-5 — Low — the NOTE test pins a five-word fragment mid-sentence

`note_qualifies_pass_through_by_output_format` asserts `stderr.contains("writing uBAM they do not")`. I confirmed the string survives the `\`-continuations and matches at runtime, so the test works today. But it is a mid-sentence fragment with no anchor, it pins the clause that needed no pinning while leaving the clause that does (F-1) unasserted, and any rewording of that sentence — including the F-1 fix, which is adjacent — breaks it for no behavioural reason.

**Recommend** asserting the two facts the test is named for separately and more loosely, e.g. that stderr mentions `--hardtrim5/3` pass-through *and* that the run succeeded, rather than one verbatim phrase.

### F-6 — Low — `integration_fastqc_phred64.rs`'s module doc is now stale

The file's `//!` header opens *"`--phred64` must reach the bundled FastQC config"* and explains the FastQC 0.13.0 encoding problem. The file now holds a block of tests about uBAM `QUAL` arithmetic in which FastQC plays no part, and the pre-flight they exercise is no longer FastQC-scoped at all — the `preflight_does_not_fire_without_fastqc` → `..._on_fastq_output` rename in this very diff is the acknowledgement of that.

Per `CLAUDE.md` (*"Current state only"*), the header should now describe the file as it is. **Recommend** either a sentence added to the header covering the pre-flight, or moving the four new tests to `integration_ubam_out.rs` where the other uBAM `--phred64` guards live. I lean to the header fix — the pre-flight tests are already here and splitting them would be worse.

### F-7 — Low — the source comment's second sentence is detail, not the fact

`src/main.rs:386-387`:

```rust
// On uBAM output the declared offset feeds the QUAL arithmetic, and bases
// below Q31 floor at 0 rather than shifting.
```

The first clause is the fact that justifies the condition and should stay. The second — *"bases below Q31 floor at 0 rather than shifting"* — is a consequence already stated in the startup NOTE, in `quality.md` and in the CHANGELOG entry, and *"rather than shifting"* is the shape `CLAUDE.md` calls out: a contrast that pre-empts an objection instead of stating current state. It also takes the comment block (lines 381–387) from five lines to seven, against the file's own one-to-two-line default.

(The arithmetic is correct, for what it is worth: a true-Q*q* base is byte *q*+33, read as *q*−31, so *q* < 31 clamps to 0.)

**Recommend** `// uBAM output: the declared offset feeds the QUAL arithmetic.` — one line, same justification, no duplication.

### F-8 — Low — redundant format detection inside the loop

```rust
for input in cli.input.iter().chain(cli.passthrough.iter()) {
    if matches!(detect_input_format(input)?, InputFormat::UnalignedBam) {
```

`input_formats` already holds the detected format of every `cli.input` entry, computed ~20 lines above precisely so it can be reused ("*Used for both validation and dispatch below*"). This re-detects each one — and, as established above, can never return `UnalignedBam` for them, because `phred64 && any_bam` has already bailed. Only the passthrough entry needs a live call.

Pre-existing, but the widened gate makes it fire on a new class of runs (`--phred64` + uBAM output, which previously paid nothing). The cost is a 24-byte peek plus at most one BGZF block per input, so this is tidiness, not performance. **Recommend**, if touched at all, `cli.input.iter().zip(&input_formats)` with the passthrough handled separately.

### F-9 — Low — two CHANGELOG entries now split one condition across one issue number

`CHANGELOG.md` carries, in the same release section, two entries both citing #532:

- *"`--phred64` is checked against the input before anything is written, when `--fastqc` is in play"*
- *"`--phred64` is also checked against the input whenever the output is uBAM"*

They describe the two disjuncts of a single `if`. A reader deciding whether the upgrade affects them has to find both and reconcile *"on FASTQ output, `--phred64` without `--fastqc` is unaffected"* in the first against the second's scope. Since the first entry is also unreleased (no `v2.4.0` tag exists — only `v2.3.0`; `Cargo.toml` is at 2.4.0 and two commits have landed under the dated heading since `59993c0 release(v2.4.0)`), **both are still editable**, and merging them into one entry stating the final condition would serve that reader better than two halves.

The new entry on its own is well-formed and close to `CLAUDE.md`'s one-sentence rule — closer than its neighbours, which run to full paragraphs. The house style in this file is visibly looser than the rule, so I am not flagging its length; only the split.

*(Noted, not flagged: the heading reads "Release on 3 November 2026" while content is still accruing on 8 November with no tag cut. Pre-existing, nothing to do with this change.)*

---

## Not findings — checked and clean

- **`quality.md` is honest about the blind spot.** The new sentences state it explicitly, name the shape of file that reaches it (*"synthetic and heavily filtered"*), and are accurate about its rarity — Q2 tails and binned low scores do both fall below the threshold. The CHANGELOG entry's *"see `docs/guide/quality.md` for what the check can and cannot settle"* points at it rather than papering over it. This is the accuracy question the caller raised, and apart from F-1 and F-2 the prose answers it correctly.
- **The reworded pass-through clause is right.** *"Writing FASTQ, the modes that do no quality arithmetic — `--hardtrim5/3` and `--clump_only` — pass the scores through untouched; writing uBAM they do not"* matches the code: `specialty::hardtrim{5,3}_to_bam` and the `--clump_only` BAM writers all receive `cli.phred_offset()`, while their FASTQ counterparts do not. The ambiguity the caller asked about is resolved.
- **Solexa input (ASCII 59–63) is now refused on the uBAM path** where it previously wrote floored scores. That is a behaviour change, but the right one, and the `bail!` already names Solexa explicitly.
- **The unguarded case is the loud one.** `-q 20 --phred64` on Phred+33 → FASTQ is still accepted, but I ran it: every read is trimmed to nothing and the report says *"Sequences removed because they became shorter than the length cutoff of 20 bp: 4 (100.0%)"* against an empty `.fq`. That announces itself; the uBAM case does not. Scoping the gate to uBAM and `--fastqc` is therefore defensible and I am not asking for it to be widened. Worth knowing that the NOTE's opening clause names *"quality trimming, and BAM QUAL on uBAM output"* together, which may lead a reader to expect symmetric protection — one more reason to make the closing sentence state its actual condition (F-1).
- `cargo fmt --all -- --check` exits 0.

---

## Priority order

1. **F-1** (Medium) — NOTE overclaims a guarantee the check does not provide; one-sentence fix.
2. **F-2** (Medium) — `outputs.md` presents supported behaviour as a v1 scope gap.
3. **F-3** (Medium) — record the ASCII-64 fixture invariant before it bites someone in two unrelated test files.
4. **F-4**–**F-9** (Low) — test tidiness, stale module doc, comment length, redundant detection, CHANGELOG split.

Nothing here blocks the change on correctness.
