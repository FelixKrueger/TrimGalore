# Code review — PR #536 `--phred64` pre-flight (Reviewer B, gen1)

**Target:** `ebfb84f`, `e4caf3c` on `fix/phred64-preflight-532`, against base `origin/dev` = `a4b67ee`
**Date:** 2026-11-06
**Method:** behaviour-first. Two pinned `git archive` copies (`e4caf3c` → `$TMPDIR/revB`, `a4b67ee` →
`$TMPDIR/revB-base`), each built `--release` with its own `CARGO_TARGET_DIR`. 38 invocations run
against **both** binaries and compared on exit code, the exact set of files left on disk, and md5 of
every data output. The shared working tree was never touched.

**Verdict: the three claims hold. Ship-able. Two Medium findings about the refusal's stated reason
and its scoping, three Low findings about message/test strength.**

---

## 1. The three claims, checked

### 1.1 "Refusal writes nothing" — CONFIRMED, including the output directory

14 refusal paths, every one left **zero** files. The comparison against base is the point:

| invocation (`--phred64 --fastqc` throughout) | base `a4b67ee` | new `e4caf3c` |
|---|---|---|
| SE, plain Phred+33 | rc=1 after `*_trimmed.fq` + .txt + .json | rc=1, **nothing** |
| SE, gzipped Phred+33 | rc=1 after `*_trimmed.fq.gz` + 2 reports | rc=1, **nothing** |
| PE, both mates Phred+33 | rc=1 after 2 mates + 4 reports | rc=1, **nothing** |
| PE, R1 good / R2 Phred+33 | rc=1 after 5 files **incl. R1's FastQC .html/.zip** | rc=1, **nothing** |
| multi-pair, pair 2 Phred+33 | rc=1 after **17 files**, pair 1 fully QC'd | rc=1, **nothing** |
| `--passthrough`, lockstep Phred+33 carrier | rc=1 after **11 files** incl. both mates' zips | rc=1, **nothing** |
| `--clump_only` (FASTQ out) | rc=1 after `*_clumped.fq` + clumping report | rc=1, **nothing** |
| boundary bytes 58 / 59 / 63 | rc=1 after output + 2 reports each | rc=1, **nothing** |
| `--fastqc_args --nogroup` (no bare `--fastqc`) | rc=1 after 3 files | rc=1, **nothing** |

The multi-pair and passthrough rows are the strongest cases for the change: base completed pair 1
end-to-end, FastQC reports and all, before dying on pair 2.

**`--output_dir` is not created either.** Pointed at a path that did not exist: base created it and
wrote three files into it; new left the parent empty. This is correct by construction — the guard sits
at `src/main.rs:384`, `naming::ensure_output_dir` at `src/main.rs:533`. Also verified that nothing
lands in the CWD (each probe ran with cwd = its own empty dir).

### 1.2 "Genuine Phred+64 still works" — CONFIRMED

SE, PE, `--passthrough`, `--clump_only` (FASTQ and uBAM out), `--output-format ubam`, boundary bytes
64 / 65 / 66, and a CRLF-terminated Phred+64 file all exit 0 with a full report set. The CRLF case
matters: `FastqReader` strips `\r` at `src/fastq.rs:433`, so a DOS-line-ended Phred+64 file is **not**
falsely refused on ASCII 13. Verified empirically, not just read.

The 10,000-read window is **exact**, no off-by-one: a sub-64 byte in read 10,000 is caught by the
pre-flight; the same byte in read 10,001 falls through to FastQC. Both directions tested with
purpose-built 10k/10k+1 fixtures.

### 1.3 "`--phred64` without `--fastqc` is byte-for-byte unchanged" — CONFIRMED

md5-compared `*.fq`, `*.fq.gz`, `*.bam` and `*_trimming_report.txt` for every rc=0 case across both
binaries. Identical throughout. The only three diffs anywhere are the embedded invocation path —
`"command_line"` in the JSON report, `@PG CL:` in the BAM, and the `Output:` line of the passthrough
report — all of which differ solely because the two builds live in different directories. No
behavioural delta.

---

## 2. Findings

### M1 (Medium) — the refusal's closing sentence is false on the uBAM-output paths

`src/main.rs:397`:

```
... those scores floor to Q0. FastQC rejects such a file either way.
```

It does not. The uBAM-output call sites pass `phred64: false` to `fastqc::run` —
`src/main.rs:2428`, `src/main.rs:2588`, `src/main.rs:2708`, and `src/clump_only.rs:917`,
`src/clump_only.rs:1177` — because the BAM declares its own encoding and the format wins upstream
(which is exactly what the pre-existing test `bam_format_encoding_wins_over_the_phred64_flag`,
`tests/integration_fastqc_phred64.rs:96`, pins). Verified on base:

```
ubamout_p33         --phred64 --fastqc --output-format ubam                 rc=0
clumponly_ubam_p33  --phred64 --fastqc --clump_only --output-format ubam    rc=0
```

Both produced a FastQC report that **passes**: Encoding `Sanger / Illumina 1.9`, per-base mean
`0.0` across all 34 bases. The new guard refuses both.

Refusing is the right outcome — I decoded the base's BAM and every `QUAL` byte is 0 — but the message
tells the user something checkable and wrong, and the CHANGELOG's framing ("instead of being trimmed,
reported, and only then rejected by FastQC") doesn't hold on those paths either. Suggest dropping that
sentence or qualifying it ("FastQC rejects such a file on the FASTQ output path").

### M2 (Medium) — scoping to `--fastqc` leaves the worst case uncovered

The same mislabelled input, same uBAM output, **without** `--fastqc`, still exits 0 on `e4caf3c`.
Decoded from `p33_trimmed.bam`:

```
read 0: l_seq=34 qual bytes=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
```

versus a genuine Phred+64 input through the same path: `[40, 40, 40, ...]`. Every quality score is
destroyed, silently, rc=0. This is precisely the failure the pre-flight exists to prevent, and the
`--fastqc` condition is what lets it through.

The `--fastqc` scoping is defensible on the FASTQ-output path (there, `--phred64` is only a claim about
how to *read* bytes that get copied through verbatim). It is not defensible on the uBAM path, where
`cli.phred_offset()` feeds real arithmetic into `BamWriter`. The pre-flight predicate is already
written and costs ~10k records; running it unconditionally whenever
`matches!(cli.output_format, OutputFormat::UBam)` would close the hole. Recommendation, not a blocker
for this PR — but note the PR currently refuses the *less* damaging variant of this pair and permits
the *more* damaging one.

### M3 (Medium) — the guard is strictly stronger than FastQC's own refusal, which the CHANGELOG implies it isn't

The pre-flight scans the **input**; FastQC sees the **output**. Trimming is only ever substring +
whole-read-drop, so the output's quality bytes are a subsequence of the input's and `min(output) >=
min(input)`. There are therefore inputs the pre-flight refuses that FastQC would have accepted.
Verified with a fixture of 30×`'h'` + 4×`'#'`:

```
tail_low  -q 20 --length 10   base: rc=0, full FastQC .html/.zip written
                              new:  rc=1, refused
```

Base trimmed the low tail away, FastQC saw only bytes ≥ 64 and passed. So this change does remove
previously-working invocations. I think that is the right call — the input genuinely is mislabelled —
but the CHANGELOG at `CHANGELOG.md:21-26` presents the pre-flight purely as an earlier version of
FastQC's verdict, and it is a broader one. One clause would cover it.

### L1 (Low) — `lowest as char` emits raw control bytes

`src/main.rs:399`. A quality line holding ASCII 11 produces, literally:

```
quality character '<VT>' (ASCII 11) within the first 10000 reads
```

with the raw byte in the terminal (hexdump-verified). ASCII 32 renders as an empty-looking `' '`.
Neither is valid FASTQ, but mangled quality lines are a real reason someone reaches for `--phred64`.
`lowest.escape_ascii()` renders `\x0b` / `\x20` instead and costs nothing.

### L2 (Low) — the Solexa branch fires for bytes the sentence itself excludes

`src/main.rs:395` offers "or it is Solexa / Illumina 1.0, which allows ASCII 59-63" for **every**
byte < 64. Verified: ASCII 58 and ASCII 35 both get that sentence, and both are below the Solexa floor
it names. The range given is correct (Solexa −5 is ASCII 59), so a careful reader can rule it out, but
the message presents an either/or where only one branch fits — and the common trigger (`'#'`=35,
`'!'`=33) is always the Phred+33 branch. Gating the Solexa clause on `(59..=63).contains(&lowest)`
makes every instance of the message true of its own input. Minor; `e4caf3c` was right to add Solexa,
this is just about when to say it.

### L3 (Low) — the "nothing written" test cannot catch a regression in the guard's placement

`tests/integration_fastqc_phred64.rs`, `preflight_refuses_mislabelled_input_before_writing_anything`
passes `--output_dir` pointing at an **already-existing** `tempdir()`, so `ensure_output_dir` is a
no-op there. The test would still pass if the guard were moved below `src/main.rs:533`. Pointing
`--output_dir` at a non-existent subdirectory and asserting `!subdir.exists()` afterwards pins the
invariant the comment at `src/main.rs:366` actually claims. (I verified the current code gets this
right; the test just doesn't hold it there.)

### L4 (Low) — no test covers paired or multi-pair, which is where the change pays most

The four new tests are all single-end. The largest observed improvement is multi-pair: 17 files on
base (pair 1 fully trimmed and QC'd) versus 0 now. Nothing locks that in. A paired test asserting an
empty output directory would be ~15 lines and would catch any future reordering that moved the guard
below the per-pair loop.

### L5 (Low) — `lowest_quality_char` is unreachable from unit tests

It lives in the binary crate (`src/main.rs:34`), so `cargo test --lib` cannot reach it and the
10,000-read boundary is only exercisable through a ~790 KB integration fixture (which is why no test
covers it). Moving it to `src/format.rs` next to `open_sync_reader` would make the boundary, the
`None`-on-no-records case, and the empty-qual case cheap to pin.

### Nits

- `src/main.rs:26-48` — the const and helper sit between the `use` block and the
  `#[global_allocator]` declaration, which conventionally leads the file.
- `bail!` inside the loop names only the first offending file; a user with both mates mislabelled
  fixes R1, re-runs, and learns about R2. Standard fail-fast, mentioned only for completeness.
- The inline comment at `src/main.rs:379-383` is five lines and carries the change's rationale, which
  `~/.claude/CLAUDE.md` puts in the commit message. The surrounding `#358` block is 18 lines in the
  same register, so this matches local convention rather than breaking it — noted, not asked for.
- The CHANGELOG entry is three sentences against the global "one sentence per entry" rule, but every
  neighbouring entry in the file is a prose paragraph too. Changing this one alone would read worse.

---

## 3. Docs and CHANGELOG accuracy

`docs/src/content/docs/guide/quality.md:40`:

- "prints a `NOTE:` on every output path when `--phred64` is given" — **true**, verified: exactly one
  NOTE on `--hardtrim5`, `--hardtrim3`, `--clump_only` and normal trim.
- "the claim is also checked against the input before any work starts" — **true**.
- "The scan covers the first 10,000 reads, so a file whose only low character sits past that window is
  still caught later, inside FastQC" — **true and exact**.
- **Gap:** that later catch still leaves `*_trimmed.fq` and both reports on disk (verified:
  `late_low` fixture, 10,001st read low → rc=1 with three files written). Since the whole entry is
  about not leaving artefacts behind, the sentence that admits the window should say the fallback
  path still part-processes.
- "a file that does is refused up front rather than part-processed and then rejected by FastQC" —
  carries the same M1 caveat: on `--output-format ubam` FastQC never rejects it.

`CHANGELOG.md:21-26`:

- "The scan covers the first 10,000 reads; `--phred64` without `--fastqc` is unaffected" — both
  **verified**, the second byte-for-byte.
- "instead of being trimmed, reported, and only then rejected by FastQC" — **true on the FASTQ-output
  paths, false on uBAM output** (M1), and narrower than what the guard actually refuses (M3).

---

## 4. Bypass hunt (brief item 2)

Every early-returning path in `main()` is covered or refused earlier — the guard at `src/main.rs:384`
precedes all of them (specialty dispatch starts at `:584`, `--clump_only` at `:700`,
`run_ubam_output` at `:969`), and precedes `ensure_output_dir` (`:533`), the output-collision
pre-flights (`:480-500`) and every `resolve_clump_layout` call (`:1486`, `:1724`, `:2011`).

| mode | `--phred64 --fastqc` outcome | refused by |
|---|---|---|
| `--hardtrim5` / `--hardtrim3` / `--clock` / `--implicon` | rc=1, nothing written | `Cli::validate()` §#421, `src/cli.rs:1076-1101` — `--fastqc` is rejected in these modes outright, on both base and new |
| `--clump_only` FASTQ out | rc=1, nothing written | new pre-flight |
| `--clump_only --output-format ubam` | rc=1, nothing written | new pre-flight (base: rc=0 — see M1) |
| `--output-format ubam` trim | rc=1, nothing written | new pre-flight (base: rc=0 — see M1) |
| uBAM **input** | rc=1, nothing written | the pre-existing `any_bam` guard at `src/main.rs:368`, which runs first |
| `--passthrough` carrier | rc=1, nothing written | new pre-flight, and correctly so — the carrier output *is* FastQC'd with `phred64` at `src/main.rs:1961-1968` |
| stdin / FIFO / process substitution | rc=1 | `format.rs` restartability check, unchanged |
| truncated `.gz` | rc=1, `incomplete deflate stream`, identical on both | pre-existing sanity check, which runs before the guard |

No bypass found. The `--passthrough` inclusion (`src/main.rs:385`) is tight, not over-broad: the
carrier's own FastQC call receives `cli.phred64`, so a Phred+33 carrier really would be rejected
downstream — base proved it by writing 11 files first.

---

## 5. Gates

- `cargo fmt --all -- --check` — clean.
- `cargo clippy --all-targets --release -- -D warnings` — clean.
- `cargo test --release` — all pass. `integration_fastqc_phred64` runs 10/10 (6 pre-existing + 4 new).

## 6. What I could not verify

- **CI's 1.99 lints.** Local toolchain is `rustc 1.95.0`; 1.99-only lints cannot reproduce here.
- **Whether the 10,000-read window is adequate in practice.** The off-by-one is exact and the fallback
  works, but I have no real Phred+64 dataset whose first 10,000 reads are clean and whose later reads
  are not — all my evidence for that path is synthetic.
- **FastQC-side behaviour is against `fastqc-rust 1.1.1`** (`Cargo.toml:39`), not the `=1.0.1` that
  `CLAUDE.md` documents. Unrelated to this PR, but the "FastQC rejects / accepts" rows above are
  1.1.1's behaviour.
- **Real-world frequency of M3's trimmed-away-tail case.** Constructed, and it reproduces reliably;
  I cannot say how often it occurs in the wild.
