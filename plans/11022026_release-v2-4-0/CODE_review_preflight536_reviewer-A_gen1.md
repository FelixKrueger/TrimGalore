# Code review — PR #536 (`--phred64` pre-flight under `--fastqc`), Reviewer A, gen 1

- **Target**: `fix/phred64-preflight-532` @ `e4caf3c`, diffed against `origin/dev` @ `a4b67ee`
- **Commits**: `ebfb84f` (pre-flight), `e4caf3c` (Solexa wording)
- **Files**: `src/main.rs` (+52), `tests/integration_fastqc_phred64.rs` (+150), `CHANGELOG.md` (+6), `docs/src/content/docs/guide/quality.md` (1 line)
- **Reviewed in place** on the shared tree; no branch switch, no worktree.

## Verdict

The guard itself is correct and correctly placed. Two user-facing claims attached to it are
false, one of them load-bearing for the design decision. One narrow combination produces a
wrong diagnosis. The tests pass for the right reason but pin almost nothing.

---

## F1 — The documented residual is false (High, documentation/comment)

`src/main.rs:26-29` and `docs/src/content/docs/guide/quality.md:40` both claim:

> a file whose only low character sits past this window still fails later, inside FastQC.

It does not. FastQC analyses the **trimmed output**, not the input, and default quality
trimming removes exactly the characters the claim relies on.

**Experiment 1** (low character past the window, at the 3' end):
10,000 reads of 34 × `h` (ASCII 104), then one read of 24 × `h` + 10 × `#` (ASCII 35).

```
trim_galore --phred64 --fastqc <file>        →  exit 0
FastQC report: Encoding = Illumina 1.5, Total Sequences = 10001, Sequence length 24-34
```

The 10 low bases were trimmed off (`34 → 24`), so FastQC never saw ASCII 35.

**Experiment 2** (realistic binned Phred+33, which is the case that matters):
10,000 reads of 34 × `F` (ASCII 70 = Q37 Phred+33), then 200 reads ending in `,` (ASCII 44).

```
trim_galore --phred64 --fastqc <file>        →  exit 0
Reads too short: 10200 (100.0%)
Reads written (passing filters): 0 (0.0%)
FastQC report: >>Basic Statistics  pass
               Encoding = Illumina 1.5
               Total Sequences = 0
```

Under `--phred64`, `F` resolves to Q6, below the default `-q 20`, so the **entire library** is
discarded and FastQC is handed a zero-read file. `fastqc-rust`'s
`phred::detect(NO_QUALITY_SEEN /* 1000 */, true)` returns `Ok(Illumina 1.5)` rather than an
error (`~/.cargo/registry/src/index.crates.io-*/fastqc-rust-1.1.1/src/utils/phred.rs:56-82`),
so nothing fails anywhere. The run exits 0 having silently produced an empty output.

The other half of the claim *is* true: `fastqc-rust` mins `lowest_char` over every record
(`src/modules/basic_stats.rs:329`) and `write_text_report` fails the whole file when
`phred::resolve` errors (`src/modules/basic_stats.rs:395-397`). The flaw is only that the file
FastQC sees is not the file the pre-flight scanned.

**Recommendation.** Drop the residual sentence from both places, or replace it with what is
actually true — a low character past the window is **not** caught; the window is a
cost/benefit bound, not a safety net backed by FastQC. The `CHANGELOG.md` entry already stops
at "The scan covers the first 10,000 reads" and needs no change.

Note this also answers the "is the bound honest?" question in the brief: the bound is honestly
*stated* (10,000, and the code does exactly that), but the justification offered for why it is
safe does not hold.

**Context, not a request to re-litigate the `--fastqc` scoping.** Experiment 2 also shows where
the real damage sits. The same file without `--fastqc` gives:

```
trim_galore --phred64 <realistic mislabelled file>   →  exit 0
Reads too short: 10200 (100.0%);  realistic_late_trimmed.fq is 0 bytes
```

So the characteristic failure of a mislabelled `--phred64` is a silently emptied library, and
that happens on exactly the runs the pre-flight deliberately does not cover. The trimming report
does say `100.0%` too short and the startup NOTE fires, so it is not invisible. If a wider net
is ever wanted without rejecting the `--hardtrim5/3` / `--clump_only` pass-through runs, the
cheap version is to reuse `lowest_quality_char` as a **warning** on every `--phred64` run and
keep the hard refusal where it is.

## F2 — "FastQC rejects such a file either way" is false on every BAM-output path (Medium)

The refusal's closing sentence (`src/main.rs:397`) asserts the counterfactual. On
`--output-format ubam` it is wrong: `fastqc::run` is called with `phred64: false` at
`src/main.rs:2428`, `src/main.rs:2587`, `src/main.rs:2707`, and `src/clump_only.rs:917` /
`src/clump_only.rs:1177`, so FastQC is never told the encoding and cannot refuse.

Verified both combinations are legal with FASTQ input and are refused *only* by the new
pre-flight:

```
trim_galore --phred64 --fastqc --output-format ubam <clean p64>              → exit 0
trim_galore --phred64 --fastqc --output-format ubam <mislabelled>            → refused (new guard)
trim_galore --clump_only --output-format ubam --phred64 --fastqc <clean p64> → exit 0
trim_galore --clump_only --output-format ubam --phred64 --fastqc <bad>       → refused (new guard)
```

This is a genuine *extension* of protection, not a bug — on those paths the pre-flight is now
the only check there has ever been. But the sentence telling the user FastQC would have caught
it is inaccurate there. Suggest trimming it to something path-independent, e.g. "No encoding
makes this file Phred+64."

## F3 — A uBAM `--passthrough` carrier gets the wrong diagnosis (Medium)

The new block (`src/main.rs:384-406`) runs **before** the carrier's own format check at
`src/main.rs:441-450`, and `any_bam` (which gates the older `--phred64` + BAM refusal at
`src/main.rs:367`) is computed from `cli.input` only — the carrier is a separate surface. The
comment at `src/main.rs:441-442` says that check exists so the carrier's format is settled
"before any output is created"; the new guard now opens the carrier as a record source ahead of
it, which is what produces the wrong message.

```
trim_galore --paired --fastqc --passthrough test_files/ubam_test.bam R1 R2
  → Error: --passthrough is not supported with uBAM input in this release. …   (correct)

trim_galore --paired --phred64 --fastqc --passthrough test_files/ubam_test.bam R1 R2
  → Error: --phred64 declares Phred+64 input, but test_files/ubam_test.bam contains
    quality character '?' (ASCII 63) … or it is Solexa / Illumina 1.0 …          (wrong blame)
```

Two problems compound: the user is sent to fix `--phred64` when the real issue is that
`--passthrough` does not accept BAM at all, and the ASCII value quoted is **synthesised by the
reader**, not a byte in the file — `bam_record_to_fastq` emits raw Phred + 33
(`src/bam.rs:1033-1038`), and a BAM with missing qual yields `'!'` (ASCII 33). This is exactly
the wrong-blame the `src/main.rs:367` guard exists to prevent.

**Fix** (one of):
- Move the new block below `src/main.rs:450` — still 80+ lines above
  `naming::ensure_output_dir` (`src/main.rs:533`), so the no-writes property is untouched; or
- skip the carrier in the loop when `detect_input_format(pt)? == UnalignedBam`.

## F4 — The new tests pin very little (Low)

All four tests pass, and they pass for the right reason (the brief's two mutants are real).
But the following mutants survive the whole suite:

| Mutant | Survives because |
| --- | --- |
| `cli.input.iter()` → `.take(1)` | all four new tests are single-input; no `--paired` case |
| drop `.chain(carrier)` | no test gives the carrier a low character the mates lack — the pre-existing `phred64_reaches_the_passthrough_carrier` uses a clean Phred+64 carrier |
| keep only the last record's min (`lowest = Some(low)`) | the refusal fixtures are 1 and 4 reads, all low; no mixed file with the low char in an early read |
| `PHRED64_PREFLIGHT_READS` → any value ≥ 4 | nothing exercises the bound |
| `lowest < 64` → `lowest < 63` | no fixture sits at ASCII 63 |
| `lowest < 64` → `lowest <= 64` | nothing asserts ASCII 64 (`'@'`, Q0) is accepted |

I checked the last two by hand against the built binary: ASCII 63 (`'?'`) is refused, ASCII 64
(`'@'`) is accepted. The code is right; the boundary just is not pinned by a test. Four cheap
additions would close most of this: a `--paired` case with the low character only in R2; a
`--passthrough` case with a Phred+33 carrier behind clean mates; a file whose low character is
in read 1 of 101; and an ASCII-64 accept assertion. `test_files/ubam_test.bam` already contains
ASCII 63 if a fixture is wanted.

## F5 — The "wrote nothing" assertion does not pin the placement (Low)

`tests/integration_fastqc_phred64.rs:259-290` asserts the output directory holds nothing but
the input. `--output_dir` is the pre-existing `tempdir()`, so the assertion also holds if the
guard moves to just *after* `naming::ensure_output_dir`. Pointing `--output_dir` at a
non-existent subdirectory and asserting it was never created would pin what the test name
claims.

## F6 — Two gratuitous clones (Nit)

`src/main.rs:385-386`:

```rust
let carrier = cli.passthrough.iter().cloned();
for input in cli.input.iter().cloned().chain(carrier) {
```

`Option::iter` and `slice::iter` both yield `&PathBuf`, so
`cli.input.iter().chain(cli.passthrough.iter())` chains directly and drops a `PathBuf` clone
per input. The `carrier` binding then becomes unnecessary too.

## F7 — `lowest_quality_char` is unreachable from unit tests (Nit)

It lives in the binary (`src/main.rs:34`), so only process-spawning integration tests can reach
it. Taking `&mut dyn RecordSource` and living in the library would make the min-accumulation and
the window bound unit-testable without a subprocess — which is what F4 needs.

---

## What I verified and found correct

- **Placement.** The guard is at `src/main.rs:384`, ahead of `naming::ensure_output_dir`
  (`src/main.rs:533`), the specialty dispatch (`src/main.rs:587+`), `--clump_only`
  (`src/main.rs:748+`), `run_ubam_output` (`src/main.rs:969`) and `--demux`
  (`src/main.rs:1635`). Nothing between `cli.validate()` (`src/main.rs:322`) and the guard
  writes: I scanned lines 322–533 for `File::create`, `create_dir`, `FastqWriter::` and
  `fs::write` and found none.
- **Scope coherence.** `--fastqc` with `--hardtrim5/3`, `--clock` or `--implicon` is rejected in
  `Cli::validate()` (`src/cli.rs:1078-1101`), which runs before the guard, so the guard can never
  fire in a specialty mode. The deliberate `&& cli.fastqc_requested()` narrowing does what the
  commit says it does.
- **BAM ordering holds.** `--phred64` + uBAM in `cli.input` is rejected at `src/main.rs:367`,
  before the guard, so the guard never opens a BAM from the input list. The carrier is the
  exception — F3.
- **No new restartability hazard.** Every input and the carrier pass `check_restartable_input`
  in `Cli::validate()` (`src/cli.rs:976`, `src/cli.rs:1048`), which rejects FIFOs and sockets,
  and `detect_input_format` additionally calls `reopen_restarts` (`src/format.rs:354`). Trim
  Galore already opens each input three or four times (format detect, sanity check, adapter
  auto-detect, trim pass), so one more open is not a new class of problem; `/dev/stdin` and
  `<(zcat …)` were already refused before this change.
- **Including the carrier is right on the FASTQ path.** The passthrough FastQC call passes
  `cli.phred64` (`src/main.rs:1968`), so the carrier is genuinely subject to FastQC's refusal
  there. The symmetry holds.
- **No false positive on empty quality lines.** `.min()` over an empty `qual` yields `None` and
  the record is skipped (`src/main.rs:38-40`); a file with no quality bytes at all returns
  `None` and is not refused. Verified with a fixture whose fourth line is empty — the run
  proceeded past the guard. `FastqReader::next_record` (`src/fastq.rs:439-469`) does not
  validate `qual.len() == seq.len()`, so this path is reachable.
- **Boundary.** `< 64` accepts ASCII 64 (`'@'`, Q0) — verified; a quality line of all `'@'` also
  does not confuse the parser.
- **`e4caf3c`'s Solexa wording is accurate.** Solexa / Illumina 1.0 is the only encoding that
  legitimately uses ASCII 59–63, and `fastqc-rust` rejects such a file under `--phred64` too
  (`fastqc-rust-1.1.1/src/utils/phred.rs:67-72`).
- **Docs claim "prints a `NOTE:` on every output path"** is true — `src/main.rs:540` is gated on
  `cli.phred64` alone. (Minor: on a refusal the NOTE is not printed, since the guard at 384
  precedes it. The refusal text stands on its own, so this seems fine.)
- `cargo fmt --all -- --check`: clean.
- `cargo test --release`: 834 passed, 0 failed, exit 0 (includes the 10 in
  `integration_fastqc_phred64`).
- `cargo clippy --all-targets --release -- -D warnings`: exit 0, zero warning/error lines, in a
  **brand-new** `CARGO_TARGET_DIR` so nothing was a cache hit. I checked this rather than
  trusting the exit code: a crate with a build script gets one `Compiling trim-galore` status
  line covering all its units rather than per-target `Checking` lines, and
  `libintegration_fastqc_phred64-*.rmeta` is present in the fresh target dir, so the new test
  file really was linted. (An earlier run of mine reported a meaningless exit code because the
  command was piped through `tail`; that result is not what is stated here.)

## What I could **not** verify

- **CI clippy at Rust 1.99.** Local toolchain is `rustc 1.95.0`; the clean cold run above is at
  1.95. Lints introduced after 1.95 (e.g. `byte_char_slices`) cannot reproduce here, so this
  does not prove CI-clean.
- **The Perl 0.6.11 validation matrix.** No Perl `trim_galore` installed. The change adds no
  code to any non-`--phred64` path, so byte-identity is unaffected by construction, but I did
  not run the matrix.
- **Non-unix behaviour.** `non_restartable_kind` returns `None` off unix
  (`src/format.rs:271-273`), so a pipe could in principle reach the pre-flight there. Pre-existing,
  not introduced by this PR.
