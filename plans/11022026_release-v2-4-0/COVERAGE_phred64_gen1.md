# Plan Coverage Report

**Mode:** B (code vs. draft)
**Plan(s):** `DRAFT_doc-sweep-commit.md` → `1b16bed`; `DRAFT_phred64-fastqc-fix.md` → `8d5567a`
**Repo / branch:** `/Users/fkrueger/Github/TrimGalore`, `chore/fastqc-rust-1.1.1` @ `8d5567a`
**Date:** 2026-11-05
**Verdict:** INCOMPLETE — 2 items unresolved (1 silent drop, 1 inaccuracy in the new copy), plus 2 uncommitted-work flags

## Summary

- Total items: 26
- DONE: 20
- DEVIATED (all documented / intentional): 6
- PARTIAL: 0
- MISSING: 0
- Gaps found outside the drafts' edit lists: 2 (G1, G4)
- Defects in the applied copy itself: 1 (G2)

Every enumerated edit in both drafts landed. The audit's findings are all in
what the drafts did **not** enumerate, plus one factual error inside text the
drafts specified verbatim.

## Verification performed

| Check | Result |
|---|---|
| `cargo check --all-targets` | clean, 0 warnings |
| `cargo check --example fastqc_only` (forced recompile after `touch`) | exit 0 — the untracked example compiles against the 5-arg signature |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0 |
| `cargo fmt --all -- --check` | **exit 1** — one diff, in the untracked `examples/fastqc_only.rs`. All tracked files clean |
| `cargo test` (full suite) | **826 passed, 0 failed, 0 ignored** (824 baseline + the 2 new cases) |
| `cargo test --test integration_fastqc_phred64` | 2 passed |
| CI `validation` smoke-test assertions re-run against a real 1.1.1 report | `summary.txt`, `fastqc_data.txt`, `Images/`, `Icons/` all PRESENT — CLAUDE.md:67's description still true |

## Coverage ledger — `DRAFT_doc-sweep-commit.md` → `1b16bed`

| # | Item | Status | Notes |
|---|------|--------|-------|
| 1 | `src/cli.rs:340` `--fastqc` help → 0.13.0 | DONE | verbatim |
| 2 | `src/fastqc.rs:4-9` module doc → 0.13.0, Docker clause dropped | DONE | verbatim |
| 3 | `README.md:35` → 0.13.0 | DONE | verbatim |
| 4 | `docs/.../quickstart.md:53` → 0.13.0 | DONE | verbatim |
| 5 | `docs/.../guide/outputs.md:149` → 0.13.0 | DONE | verbatim |
| 6 | `docs/.../guide/overview.md:10` → 0.13.0 | DONE | verbatim |
| 7 | `Cargo.toml:36-42` pin rationale rewritten | DONE | verbatim |
| 8 | `Cargo.toml:43-44` noodles comment `1.0.1` → `1.1.1` | DONE | verified: 1.1.1 still resolves `noodles` to 0.88.0, so the claim is true as well as current |
| 9 | `CLAUDE.md:67` pin + CI-scope clause | DONE | verbatim |
| 10 | `AGENTS.md:67` identical text | DEVIATED | applied to the working tree, **not committable** — `AGENTS.md` is untracked. See G5b |
| 11 | `.github/dependabot.yml:50-56` rationale rewritten | DONE | verbatim |
| 12 | `src/cli.rs:348-351` `--fastqc_args` help gains `--png`, `--svg` no-op note | DONE | text verbatim — **but the text is wrong**, see G2 |
| 13 | `src/fastqc.rs:87-89` `Flags handled:` list | DONE | verbatim — see G2 |
| 14 | `docs/.../guide/outputs.md:155` flag subset | DONE | verbatim — see G2 |
| 15 | `tests/integration_preflight_tripwire.rs:29-32` stale `--svg` claim | DONE | verbatim. The new wording ("neither `--png` nor the no-op `--svg` writes an image file beside the primary artefacts") is **accurate** — images live inside the `.zip`, not beside it |
| D | `test_files/README.md:14` deliberately deferred | DONE | confirmed untouched by both commits (last touched by `4e4a762`, no working-tree diff) and confirmed **still inaccurate**: its parenthetical describes the 0.12.1 heuristic, which 0.13.0 changed, and the `--phred64` fix now means TrimGalore declares the encoding rather than letting it be guessed |

### Independent enumeration of every remaining `0.12.1` in the repo

Excluding `plans/`, `target/`, `target_repro_b/`, `optimus_prime/`, `SESSION_HANDOFF*`:

| Site | Should it have changed? |
|---|---|
| `Cargo.lock:899` | **No** — that is `linux-raw-sys 0.12.1`, unrelated to FastQC |
| `CHANGELOG.md:1398`, `:1420` | **No** — historical v2.1.0 entries |
| `docs/.../reference/changelog.md:654`, `:676` | **No** — same entries on the synced page |
| `test_files/README.md:14` | Yes, eventually — **deliberately deferred**, see D and the draft's "Deferred" section |

The draft's count of **ten** sites asserting `FastQC 0.12.1` is correct, and all
ten are now 0.13.0. No eleventh site was missed. A `grep` for `FastQC 0.` across
all live files returns only 0.13.0 plus the deferred fixture note; no workflow,
Dockerfile or recipe asserts a FastQC version number.

### `--png` / `--svg` coverage

`--png` is now named at all three live sites that list the `--fastqc_args`
subset (`src/cli.rs:350`, `src/fastqc.rs:91`, `docs/.../outputs.md:155`), plus
the test module comment. The only remaining `--svg` mentions that do **not**
name `--png` are `CHANGELOG.md:1427` and `docs/.../reference/changelog.md:683`
— both historical v2.1.0 entries, correctly left alone.

No remaining `--svg` mention is misleading. The **new `--png` wording is**, which
is G2.

## Coverage ledger — `DRAFT_phred64-fastqc-fix.md` → `8d5567a`

| # | Item | Status | Notes |
|---|------|--------|-------|
| 16 | `fastqc::run` gains `phred64: bool` | DONE | `src/fastqc.rs:61`; set as `phred64` in the `FastQCConfig` literal at `:67` |
| 17 | `run()` doc line for `phred64` | DEVIATED | wording changed to match the dropped guard (correct). But it landed in the doc block above `pub fn resolved_threads` — see G3 |
| 18 | `.bam` guard inside `run()` | DEVIATED (intentional #2) | **Confirmed real and correctly reasoned.** `fastqc-rust-1.1.1/src/runner.rs:396` resolves `seq_file.known_phred_encoding().or(config.phred64.then_some(PHRED64))`, so a format's own declared encoding wins. The guard would have been inert |
| 19 | `src/main.rs` — all 7 direct call sites | DONE | FASTQ: `:1566` (SE), `:1884`/`:1891` (PE), `:1904` (passthrough) pass `cli.phred64`. BAM: `:2364`, `:2523`, `:2643` pass literal `false` — a consequence of dropping the guard, and consistent with it |
| 20 | `clump_only.rs` — 4 fns gain `phred64: bool` | DEVIATED (intentional #1) | **Confirmed.** Only `clump_only_single` (`:286`) and `clump_only_paired` (`:434`) took the parameter, placed immediately after `fastqc_args` as specified. `clump_only_single_to_bam` (`:792`) and `clump_only_paired_to_bam_one_pair` (`:954`) already carry `input_phred_offset: u8` and pass literal `false`, each with a one-line comment ("BAM QUAL is raw Phred, so no ASCII offset is declared"). Correctly done |
| 21 | `clump_only.rs` — 5 inner `fastqc::run` calls | DONE | `:421`, `:612`, `:613` pass `phred64`; `:917`, `:1177` pass `false` |
| 22 | `src/main.rs` — 5 `clump_only_*` dispatch sites | DONE | Independently mapped: `:689` (paired, FASTQ) and `:726` (single, FASTQ) gained `cli.phred64`; `:771`, `:834`, `:875` call the `_to_bam` functions, which correctly take no such argument |
| 23 | `Cargo.toml` `zip` dev-dependency | DEVIATED (intentional #4) | **Confirmed correct.** `zip = { version = "9", default-features = false, features = ["deflate", "time"] }` is a byte-for-byte match for `fastqc-rust-1.1.1/Cargo.toml:124-130`. `Cargo.lock`'s only change is `+ "zip",` under `trim-galore`'s dependency list — **no new packages**. The rationale is documented in both the Cargo.toml comment and the commit message |
| 24 | `tests/integration_fastqc_phred64.rs` | DONE + EXCEEDED | The draft left the BAM companion case out ("needs a check I have not run"). The commit **added it** (`ubam_output_is_not_told_about_phred64`), so both halves of the contract are pinned. The commit message records a mutation check in both directions |
| 25 | `clump_only.rs` unit-test call sites | DONE | twelve in-module test call sites updated for the new arity |
| 26 | `--fastqc` on specialty modes needs no change | DONE | confirmed: no `fastqc::run` in `specialty.rs`; `cli.rs:2397` holds the #421 refusals |

### Every `fastqc::run` call site in the tree

Exhaustive `grep` over `src/`, `tests/`, `examples/` (no `benches/` exists):

- `src/main.rs` — 7 sites, all updated
- `src/clump_only.rs` — 5 sites, all updated
- `examples/fastqc_only.rs:11` — 1 site, updated (untracked; see G5a)

**13 call sites, 0 missed.** Note the draft says "all 11 call sites" twice; its
own table lists 12. An arithmetic slip in the draft, not a code gap.

## Gaps (detail)

### G1 — `CLAUDE.md:58` and `AGENTS.md:58` still say `=1.0.1` (HIGH)

**Expected:** the doc sweep's purpose was to remove stale `fastqc-rust` /
FastQC version assertions. The draft enumerated `CLAUDE.md:67` and
`AGENTS.md:67` only.

**Found:**

```
- `fastqc.rs` — bundled FastQC integration via the `fastqc-rust` crate
  (exact-pinned to `=1.0.1`; treat upstream bumps as deliberate test events).
```

This is the module-map bullet, nine lines above the "Bundled FastQC dependency"
section the draft did fix. The committed `CLAUDE.md` therefore now
**contradicts itself**: line 58 says `=1.0.1`, line 67 says `=1.1.1`.

**Gap:** update `=1.0.1` → `=1.1.1` at `CLAUDE.md:58` (tracked, needs a commit)
and `AGENTS.md:58` (untracked). This is the one genuine silent drop — the site
was never enumerated, so it was never considered and dropped, which is why it
survived both review passes on the draft.

### G2 — the new `--png` copy is inaccurate (MEDIUM, user-visible `--help`)

**Expected:** accurate `--fastqc_args` documentation. The drafts' whole premise.

**Found:** the text was applied verbatim at all three sites, so this is **not a
coverage gap** — the specified text is itself wrong.

`fastqc-rust-1.1.1/src/config.rs:49` documents the field the flag sets as:

```rust
/// Embed PNG rather than SVG charts in the HTML report.
pub png_output: bool,
```

I generated a real report with neither flag passed. The archive contains
**both formats for every chart** — `per_base_quality.svg` *and*
`per_base_quality.png`, and so on for all nine. `--png` does not cause PNGs to
be written; it only changes which of the two the HTML `<img>` tags reference.

So `src/cli.rs:350-351` — which a user reads in `--help` — presents `--png`
alongside `--nogroup` and `--nofilter` as a flag "translated to the bundled
engine", and explains `--svg` as a no-op "since SVG graphs are the default".
The natural reading is "`--png` switches the output to PNG files"; the truth is
"both file sets are always written; `--png` switches which one the HTML
embeds". Same wording at `src/fastqc.rs:91` and `docs/.../outputs.md:155`.

**Gap:** reword the three sites to say what `--png` does — e.g. "`--png` embeds
the PNG charts in the HTML report instead of the SVGs; both are written to the
archive either way." Flagging rather than deciding the phrasing.

(`tests/integration_preflight_tripwire.rs:31` is unaffected — "writes no image
file beside the primary artefacts" is true for both flags.)

### G3 — the new `phred64` doc line documents the wrong function (LOW)

**Expected:** draft §1 — "`run()` gains the parameter and the BAM guard",
showing the `phred64` doc line added to `run()`'s parameter list.

**Found:** `src/fastqc.rs:36-37`. The doc block headed
`/// Run FastQC on a single trimmed output file.` — which documents
`output_path`, `fastqc_args`, `output_dir`, `cores` and now `phred64` — sits
above `pub fn resolved_threads(fastqc_args, cores)`, not above `pub fn run`.
`pub fn run` has **no doc comment at all**.

**Attribution:** pre-existing. `resolved_threads` was inserted between the block
and `run` by `fdade7e` (2026-08-17, `#446`). `8d5567a` did not cause it, but
added the new line into the misplaced block, so rustdoc now renders a `phred64`
parameter description on a public function that has no such parameter.

**Gap:** move the parameter block down onto `pub fn run` and give
`resolved_threads` its own two lines. Cosmetic, no behaviour impact, worth a
line in whatever commit fixes G1.

### G4 — no CHANGELOG entry for any of the three FastQC commits (MEDIUM, out-of-draft)

**Expected:** neither draft lists a CHANGELOG edit, so strictly this is out of
scope. Raising it because it is a release-readiness gap on a release branch.

**Found:** `CHANGELOG.md` was touched by **none** of `d5f41c5`, `1b16bed`,
`8d5567a`. `### Unreleased` is empty. The `### Version 2.4.0 (Release on
3 November 2026)` section was finalized by `59993c0` on 2026-11-03, two days
before this work.

So a user deciding whether to upgrade is told nothing about: the bundled engine
moving from Java FastQC 0.12.1 to 0.13.0; `--svg` becoming an accepted no-op;
`--png` becoming an accepted flag; or — the behaviour change — a
`--phred64 --fastqc` report that previously came out right and, between
`d5f41c5` and `8d5567a`, came out 31 too high.

**Gap:** one sentence per item under `### Unreleased`, or folded into the 2.4.0
section if it has not shipped. Per the house rule, one sentence each.

### G5 — uncommitted work

**a) `examples/fastqc_only.rs`** (untracked). Edited to the 5-arg signature and
**confirmed to compile** (forced recompile, exit 0). Two flags:

- It is **not rustfmt-clean** — the call now exceeds 100 columns and
  `cargo fmt --all -- --check` exits 1 on it. CI is unaffected (untracked files
  are not in the checkout), but the local `fmt` gate is now red for anyone with
  the file present, and it is the only `fmt` diff in the tree.
- Being untracked, the signature fix **will be lost** on a fresh clone or
  worktree, where the example then fails to compile.

**b) `AGENTS.md`** (untracked). Working-tree content is identical to the
committed `CLAUDE.md` except the two intentional lines (`# AGENTS.md` vs
`# CLAUDE.md`, and "Codex" vs "Claude Code"), so edit 10 **did** land. The
divergence risk is structural: `AGENTS.md` is a verbatim mirror of a tracked
file but is itself untracked, so every `CLAUDE.md` commit silently widens the
gap for anyone who does not also have the local copy — and a fresh checkout has
no `AGENTS.md` at all. Note G1 affects line 58 of **both** files.

## Observations (no action implied)

- **The draft's stated reason for sparing the docs changelog page is false.**
  It says "the docs changelog page is synced". `docs/.../reference/changelog.md`
  was last synced on 2026-07-25 (`4f86c0f`), has no `### Version 2.4.0` section
  at all, and its `### Unreleased` still holds `--clump_only` content. The
  *decision* to leave its 0.12.1 lines alone is still right — they are
  historical v2.1.0 entries either way — but it holds on different grounds than
  the ones given.
- **`src/fastqc.rs:101-102`** carries "`--svg` is kept as a no-op so the flag
  that used to select them is not an error" — historical narrative of the kind
  CLAUDE.md's comment rule forbids. Introduced by `d5f41c5`, outside both
  audited commits, but three lines from text the sweep did touch.
- On-disk the docs directory is `Docs/` while git tracks `docs/`; `core.ignorecase
  = true` hides it. Pre-existing, unrelated to these commits, and the reason the
  drafts' `Docs/...` paths and the commit's `docs/...` paths both resolve.

## Verdict

**INCOMPLETE.** Both drafts were implemented faithfully — every one of the 26
enumerated items landed, the four declared deviations are each real, correctly
done and documented, the deferred fixture note is genuinely untouched and
genuinely still wrong, and 826 tests pass with clippy clean. Nothing was
silently dropped from the drafts' edit lists.

What remains:

1. **G1 — `CLAUDE.md:58` and `AGENTS.md:58` still say `=1.0.1`.** The committed
   `CLAUDE.md` self-contradicts. This is the one real silent drop; the site was
   never enumerated in the draft, which is why no review caught it.
2. **G2 — reword `--png` at `src/cli.rs:350-351`, `src/fastqc.rs:91` and
   `docs/.../guide/outputs.md:155`.** `png_output` only switches which image
   format the HTML embeds; both SVG and PNG files are written to the archive
   regardless. The user-visible `--help` text implies otherwise. The drafts
   specified this text verbatim, so it is a defect in the corrective copy rather
   than a coverage gap — but it is the same class of inaccuracy the sweep
   existed to remove.
3. **G4 — add a CHANGELOG entry** for the 0.13.0 retarget, the `--svg`/`--png`
   flag change, and the `--phred64` FastQC fix. One sentence each.
4. **G5 — commit or discard** `examples/fastqc_only.rs` (and `rustfmt` it; it is
   the only `fmt` diff in the tree), and decide whether `AGENTS.md` should be
   tracked.
5. **G3 — optional**, fold into the G1 commit: move the `run()` parameter doc
   block off `pub fn resolved_threads`.

None of these is a behaviour defect. G1, G2 and G4 are all things a user or a
future agent reads and is misled by, on a release branch.
