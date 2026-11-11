# Plan Coverage Report

**Mode:** B (code vs. plan), post-implementation
**Plan(s):** `plans/08112026_ubam-qname-whitespace/PLAN.md` (v3)
**Date:** 2026-08-12
**Verdict:** INCOMPLETE — 1 item unresolved (a one-line source comment; no code or test gap)

Audited in an isolated worktree at `1fc00f6` (`$TMPDIR/coverage-415-gen1`) with the uncommitted
`git diff` applied and all seven untracked `test_files/` additions copied in and md5-verified
against the working tree. Nothing in `/Users/fkrueger/Github/TrimGalore` was modified.

## Summary

- Total items: 48
- DONE: 41
- PARTIAL: 1
- MISSING: 0
- DEVIATED (all five documented in the plan's Implementation notes): 5
- DEFERRED by explicit maintainer decision (not a gap): 1 — §5 step 8

All three claimed test runs reproduce exactly: `cargo test`, `cargo test --release` and
`cargo nextest run` each report **659** tests, 0 failed, 0 ignored, 0 skipped. The +16 delta is
structurally confirmed (6 new unit tests in `src/bam.rs`, 9 in `tests/integration_ubam.rs`, 1 in
`tests/integration_clump_only_ubam.rs`), so the 643 baseline claim holds. `cargo fmt --all --
--check`, `cargo clippy --all-targets --release -- -D warnings` and `cargo test --no-run` all exit
0. The docs site builds clean and renders the new section.

## Coverage ledger

### §3 Behavior (6 numbered behaviours)

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `reject_whitespace_in_qname` called from `bam_record_to_fastq` after the empty-name check, before the id is built | §3.1 | DONE | `src/bam.rs:933`, between the empty check and `from_utf8` at `:934` |
| 2 | Error names the record, the offending name with whitespace escaped, and a remediation | §3.2 | DONE | Error text matches §5 step 1 **verbatim**, em-dash included; `escape_debug` on a `from_utf8_lossy` string. Record index arrives from the caller wrappers as §5 step 1 predicted — measured: `Error: test_files/ubam_ws_qname_late.bam: BAM record 2` |
| 3 | Each `Z` value scanned in the tag-appending loop before it is appended; space stays legal | §3.3 | DEVIATED | Documented deviation 1. Check sits in `append_tag_type_and_value`'s `Value::String` arm (`src/bam.rs:1098`); the tag name comes from `.with_context(\|\| format!("aux tag '{}'", tag_name))` at the loop call site (`:957`). Behaviour equivalent; space still legal |
| 4 | Otherwise proceed exactly as today; no other behaviour changes | §3.4 | DONE | Diff touches only the two predicates, the two call sites, the `with_context`, and tests |
| 5 | Every BAM input path inherits both checks — SE and paired, FASTQ and uBAM out, `--clump_only` and specialty modes | §3.5 | DONE | Five paths covered by tests. Specialty inheritance verified empirically beyond the plan: `--hardtrim5 5` and `--hardtrim3 5` on `ubam_ws_qname.bam` both exit 1 with the whitespace message. (`--clock`/`--implicon` need an even number of *files*, so a single interleaved BAM cannot reach them — pre-existing, unrelated to #415) |
| 6 | Nothing on the FASTQ-input side is touched; a space in a FASTQ header stays a description, dropped with the existing notice | §3.6 | DONE | No FASTQ-side code changed; pinned by `description_notice_wording_is_fastq_input_only` |

### §3 Edge-case table (12 rows)

| # | Item | Status | Notes |
|---|------|--------|-------|
| 7 | QNAME with a space → error | DONE | Unit (`qname_boundary_bytes_rejected`) + 6 integration cases |
| 8 | QNAME with a tab → error | DONE | Unit; escaping pinned by `qname_error_escapes_invisible_bytes` |
| 9 | QNAME with `\n`, `\r`, `\f` → error | DONE | All three in the unit loop; `\n` also end-to-end |
| 10 | QNAME with `\v` accepted | DONE | Pinned by `qname_non_boundary_bytes_accepted`. Deliberately not integration-tested (documented deviation 5) |
| 11 | QNAME with non-ASCII UTF-8 accepted by this check | DONE | Unit case `"read\u{e9}one"` |
| 12 | QNAME containing `@` accepted by this check | DONE | Unit case |
| 13 | QNAME longer than 254 bytes unchanged | DONE | Predicate tests whitespace only; no length logic added |
| 14 | Leading/trailing whitespace → error, no trimming | DONE | Guard present by construction: `name.iter().any(u8::is_ascii_whitespace)` is position-independent and no `trim` is called |
| 15 | Empty QNAME unchanged, already an error | DONE | `:931-932` empty check preserved *above* the new call |
| 16 | Invalid-UTF-8 QNAME unchanged; the new check runs **before** it so a name that is both reports the whitespace | DONE | Ordering verified: `:933` precedes `from_utf8` at `:934`, and the message uses `from_utf8_lossy`, so an invalid name cannot degrade into an encoding error |
| 17 | A `Z` tag value containing a space accepted | DONE | Unit (`tag_value_space_accepted`) + integration (`clean_tag_value_with_space_still_accepted`, fixture `CB:Z:has a space`) |
| 18 | A `Z` tag value containing `\n`, `\r`, `\f` or `\t` → error | DONE | Unit (`tag_value_framing_bytes_rejected`); `\n` also end-to-end |

### §4 Signature (2 signatures)

| # | Item | Status | Notes |
|---|------|--------|-------|
| 19 | `fn reject_whitespace_in_qname(name: &[u8]) -> Result<()>` with the two-line docstring | DONE | `src/bam.rs:869`; signature and docstring reproduced **verbatim** from §4 |
| 20 | `fn reject_framing_in_tag_value(tag: &str, value: &str) -> Result<()>` with the two-line docstring | DEVIATED | Documented deviation 1 (drops the `tag` parameter) and deviation 2 (docstring rewritten around the derived predicate). `src/bam.rs:889`. Both are free functions and unit-tested per character class, as §4 required |

### §2.8 (folded into scope)

| # | Item | Status | Notes |
|---|------|--------|-------|
| 21 | The `Z`-tag-value predicate is distinct from the QNAME predicate — space accepted, `\t`/`\n`/`\r`/`\f` rejected | DONE | Implemented as `b != b' ' && b.is_ascii_whitespace()`. Documented deviation 2 makes it derived rather than enumerated and **corrects §2.8's own table**, which had justified rejecting `\f` as "breaks FASTQ line framing" |

### §5 Implementation outline (8 steps)

| # | Item | Status | Notes |
|---|------|--------|-------|
| 22 | Step 1 — both predicates, both call sites, the exact error text, `escape_debug`, `from_utf8_lossy` | DONE | Error string is a character-for-character match with the plan |
| 23 | Step 2 — unit tests per character class; accepted classes stay unit-only | DONE | 6 unit tests at `src/bam.rs:1449-1530`; all five rejected QNAME bytes, all four accepted classes, both tag directions, plus the escaping assertion |
| 24 | Step 3 — fixtures and a committed generator, recorded in `test_files/README.md` | DEVIATED | Documented deviation 3: **six** fixtures, not five. All six shapes verified with samtools (see Fixture verification below). `test_files/mk_bam_fixture.py` reproduces both LF fixtures **byte-identically** (md5 match) via the documented `\| bgzip` recipe |
| 25 | Step 4 — nine integration cases, each naming its entry point | DEVIATED | Eight implemented as specified; bullet 8 (acceptance twin) is documented deviation 4 — see item 40 |
| 26 | Step 5 — `--clump_only` guard | DONE | `tests/integration_clump_only_ubam.rs:816`; verified load-bearing |
| 27 | Step 6 — one docs addition in the uBAM **input** section: refused on every output format, with the reason, plus the samtools recipe and its space-only caveat; no `--help` change | DONE | `### Read names containing whitespace` at `docs/src/content/docs/quickstart.md:70`, inside `## Unaligned BAM (uBAM)` (`:55`) and above `## Common combinations` (`:84`). States "**on every output format**", carries the `FS="\t"` awk recipe and the tab/newline caveat. No `--help` change in the diff. `npx astro build` exits 0 and the section renders in `dist/quickstart/index.html` |
| 28 | Step 7 — one `CHANGELOG.md` entry in the file's existing format | DONE | First entry under `### Unreleased` / `#### Changes`, correct link format, and it names the behaviour change. Nit: it is **two** sentences where §5 step 7 said one — see Gaps |
| 29 | Step 8 — file the two follow-up GitHub issues | DEFERRED | **By decision, not oversight.** Recorded in the plan's "Still open, by decision" paragraph; filing is outward-facing and reserved for the maintainer. Not counted as a gap |

### §2.6 decision record

| # | Item | Status | Notes |
|---|------|--------|-------|
| 30 | "A test pins that a whitespace-QNAME BAM now errors, **and a comment records why the wording is unconditionally correct**" | PARTIAL | The test exists (item 48). The comment does not exist at the notice — see Gaps. This is the one gap I found that is **not** among the five documented deviations |

### §9 Validation (18 V-rows)

| # | Item | Status | Notes |
|---|------|--------|-------|
| 31 | V1 — the defect exists before the fix | DONE | **Independently reproduced**, not taken on trust: with both call sites disabled, `--output-format ubam --preserve-tags CB` on `ubam_ws_qname.bam` gives exit 0, QNAME `name`, tags `[]` |
| 32 | V1b — LF QNAME, FASTQ out, pre-fix | DONE | **Independently reproduced**: exit 0, **13 lines for 3 records**. (The plan's V1c, LF in a `CB:Z` value, also reproduces at exit 0 / 13 lines) |
| 33 | V2 — refused at the sanity-check arm, nothing written | DONE | `ws_qname_refused_at_sanity_check_entry_point`; asserts non-zero exit, both the message and the offending name, and an empty `*_trimmed*`/`*_val*` glob |
| 34 | V3 — refused in the trimming loop, uBAM out, **partial `*_trimmed.bam` present** | DONE | Asserts `BAM record 2` (proving it is not the sanity check) and that the partial exists with non-zero length — the A7 residue asserted rather than wished away |
| 35 | V4 — the widened scope is real: same fixture, FASTQ out, refused | DONE | Reads plain `.fq`; no `*.fq.gz` glob anywhere in the new tests |
| 36 | V5 — the `:380` threaded reader is covered | DONE | `--cores 2` arm |
| 37 | V6 — the `:261` de-interleaver is covered, with its wrapper shape | DONE | Asserts `BAM record 3` — the first read of the second pair, which only the de-interleaver can report (the sanity check sees record 1 only). Not the vacuous pass v2 would have produced |
| 38 | V7 — LF in QNAME end-to-end | DONE | Asserts the escaped `readB\nEVIL` in the message |
| 39 | V8 — the tag-value scan | DONE | Asserts both `aux tag` (the `with_context` from deviation 1) and the escaped `AAA\nCCC` |
| 40 | V9 — refusal is not over-broad; clean-fixture acceptance twin | DEVIATED | Documented deviation 4. Verified real: `preserve_tags_roundtrip_matches_golden` (`tests/integration_ubam.rs:173`) round-trips `ubam_test_with_tags.bam` against a committed golden and spot-checks `CB:Z:ATCGATCG-1` / `UB:Z:GCTAGCTA` — both whitespace-free, so the new tag predicate cannot disturb it, exactly as the notes claim. Not duplicated |
| 41 | V10 — every QNAME boundary byte | DONE | space, `\t`, `\n`, `\r`, `\f` |
| 42 | V11 — accepted classes stay accepted, unit-only | DONE | plain, `\v`, `@`, non-ASCII |
| 43 | V12 — the tag predicate is not the QNAME predicate | DONE | Both directions asserted |
| 44 | V13 — `--clump_only` inherits the refusal | DONE | On the space fixture, per the V-row |
| 45 | V14 — the multi-input partial-run shape is **recorded** | DONE | Recorded in plan §7 as the row asks. I also measured it: two SE inputs with the second bad leaves input 1 complete (`ubam_test_trimmed.fq` + both reports) and input 2 partial (`ubam_ws_qname_late_trimmed.fq`, 91 B, no report). The documented shape is accurate |
| 46 | V15 — §2.6's unreachability | DONE | `description_notice_wording_is_fastq_input_only`: FASTQ input succeeds and stderr still contains `FASTQ header text` |
| 47 | V16 — no regression | DONE | 659/659 in all three runners; clippy and fmt exit 0. See Summary |
| 48 | V17 — compile check after the signature change | DONE | `cargo test --no-run` exit 0 |

## Gaps (detail)

### Item 30: §2.6's source comment (PARTIAL — the only undocumented gap)

**Expected:** §2.6 resolves the issue's second half ("the #406 notice says *FASTQ* header text even on
a uBAM-input run") by deleting the possibility rather than rewording, and specifies two artifacts:
"A test pins that a whitespace-QNAME BAM now errors, **and a comment records why the wording is
unconditionally correct**."

**Found:** The test is present and passing. `emit_description_dropped_once`
(`src/bam.rs:1068-1082`) is untouched by the diff and carries no comment. The reasoning lives only
inside the body of `description_notice_wording_is_fastq_input_only`
(`tests/integration_ubam.rs:479-483`).

**Gap:** A maintainer reading the notice string `"...so FASTQ header text after the first space is
not carried into uBAM output"` gets no pointer to why "FASTQ" is now unconditionally accurate, and
so may re-open the #406 wording question or "fix" it per-direction — the outcome §2.6 explicitly
decided against. One line above `emit_description_dropped_once` closes it. Note this also sits
squarely in the project's own comment convention ("put an invariant on the thing it constrains,
not in the function that relies on it"): the invariant is about the notice, and the comment is in
the test.

### Item 28: CHANGELOG sentence count (nit, not a coverage gap)

§5 step 7 said "**One sentence**, in the file's existing format", and the global convention is one
sentence with no exceptions. The entry is two: a bolded headline sentence with the issue link, then
a second sentence covering the read-side scope, the `Z`-tag-value half, and the docs pointer. It is
far shorter than its `#398` neighbour (which runs to three paragraphs), so it is consistent with
the file as it stands, and everything the plan asked for is present. Flagged only so the choice is
deliberate.

### Item 29: §5 step 8 (deferred, not missing)

The two follow-up issues — partial output on any per-record bail, and the opaque noodles write
error — are not filed. The plan's Implementation notes record this as reserved for the maintainer
because filing is outward-facing. Recorded here as DEFERRED so it stays visible; the residue it
covers ships documented in §7 and asserted by V3, so it cannot change unnoticed.

### Minor documentation imprecision (no action required)

`test_files/README.md`'s verify line is commented `# shows readB \n EVIL`. Run as written
(`samtools view … | sed -n 2p | cut -f1 | od -c`) it shows `readB \n` — the embedded newline puts
`EVIL` on line 3, so `sed -n 2p` cannot reach it. The command does demonstrate the raw LF, which is
its purpose.

## Fixture verification

All six verified with samtools 1.21 against the shapes §5 step 3 and `test_files/README.md`
specify. The multi-record fixtures all place the offender past record 1, which is what makes V3–V8
non-vacuous.

| Fixture | Verified shape | Matches spec |
|---|---|---|
| `ubam_ws_qname.bam` | 1 rec, QNAME `name with space`, `CB:Z:AAACCC` | yes |
| `ubam_ws_qname_late.bam` | 3 recs, clean `readA` first, `readB with space` at rec 2 | yes |
| `ubam_ws_qname_paired.bam` | 4 recs = 2 interleaved pairs, flags 77/141, offender at rec 3 (`pairB with space`) | yes |
| `ubam_tagvalue_space.bam` | 1 rec `readA`, `CB:Z:has a space` | yes (the deviation-3 sixth fixture) |
| `ubam_lf_qname.bam` | 3 recs, raw LF inside rec 2's QNAME (`readB\nEVIL`) | yes |
| `ubam_lf_tagvalue.bam` | 3 recs, raw LF inside rec 2's `CB:Z` value (`AAA\nCCC`) | yes |
| `mk_bam_fixture.py` | regenerates both LF fixtures **md5-identical** to the committed bytes | yes |

## Test verification

All 16 new tests plus the reused acceptance twin, observed by name in the `cargo test` output.

The "load-bearing" column comes from a tripwire: I disabled both call sites in my worktree
(`reject_whitespace_in_qname(name)?` and `reject_framing_in_tag_value(s)?`), rebuilt, and re-ran the
two integration binaries with `--no-fail-fast`. **All 8 refusal tests went red** and the 2
acceptance tests correctly stayed green. `src/bam.rs` was then restored and md5-verified against
the working tree, and the binaries re-run green (16 + 17 passed). This proves the suite can fail
for the reason it claims — the failure mode the plan's own §11 item 5 flags as v2's defect.

| Test | File | Status | Load-bearing |
|------|------|--------|--------------|
| `bam::tests::qname_boundary_bytes_rejected` | `src/bam.rs:1449` | PASS | predicate-level |
| `bam::tests::qname_non_boundary_bytes_accepted` | `src/bam.rs:1466` | PASS | predicate-level (the §2.7 `\v` pin) |
| `bam::tests::qname_error_escapes_invisible_bytes` | `src/bam.rs:1484` | PASS | predicate-level |
| `bam::tests::tag_value_framing_bytes_rejected` | `src/bam.rs:1495` | PASS | predicate-level |
| `bam::tests::tag_value_space_accepted` | `src/bam.rs:1511` | PASS | predicate-level |
| `bam::tests::read_side_emitter_rejects_newline_in_string_tag` | `src/bam.rs:1518` | PASS | predicate-level (covers deviation 1's moved call site) |
| `ws_qname_refused_at_sanity_check_entry_point` | `tests/integration_ubam.rs:273` | PASS | yes — red without the fix |
| `ws_qname_refused_in_trimming_loop_ubam_out_leaves_partial` | `tests/integration_ubam.rs:304` | PASS | yes — red without the fix |
| `ws_qname_refused_in_trimming_loop_fastq_out` | `tests/integration_ubam.rs:333` | PASS | yes — red without the fix |
| `ws_qname_refused_via_threaded_reader` | `tests/integration_ubam.rs:357` | PASS | yes — red without the fix |
| `ws_qname_refused_via_interleaved_deinterleaver` | `tests/integration_ubam.rs:378` | PASS | yes — red without the fix |
| `lf_in_qname_refused_end_to_end` | `tests/integration_ubam.rs:407` | PASS | yes — red without the fix |
| `lf_in_tag_value_refused_end_to_end` | `tests/integration_ubam.rs:429` | PASS | yes — red without the fix |
| `clean_tag_value_with_space_still_accepted` | `tests/integration_ubam.rs:452` | PASS | acceptance — green either way, correctly |
| `description_notice_wording_is_fastq_input_only` | `tests/integration_ubam.rs:479` | PASS | acceptance — green either way, correctly |
| `clump_only_ubam_refuses_whitespace_qname` | `tests/integration_clump_only_ubam.rs:816` | PASS | yes — red without the fix |
| `preserve_tags_roundtrip_matches_golden` (reused for V9) | `tests/integration_ubam.rs:173` | PASS | acceptance twin, pre-existing |

### Commands run

| Command | Result |
|---|---|
| `cargo test --no-run` | exit 0 (V17) |
| `cargo test` | exit 0 — 15 binaries, **659 passed**, 0 failed, 0 ignored |
| `cargo test --release` | exit 0 — **659 passed** |
| `cargo nextest run` | exit 0 — **659 tests run: 659 passed, 0 skipped** |
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --release -- -D warnings` | exit 0 |
| `npx astro build` (docs, in-worktree) | exit 0; `dist/quickstart/index.html` contains the new section and the awk recipe |

Two harness false negatives were caught and corrected rather than reported as findings: a `tail -80`
in my first log wrapper truncated 13 of the 15 per-binary `test result:` lines, and a `for` loop
relying on word-splitting of `$mode` (which zsh does not do) made the specialty modes look like
they had exit-2 CLI errors. Both were re-run correctly before anything was recorded above.

## Verdict

**INCOMPLETE — 1 item unresolved.** The implementation is substantively complete: every behaviour,
every edge-case row, both signatures, all 18 V-rows and 7 of the 8 outline steps are satisfied, the
test counts and all six quality gates verify exactly as claimed, and the five deviations the plan
documents are all real and correctly described. The pre-fix evidence (V1, V1b, V1c) reproduced
independently, and a disable-the-fix tripwire confirms all eight refusal tests are load-bearing.

One thing remains for the implementer:

1. **Add the §2.6 comment.** One line above `emit_description_dropped_once` (`src/bam.rs:1068`)
   recording why the notice's "FASTQ header text" wording is now unconditionally correct — no space
   can reach `parse_name_and_data` from BAM input once whitespace QNAMEs are refused, so the notice
   is reachable only from FASTQ input. §2.6 specified this alongside the test; the test landed and
   the comment did not. It is the only requirement in the plan that is neither met nor documented
   as a deliberate deviation.

Two items need a decision rather than work, and neither blocks:

2. **§5 step 8's two follow-up issues are unfiled by design**, reserved for the maintainer. Nothing
   in the code depends on them; the residue they cover is asserted by V3.
3. **The CHANGELOG entry runs to two sentences** where the plan and the global convention say one.
   Consistent with its neighbours in the file, so this is a style call, not a gap.
