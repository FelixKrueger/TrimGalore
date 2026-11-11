# Plan Coverage Report

**Mode:** B (code vs. plan; no `IMPL.md` exists — `PLAN.md` §5 is the task-ordered
implementation plan, audited together with §3's assertions, §5's case table, §5 step 5's
self-tests, §9's V1–V10 and §8's A1–A11)
**Plan(s):** `plans/08102026_preflight-tripwire/PLAN.md` (v2)
**Code audited:** branch `feature/414-preflight-tripwire`; `src/io.rs` (+68, uncommitted)
and `tests/integration_preflight_tripwire.rs` (untracked)
**State pinned to:** `tests/integration_preflight_tripwire.rs` 1294 lines / 51 tests,
sha256 `148d0146…`; `src/io.rs` sha256 `714b36ee…`
**Date:** 2026-08-10
**Verdict:** **INCOMPLETE — 4 items unresolved**

No unresolved item is a false-pass risk in the suite as it stands today. All four are
either a documentation gap or a decay path inside the anti-decay mechanism — which is
the class the plan itself calls fatal, so they are reported rather than waived.

**The test file changed during this audit.** I first audited it at 1263 lines / 50 tests;
Reviewer A then applied four fixes (F1–F4, `CODE_REVIEW_A.md` §4) taking it to 1294 / 51.
Everything below was re-verified against the pinned state above. F1 resolved one of my five
gaps outright — see **G2 (resolved)**.

## Summary

- Total items: 107
- DONE: 99
- PARTIAL: 4
- MISSING: 1
- DEVIATED: 3

Verified by execution, not by reading §13 — all against the pinned state: `cargo test`
**641 passed**; `cargo test --test integration_preflight_tripwire` **51/51** (debug, 10.6 s);
`cargo test --release --test …` **51/51** (1.4 s); `cargo clippy --all-targets --release
-- -D warnings` clean; `cargo fmt --all -- --check` clean; `cargo test --no-run` clean.
Before the F1–F4 fixes the same commands gave 640 / 50 / 50, and `cargo nextest run` gave
**640/640** in 17.8 s (12 separate release runs of the tripwire in total, including one in a
clean-room target dir, one serial, and three at `--test-threads=32` — no flakiness).

### Tree state at audit time — read this first

`git diff --stat dev` → `src/io.rs | 68 +++…`, one file, matching §13's stated shape.
`grep -c MUTANT src/main.rs` → 0.

My **first** release run of the suite failed 11 tests. That was not flakiness and not a
defect: Reviewer B had a live negative-control line in the shared tree
(`std::fs::write(input.with_extension("MUTANT"), b"x")` at `src/main.rs:1014`, recorded in
`CODE_REVIEW_A.md` §0), and the tripwire went red on exactly the 11 `se_trim_*` cases —
the arm that writer sits on. It was reverted before my second run. Two things follow:

1. **This is an unplanned, independent V2-class negative control, and the tripwire caught
   it.** A genuinely new writer on site `1006` produced the intended failure naming the
   unplanned paths, the site, and the candidate list. I observed the red independently of
   the reviewer who planted it.
2. **A concurrent agent was mutating `src/main.rs` during this audit.** Every result above
   was re-run against the clean tree afterwards; treat any earlier verification in this
   session that was not re-run as suspect.

---

## Coverage ledger

### §5 step 1 — the `src/io.rs` hook

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 1 | `DUMP_PLANNED_ENV` pub const beside `collision_key` | §5.1, §4 | DONE | `io.rs:87`, value `TRIM_GALORE_DUMP_PLANNED_OUTPUTS` |
| 2 | `DUMP_PLANNED_PREFIX` pub const | §5.1, §4 | DONE | `io.rs:91`, `"TRIM_GALORE_PLANNED\t"` |
| 3 | `dump_planned_lines(&[PlannedOutput], &str) -> Vec<String>`, pure | §4, §5.1 | DONE | `io.rs:100`; signature matches §4 exactly; no env access, no I/O |
| 4 | `#[track_caller]` on `preflight_output_collisions` | §3.1.1, §5.1 | DONE | `io.rs:172` |
| 5 | Gated dump **before** `input_keys`, per §5.1's snippet | §3.1.1, §5.1 | DONE | `io.rs:178-184`, `input_keys` at `:186`; code is the snippet verbatim |
| 6 | No `cfg` gate (A4 — the most important assumption) | §3.1.4, A4 | DONE | No `#[cfg]` on hook, consts or helper; release run dumps (verified by hand) |
| 7 | Paths dumped as constructed; `OutputSource` omitted | §3.1.5 | DONE | `.map(\|(path, _)\| …)` discards the source |

### §5 step 2 — the `src/io.rs` unit test

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 8 | `dump_planned_lines_renders_one_line_per_candidate` | §5.2 | DONE | `io.rs:977`; pins wire format, tab count, and that an empty list dumps nothing |
| 9 | **No** `env::set_var` (preserves `ci.yml:52-55` nextest-clean invariant) | §5.2 | DONE | `grep -rn set_var src/ tests/` → one doc-comment mention only; ci.yml comment still true |

### §5 step 3 — test-file helpers

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 10 | `case_root(tag)` — unique tag, canonicalised | §3.2.1, §5.3 | DONE | 51 tags (41 `Tripwire` + 10 direct), **zero duplicates** (checked mechanically); all 41 equal their test fn name, satisfying "derive the tag from the test function name" |
| 11 | `write_fastq(path, id_prefix)`; mates share the prefix | §3.2.2, §5.3 | DONE | `:91`; every paired case stages both mates with one prefix |
| 12 | `write_fastq_gz` — real gzip bytes | §3.2.2, §5.3 | DONE | `:100`, `flate2::write::GzEncoder`; `flate2` is a normal dependency so no dev-dep was needed |
| 13 | `copy_fixture(name, dest)` from `CARGO_MANIFEST_DIR` | §5.3, A9 | DONE | `:61`, `:118`; both fixtures present in `test_files/` |
| 14 | `snapshot(root)` — recursive, **files only** | §5.3, A7 | DONE | `:133`; directories pushed, not recorded |
| 15 | `hash_guarded_inputs` | §5.3 | DEVIATED | Implemented as `hash_file` (`:152`) plus an inline map in `Tripwire::run` (`:259-271`). Same effect; the rename is not in §13's deviation list |
| 16 | `run_arm(root, args)` — `current_dir` + env var on the **child only** | §3.2.4, §5.3 | DONE | `:159`; `.env(DUMP_PLANNED_ENV, "1")` on the `Command` |
| 17 | `parse_planned(stderr) -> Vec<(String, PathBuf)>` | §5.3 | DONE | `:174`; `split_once('\t')`, panics on a malformed line |
| 18 | `compare(...)` with §4's signature | §4, §5.3 | DEVIATED | See **G1**. `guarded_inputs` dropped, `fastqc_capable` added, `created` is `&BTreeSet` not `&[PathBuf]`. Undocumented in §13 |
| 19 | `assert_tripwire(...)` driving T1/T3/T4/T5 | §5.3 | DEVIATED | Implemented as the `Tripwire` builder + `.run(stage)` (`:216-309`). Functionally equivalent and arguably better; the restructure is not in §13's deviation list |

### §5 step 4 — the 40 arm cases

All 40 rows exist as `#[test]`s **with the flags the table specifies**. Four rows differ in
positional count or path depth only (noted); one row's test name differs.

| # | Plan row | Test | Status | Notes |
|---|------|--------|--------|-------|
| 20 | 1 `se_trim_two_inputs` | same | DONE | site 1006 |
| 21 | 2 `se_trim_output_dir` | same | DONE | flags `-o out` as specified; inputs staged under `in/` rather than the root — strictly stronger (`-o` relocation from a different dir) |
| 22 | 3 `se_trim_fastqc` | same | DONE | exercises the allowance |
| 23 | 4 `se_trim_fastqc_args_outdir` | same | DONE | **`qc/` pre-created** at `:504` as §2.2(b) requires; comment records why |
| 24 | 5 `se_trim_demux` | same | DONE | real tabs, two 8 bp barcodes (`:128`); `bc.txt` in `guarded` |
| 25 | 6 `se_trim_demux_fastqc` | same | DONE | |
| 26 | 7 `se_trim_no_report_file` | same | DONE | |
| 27 | 8 `se_trim_gzipped_input` | same | DONE | real gzip via `write_fastq_gz` |
| 28 | 9 `se_trim_clumpify` | same | DONE | `--clumpify --cores 2`; **asserts `clumpify:` on stderr** (`:573`) so a memory-floor fallback cannot void the case |
| 29 | 10 `se_trim_from_ubam_input` | same | DONE | `ubam_test.bam` copied in |
| 30 | 11 `paired_trim_split_dirs` | same | DONE | **asserts `r2/` gained nothing** (`:598`) |
| 31 | 12 `paired_trim_multipair` | same | DONE | 2 pairs |
| 32 | 13 `paired_trim_output_dir` | same | DONE | `--paired -o out` over split dirs |
| 33 | 14 `paired_trim_retain_unpaired` | same | DONE | **short R2** via `write_fastq_short` so reads really strand |
| 34 | 15 `paired_trim_retain_unpaired_cores` | same | DONE | `--cores 2` → `parallel.rs` writer block |
| 35 | 16 `paired_trim_passthrough` | same | DONE | carrier in `guarded` |
| 36 | 17 `paired_trim_passthrough_cores` | same | DONE | |
| 37 | 18 `paired_trim_fastqc` | same | DONE | |
| 38 | 19 `paired_trim_passthrough_fastqc` | same | DONE | |
| 39 | 20 `paired_trim_basename` | same | DONE | `--basename fixed` |
| 40 | 21 `paired_trim_no_report_file` | same | DONE | |
| 41 | 22 `hardtrim5_fastq` | same | DONE | `fastqc_incapable`, site 466 |
| 42 | 23 `hardtrim5_ubam` | same | DONE | `+ --output-format ubam` |
| 43 | 24 `hardtrim3_fastq` | same | DONE | two inputs as the row specifies |
| 44 | 25 `hardtrim3_ubam` | same | DONE | one input (row's "two inputs" note is carried by #43) |
| 45 | 26 `clock_multipair` | same | DONE | 2 pairs, site 2670 |
| 46 | 27 `implicon` | `implicon_pair` | DONE | test **renamed**; `.fastq` inputs as the row specifies. Rename not in §13 |
| 47 | 28 `clump_only_se` | same | DONE | two inputs rather than one — additive |
| 48 | 29 `clump_only_paired_split_dirs` | same | DONE | site 2670 |
| 49 | 30 `clump_only_se_fastqc` | same | DONE | |
| 50 | 31 `clump_only_paired_no_report_file` | same | DONE | `clump_report_candidates` early return |
| 51 | 32 `clump_only_ubam_se` | same | DONE | site 809 |
| 52 | 33 `clump_only_ubam_paired_interleaved` | same | DONE | 1 BAM, site 706 |
| 53 | 34 `clump_only_ubam_paired_two_files` | same | DONE | 2 FASTQ, site 752 |
| 54 | 35 `ubam_out_se` | same | DONE | two inputs rather than one — additive |
| 55 | 36 `ubam_out_se_fastqc` | same | DONE | FastQC on a `.bam` primary |
| 56 | 37 `ubam_out_paired_two_files` | same | DONE | site 2047 |
| 57 | 38 `orphan_paired_fastq_from_interleaved_bam` | same | DONE | `Arm::Unplanned` + `fastqc_incapable` |
| 58 | 39 `orphan_paired_fastq_retain_unpaired` | same | DONE | legal FASTQ-out shape, records the arm's write set |
| 59 | 40 `orphan_paired_ubam_from_interleaved_bam` | same | DONE | `Arm::Unplanned`; correctly **not** marked incapable — see item 77 |
| 60 | *(extra)* `se_trim_fastqc_args_only` | — | DONE | §13 deviation 3; puts §3.3's contested `--fastqc_args`-alone evidence in the suite |

**Count check:** 40 plan rows → 40 tests, plus 1 extra = **41 arm cases**, matching §13.

### §5 step 5 — self-tests

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 61 | `comparison_reports_an_unplanned_path` | §5.5 | DONE | asserts the path, the site, **and** the "do NOT narrow" steer |
| 62 | `comparison_accepts_overplanning` | §5.5 | DONE | |
| 63 | `hook_is_inert_when_env_unset` | §5.5, V1 | DONE | merged per §13 deviation 2 into `hook_is_inert_when_env_unset_and_that_goes_red`; prefixed-line count `>0` set, `==0` unset, through the real binary |
| 64 | `end_to_end_goes_red_without_the_dump` | §5.5, V3 | DONE | same test; `compare(…, planned = &[], …)` → `Err` |
| 65 | `every_site_is_reached_by_some_case` — 11 distinct sites | §5.5, V5 | DONE | 11 shapes, no hard-coded line numbers, panics if two shapes share a site |
| 66 | Refused-run case — non-zero exit **and** `created.is_empty()` | §5.5, V6 | DONE | `a_refused_run_writes_nothing`; also asserts `"Output path collision"` |
| 67 | Four extra `compare` self-tests | §13 dev. 1 | DONE | empty run, missing dump, two sites, stray FastQC artefact — one per assertion inside `compare`, **except T3** (see G1) |

### §5 step 6 — module doc

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 68 | What the file guards | §5.6 | DONE | `:1-13` |
| 69 | The five prior instances | §5.6 | DONE | `:4-5` names #383, #388, #391, #409, #398 |
| 70 | Expectations come from the binary, not the author | §5.6 | DONE | `:10-13` |
| 71 | Weakening the allowance is the wrong fix | §5.6, §3.4 | DONE | `:15-21`, including #398's mirror-image lesson |
| 72 | The five arms where `--fastqc` is a no-op | §5.6, §2.2c | DONE | `:30-34`, links #421 |
| 73 | `collision_key`'s case-folding limit | §5.6, A8 | DONE | `:36-42` |

### §3 — assertions and mechanics

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 74 | **T1** `created \ allowed ⊆ planned`, keyed on `collision_key(root.join(p))` | §3.3 | DONE | `:412-435`; `root.join` handles absolute `-o` candidates |
| 75 | Converse (`planned ⊆ created`) **not** asserted | §3.3, A5 | DONE | absent by design; item 62 pins it |
| 76 | Allowance derived **once**, inside `compare`, from args as `--fastqc \|\| --fastqc_args` | §3.3 | DONE | `fastqc_active` (`:211`); single source of truth; disjunction preserved and now proven by item 60 |
| 77 | Allowance scoped to the FastQC-capable arms | §3.3 | DONE | Per-case `fastqc_capable` flag rather than a named `FASTQC_CAPABLE_ARMS` const. I re-derived the set from the code: `specialty.rs` has **0** `fastqc` refs; `run_paired_ubam_single_file` (`main.rs:1811-2004`) has none; `run_ubam_output_paired_single_file` (`:2414-2559`) **does** gate at `:2522`. So the incapable set is exactly {hardtrim5, hardtrim3, clock, implicon, orphan-14} — and the 8 `.fastqc_incapable()` tests cover exactly those five arms, with orphan-15 correctly left capable |
| 78 | Allowance self-policing: non-empty when active+capable, empty when inactive | §3.3 | DONE | `:338-359`; plus a third arm erroring if an `incapable` arm ever produces artefacts (i.e. #421 fixed) |
| 79 | **T3** no pre-existing file's bytes changed (all 15 arms) | §3.3, §3.2.3 | PARTIAL | Scope now matches §3.3's headline exactly — Reviewer A's F1 rebased the hash map onto `before` (`:258-274`), so every pre-existing file is covered, not only `guarded` (which is retained as the staging assertion). Remaining gap is **G1**: T3 lives in `Tripwire::run`, outside the pure `compare`, so it is the only assertion with no *permanent* self-test |
| 80 | **T4** `created` non-empty; `planned` non-empty on `Planned` | §3.3 | DONE | `:325-327`, `:373-379`; items 67's two tests pin both |
| 81 | **T5** every dumped line carries the same site | §3.3 | DONE | `:384-391`, implemented exactly as §3.3 words it (see item 99 for the A10 claim built on top of it) |
| 82 | `Arm::Unplanned` skips T1/T5 and asserts `planned` **empty** | §3.3 | DONE | `:361-371`; message tells the implementer to promote the case |
| 83 | §3.4 failure message: arm, dumped site, each unplanned path, full candidate list, correct fix, forbidden fix | §3.4 | DONE | `:427-434` + `panic!("{}: {msg}", self.tag)`; wording matches §3.4 including "Do NOT add them to the FastQC allowance" |
| 84 | Planned paths must resolve inside the case root | §3.3 (Rev. A) | DONE | `:409-419`; Reviewer A's F4 made it a path-boundary test (`root_key` + `MAIN_SEPARATOR`) rather than a raw string prefix |
| 85 | `guarded_inputs` = positionals + `--passthrough` + `--demux` barcode file, hashed | §3.2.3 | DONE | Checked case by case against `main.rs:74-83`: every positional is declared in all 41 cases; carrier declared in both passthrough cases; `bc.txt` in both demux cases. Staging is itself asserted (`:259-265`) |
| 107 | *(beyond plan)* The allowance may not swallow the whole write set, + its self-test | Rev. A F2/F3 | DONE | `:355-361` errors if `rest` is empty, closing the "filter that removes everything" hole that would leave T1 comparing two empty sets; `comparison_rejects_an_allowance_that_takes_everything` (`:1125`) proves it fires. Unreachable today, but it is a decay path the plan's self-policing argument did not cover |

### §9 — validations

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 86 | V1 hook fires and is inert unset | V1 | DONE | item 63; also confirmed by hand with the release binary |
| 87 | V2 whole chain can go red (one-off control) | V2 | DONE | §13 records it at `main.rs:918-924` naming all four report paths; I verified the range is exactly the `if !cli.no_report_file { … }` block with `planned.extend(candidates)` intact at `:925`, and that `git diff dev -- src/main.rs` is empty. I did **not** re-execute it (it requires editing source, out of scope here) — but Reviewer B's independent `MUTANT` probe reproduced the property live, and I saw the red myself |
| 88 | V3 same, permanently | V3 | DONE | items 63, 64, 61, 62 — no rebuild needed |
| 89 | V4 the allowance cannot silently widen | V4 | DONE | promoted to a permanent `compare` self-test (`comparison_rejects_stray_fastqc_artifacts`) rather than a binary run — §13 deviation 1, an improvement over the plan |
| 90 | V5 site coverage is real | V5 | DONE | `every_site_is_reached_by_some_case` → 11 distinct sites; passes. The 41 arm cases independently span all 11 sites |
| 91 | V6 the pre-flight's purpose | V6 | DONE | item 66; passes |
| 92 | V7 all cases green in **both** profiles | V7 | DONE | Pinned state: debug **51/51**, release **51/51**. Pre-F1–F4: 50/50 debug and 50/50 release × 12 runs (incl. clean-room target dir, serial, 32-thread), nextest 640/640 — no flakiness in any run |
| 93 | V8 no regression: full test, clippy, fmt | V8 | DONE | `cargo test` **641 passed / 0 failed** at the pinned state (640 pre-fix); clippy `-D warnings` clean (also run explicitly against the `integration_preflight_tripwire` target); fmt clean |
| 94 | V9 `cargo test --no-run` | V9 | DONE | clean; all test executables built |
| 95 | V10 orphan arms stay honest | V10 | DONE | 3 `orphan_*` cases pass with `planned` empty and T3 satisfied |

### §8 — assumptions that oblige the code

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 96 | A2 residual: probe `--fastqc_args "--svg"` once during implementation | A2 | MISSING | See **G4**. No record in §13, PROGRESS.md or `spikes/`; no case passes `--svg`. **I ran the probe** (release binary, 40-read fixture): write set is exactly `a_trimmed_fastqc.html` + `a_trimmed_fastqc.zip` — no `.svg` beside the primary, no directory. A2 holds; only the record is missing |
| 97 | A3 newline / non-UTF-8 documented as a limitation | A3 | DONE | `io.rs:97-98` on `dump_planned_lines` |
| 98 | A9 BAM fixtures copied, never read in place; both committed | A9 | DONE | `copy_fixture` used for all BAM cases; both files present in `test_files/` |
| 99 | A10 one pre-flight invocation per run, **machine-checked** by T5 | A10 | PARTIAL | See **G3**. T5 is implemented as specified, but it compares *distinct sites*; N invocations from one site emit byte-identical stderr (the format has no invocation boundary), so the one-per-run invariant is **not** machine-checked. A10's "Fixed, machine-checked" overstates what the specified assertion can deliver |
| 100 | A11 the `--basename` dependency recorded in the `Arm::Unplanned` doc **and** the follow-up issue | A11 | PARTIAL | Doc comment done and precise (`:191-198`). The follow-up issue is **not filed** — deliberate and documented (§13 "Not done, deliberately"; PROGRESS.md follow-up 3), pending Felix's call |

### §13 — the implementation notes' own claims

| # | Item | Source | Status | Notes |
|---|------|--------|--------|-------|
| 101 | "50 tests (41 arm cases + 9 self-tests)" | §13 | DONE | Verified exactly as claimed: 50 `#[test]`, 41 `Tripwire` cases, 9 self-tests. **Now 51** after Reviewer A's F3 |
| 102 | "Test total 589 → 640" | §13 | DONE | Verified **640**; 640 − 50 − 1 = 589, consistent with the stated baseline (V8's own caveat about counting `cfg`-gated blocks applies). **Now 641** after F3 |
| 103 | "68 lines, no candidate-assembly code touched" | §13, §7 | DONE | `git diff --stat dev` → `src/io.rs | 68 +`, one file. Candidate order therefore unchanged |
| 104 | V2 reverted; `git diff src/main.rs` empty | §13 | DONE | Empty now, and `grep -c MUTANT src/main.rs` → 0 |
| 105 | The four stated deviations are real | §13 | DONE | (1) six `comparison_*` tests present; (2) the two hook tests are one test; (3) `se_trim_fastqc_args_only` present; (4) `Outcome` has only `stderr` + `created` |
| 106 | "Nothing in the plan was skipped" | §13 | PARTIAL | See **G5**. Three real deviations are absent from the list: the `compare` signature (items 18/79), `hash_guarded_inputs` → `hash_file` + inline loop (15), `assert_tripwire` → the `Tripwire` builder (19). Plus the unrecorded A2 probe (96). Also cosmetic: §13 records V2's site as `src/main.rs:928`; the unmodified tree dumps `src/main.rs:927` (I confirmed by hand), consistent with the block having been edited at the time |

---

## Gaps (detail)

### G1 — `compare` lost `guarded_inputs`, so T3 sits outside the self-test net

**Expected** (§4): `compare(root, created, planned, guarded_inputs, args, arm) -> Result<(), String>`,
described as "Pure, so it can be unit-tested in both directions"; §3.3 lists T3 among the
assertions the comparison makes.

**Found:** `compare(root, created: &BTreeSet<PathBuf>, planned, args, arm, fastqc_capable: bool)`
(`tests/integration_preflight_tripwire.rs:313-320`). T3 was moved out, into
`Tripwire::run` (`:284-294`).

**Gap:** every other assertion has a *permanent* self-test proving it can go red — T1
(`comparison_reports_an_unplanned_path`), T4 twice, T5, the FastQC allowance, over-planning,
and now the allowance-takes-everything guard. T3 has none, and cannot get one while it lives
outside `compare`. Across all 51 tests the hash loop never observes a change, so if `hashes`
were empty, keyed on the wrong root, or `hash_file` were degenerate, T3 would pass everywhere
and nothing in the suite would notice. That matters more than usual because §3.3 says T3 "is
what gives the two orphan arms their coverage": on the three `orphan_*` cases T1 and T5 are
skipped by construction, so T3 is the substantive assertion.

T3 *does* fire today — Reviewer A planted a clobber-the-input mutation and confirmed it goes
red with the right message while T1 stays green (`CODE_REVIEW_A.md` §5). That is a one-off
manual check of the kind §5 step 5 exists to replace, so the gap is the missing permanent
proof, not a doubt about the assertion. `fastqc_capable` being added is a sound consequence
of §3.3's scoping and needs only documenting; the dropped parameter is the substantive part.

**To resolve:** either pass the before-hash map into `compare` and add
`comparison_detects_a_modified_input`, or extract the loop as a pure
`fn unchanged(before: &BTreeMap<PathBuf, u64>) -> Result<(), String>` and self-test that.
Then record the signature change in §13.

### G2 — T3's scope was `guarded`, not "every pre-existing file" — **RESOLVED during this audit**

**Expected:** §3.3 heads T3 "no pre-existing file's bytes changed (all 15 arms)"; §3.2
step 3 says "a content hash of every member of **`guarded_inputs`**" — the plan states the
assertion two ways.

**Found initially:** the narrower of the two — the hash map was built from `self.guarded`.
A no-op today (I checked every case: `guarded` lists exactly the staged files, and the one
non-input staging, `se_trim_fastqc_args_outdir`'s `qc/`, is a directory and so absent from
`before`), but a coincidence maintained by hand per case — the decay shape this file exists
to police.

**Resolution:** Reviewer A's F1 rebased the map onto `before` (`:258-274`) and kept `guarded`
as the staging assertion. T3 now matches §3.3's wording, and I re-verified all 51 tests green
in both profiles afterwards. Nothing further is required; §3.2 step 3's narrower wording
could be aligned with §3.3 for tidiness.

### G3 — A10 is not machine-checked

**Expected:** A10 — "Exactly one pre-flight invocation per run … `parse_planned` relies on
this; **T5 asserts it**", tabled as "Fixed, machine-checked". §3.3 adds: "The day someone
moves a pre-flight inside a loop, the parse would start unioning two arms' lists … T5 makes
that a visible failure."

**Found:** T5 is `planned.iter().find(|(s, _)| s != site)` (`:385`) — it detects two
*distinct* sites, which is exactly what §3.3's sentence "every dumped line carries the same
site" specifies, so the assertion itself is DONE. But a pre-flight moved inside its
per-pair loop emits N invocations from the **same** `file:line`; the dump format carries no
invocation boundary, so for a 2-pair run one call of 12 candidates and two calls of 6 are
byte-identical on stderr and `find` returns `None`.

**Gap:** the scenario A10 and T5's own comment claim to guard is undetectable. Practical
impact is limited — the union of per-pair lists equals the single list, so T1's strength is
unchanged for any path in its own pair's list — but the claim is stronger than the code,
inside the one file whose value rests on its claims being exact.

**To resolve:** either give the dump an invocation boundary (one
`TRIM_GALORE_PLANNED_BEGIN\t<site>` line per call, or an invocation counter in the site
field — ~4 lines in `io.rs`) and assert exactly one, or downgrade A10 and the T5 comment to
what is actually asserted: *no run reaches two distinct pre-flight sites*.

### G4 — A2's `--svg` residual was never recorded as probed

**Expected:** A2 — "Residual: `--fastqc_args \"--svg\"` (`fastqc.rs:81`) was not probed;
**probe it once during implementation**."

**Found:** no record in §13, PROGRESS.md or `spikes/`, and no case passes `--svg`.
`apply_fastqc_args` does accept `--svg` → `config.svg_output = true` (`fastqc.rs:81`), so
the residual was real.

**Gap:** documentation only. **I ran the probe during this audit**: `trim_galore
--fastqc_args "--svg" a.fastq` (release binary, 40 reads) exits 0 and writes exactly
`a_trimmed.fq`, both trimming reports, and `a_trimmed_fastqc.{html,zip}` — no `.svg` file
beside the primary and no directory. A2's two-suffix rule is therefore complete and the
allowance covers it.

**To resolve:** one line in §13 (or the module doc) recording the probe and its result.

### G5 — §13's deviation list is incomplete

**Expected:** §13 — "Nothing in the plan was skipped", followed by four deviations "all
additive".

**Found:** four further deviations, none listed: the `compare` signature (G1),
`hash_guarded_inputs` becoming `hash_file` plus an inline loop, `assert_tripwire` becoming
the `Tripwire` builder, and the `implicon` → `implicon_pair` rename. The middle two are
cosmetic — the builder is arguably the better shape — but a reader auditing §4 against the
code will not find what §4 describes.

**To resolve:** add them to §13's list. The `compare` signature entry should carry G1's
consequence, not just the shape change. If §13 is to stay the single record of what the
branch contains, it should also gain a line for Reviewer A's F1–F4, which are currently
recorded only in `CODE_REVIEW_A.md` §4.

---

## Test verification

Profiles: `cargo test` (debug), `cargo test --release`, `cargo nextest run`. All 51 tests in
the tripwire file plus the `io.rs` unit test pass in every profile run, at the pinned state.

| Test | File | Status |
|------|------|--------|
| `dump_planned_lines_renders_one_line_per_candidate` | `src/io.rs:977` | PASS |
| `se_trim_two_inputs` | `tests/integration_preflight_tripwire.rs` | PASS |
| `se_trim_output_dir` | same | PASS |
| `se_trim_fastqc` | same | PASS |
| `se_trim_fastqc_args_outdir` | same | PASS |
| `se_trim_fastqc_args_only` | same | PASS |
| `se_trim_demux` | same | PASS |
| `se_trim_demux_fastqc` | same | PASS |
| `se_trim_no_report_file` | same | PASS |
| `se_trim_gzipped_input` | same | PASS |
| `se_trim_clumpify` | same | PASS |
| `se_trim_from_ubam_input` | same | PASS |
| `paired_trim_split_dirs` | same | PASS |
| `paired_trim_multipair` | same | PASS |
| `paired_trim_output_dir` | same | PASS |
| `paired_trim_retain_unpaired` | same | PASS |
| `paired_trim_retain_unpaired_cores` | same | PASS |
| `paired_trim_passthrough` | same | PASS |
| `paired_trim_passthrough_cores` | same | PASS |
| `paired_trim_fastqc` | same | PASS |
| `paired_trim_passthrough_fastqc` | same | PASS |
| `paired_trim_basename` | same | PASS |
| `paired_trim_no_report_file` | same | PASS |
| `hardtrim5_fastq` | same | PASS |
| `hardtrim5_ubam` | same | PASS |
| `hardtrim3_fastq` | same | PASS |
| `hardtrim3_ubam` | same | PASS |
| `clock_multipair` | same | PASS |
| `implicon_pair` | same | PASS |
| `clump_only_paired_split_dirs` | same | PASS |
| `clump_only_paired_no_report_file` | same | PASS |
| `clump_only_se` | same | PASS |
| `clump_only_se_fastqc` | same | PASS |
| `clump_only_ubam_se` | same | PASS |
| `clump_only_ubam_paired_interleaved` | same | PASS |
| `clump_only_ubam_paired_two_files` | same | PASS |
| `ubam_out_se` | same | PASS |
| `ubam_out_se_fastqc` | same | PASS |
| `ubam_out_paired_two_files` | same | PASS |
| `orphan_paired_fastq_from_interleaved_bam` | same | PASS |
| `orphan_paired_fastq_retain_unpaired` | same | PASS |
| `orphan_paired_ubam_from_interleaved_bam` | same | PASS |
| `comparison_reports_an_unplanned_path` | same | PASS |
| `comparison_accepts_overplanning` | same | PASS |
| `comparison_rejects_an_empty_run` | same | PASS |
| `comparison_rejects_a_missing_dump` | same | PASS |
| `comparison_rejects_two_sites_in_one_run` | same | PASS |
| `comparison_rejects_stray_fastqc_artifacts` | same | PASS |
| `comparison_rejects_an_allowance_that_takes_everything` | same | PASS (added by Rev. A F3) |
| `hook_is_inert_when_env_unset_and_that_goes_red` | same | PASS |
| `every_site_is_reached_by_some_case` | same | PASS |
| `a_refused_run_writes_nothing` | same | PASS |
| `comparison_detects_a_modified_input` (T3) | — | **MISSING** (G1) |

Suite runtime: 10.6 s debug, 1.4 s release — consistent with §13's figures. Whole-repo
`cargo test` **641 passed / 0 failed**; `cargo nextest run` was 640/640 in 17.8 s at the
pre-F1–F4 state.

---

## Verdict

**INCOMPLETE — 4 items unresolved.** The plan's substance is implemented: all 40 case rows
exist with the flags specified (plus one extra), all 5 specified self-tests exist (plus
five), T1/T4/T5 and the FastQC allowance are implemented as written and each has a permanent
test proving it can fail, the hook is exactly §5 step 1's snippet with no `cfg` gate, the
module doc carries all six required contents, and V1/V3–V10 are all satisfied with V7 and V8
re-verified from scratch. Nothing found makes a current case pass vacuously, and the one
genuinely new writer that appeared in the tree during the audit was caught.

What remains, in the order an implementer should take it:

1. **G1 — restore T3 to the self-test net.** Pass the before-hash map into `compare` (§4's
   `guarded_inputs` slot) and add `comparison_detects_a_modified_input`, or extract a pure
   `unchanged(...)` and self-test that. T3 carries #409's data-loss harm and is the only
   substantive assertion on the three `orphan_*` cases, yet it is the only assertion with no
   permanent proof it can fire. Its scope is now correct (F1); only its location is not.
2. **G3 — make A10 true or make the claim true.** Either add an invocation boundary to the
   dump and assert exactly one, or reword A10 and T5's comment to "no run reaches two
   distinct pre-flight sites".
3. **G5 — complete §13's deviation list** with the `compare` signature (including G1's
   consequence), `hash_guarded_inputs`, `assert_tripwire`, the `implicon` → `implicon_pair`
   rename, and a pointer to Reviewer A's F1–F4.
4. **G4 — record the `--svg` probe.** Already performed here: write set is exactly the two
   allowed suffixes, so A2 holds. One line in §13 closes it; no code change.

**G2 was resolved mid-audit** by Reviewer A's F1 and needs no further work.

Two further items need a decision rather than work: A11's follow-up issue and §10's namer
issue are unfiled by explicit deliberation (§13, PROGRESS.md), pending Felix. #421 is filed.

One caveat on this audit's own evidence: another agent mutated `src/main.rs` mid-audit
(`CODE_REVIEW_A.md` §0). Every result reported here was produced against the clean tree —
`git diff --stat dev` = `src/io.rs` only, `grep -c MUTANT src/main.rs` = 0 — and re-run
after the revert. That mutation also served as an unplanned live negative control, and the
tripwire caught it on the 11 `se_trim_*` cases, which is the strongest single piece of
evidence in this report that the mechanism works.
