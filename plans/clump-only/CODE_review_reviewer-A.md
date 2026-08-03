# Code review — `--clump_only` — Reviewer A

**Target:** `feature/clump-only` (working tree; no commits yet on branch)
**Scope:** `src/clump_only.rs` (new), `src/cli.rs`, `src/main.rs`, `src/io.rs`, `src/lib.rs`, `tests/integration_clump_only.rs` (new), `.github/workflows/ci.yml`, `docs/src/content/docs/modes/clump-only.md`, `CHANGELOG.md`.

## 1. Verdict

**APPROVE WITH REVISIONS.** One high-severity crash path on a specific input combination; several minor issues. The byte-identity contract, rejection matrix, and CI enforcement are otherwise sound.

## 2. Critical findings (ranked)

### C-1 (High) — Panic on `--paired --clump_only <single-uBAM>`

**Where:** `src/main.rs:333` → `src/main.rs:2076` (`run_specialty_paired`).

**What's wrong:** `Cli::validate()`'s paired-input check (`src/cli.rs:503-505`) intentionally allows `input.len() == 1` when `--paired` is set (to reserve the slot for interleaved uBAM). Line 193 of `main.rs` only bails for N=1 **when the input is NOT uBAM**, so a uBAM slips through. Then `--clump_only` reaches `run_specialty_paired`, which does:

```rust
for chunk in cli.input.chunks(2) {
    let (o1, o2) = output_names(&chunk[0], &chunk[1]);   // chunk[1] on len=1 → panic
```

`reject_ubam` inside `clump_only_paired` would have caught this, but the pre-flight fires first and indexes out of bounds. Verified reachable from user CLI. This is a real panic on user input.

**Fix:** hoist uBAM rejection to the top of the `if cli.clump_only { … }` block in `main.rs:326`, using the already-computed `any_bam` (line 167):

```rust
if cli.clump_only {
    if any_bam {
        anyhow::bail!(
            "uBAM input is not yet supported under --clump_only \
             (v1 is FASTQ in / FASTQ out only; see #353)."
        );
    }
    // …existing dispatch
}
```

This also removes the duplicated `detect_input_format` call inside `reject_ubam` (which currently re-detects for every input — see N-3).

## 3. Notable but non-blocking

### N-1 (Medium) — Misleading comment about `--no_report_file`

**Where:** `src/clump_only.rs:355`.

Comment claims "skipped if user passed `--no_report_file` at the caller", but the code unconditionally writes the report. `--no_report_file` is not threaded into `clump_only_single`/`_paired`. Either honor the flag (add a `write_report: bool` parameter, gated in `main.rs` on `!cli.no_report_file`) or delete the misleading comment. Given other specialty modes respect `--no_report_file` (`src/main.rs:890, 1209, 1429, 1647, 1787, 1902`), honoring it here is the consistent choice.

### N-2 (Low) — Silent-accept regression guard covers only `-q`, not `--stringency` or `-e`

**Where:** `tests/integration_clump_only.rs:331` (`silently_accepts_quality_flag`).

Plan resolution #9 silent-accepts three flags (`-q`, `--stringency`, `-e`) but only one is regression-guarded. If a future refactor accidentally routed `--stringency` into the trim path under `--clump_only`, no test would catch it. Add two more cases mirroring the existing one (same fixture, `--stringency 5` and `-e 0.3`, assert byte-identity).

### N-3 (Low) — Duplicate format detection

**Where:** `src/clump_only.rs:234` inside `reject_ubam`.

`detect_input_format` is called once in `main.rs:163` (result stored in `input_formats`) and again inside every `clump_only_single`/`_paired` call. Cheap (~24-byte peek), but not free. If C-1's fix is applied (`any_bam` gate in `main.rs`), `reject_ubam` becomes reachable only if the caller violated the invariant, so it can be downgraded to a `debug_assert!` or removed.

### N-4 (Low) — Duplicated case-fold normalization

**Where:** `src/main.rs:359` (`--clump_only` SE pre-flight) defines an inline `norm` closure that duplicates `io::norm_path` (`src/io.rs:34`, already `pub(crate)`). Trivial to swap.

### N-5 (Low) — Test coverage: sort-order tiebreaker cascade is only implicitly exercised

The deterministic tiebreak on `(minimizer, seq, qual, id)` is only exercised by inputs with distinct minimizers (`synth_records`). Neither unit nor integration tests construct records with identical minimizers-but-different-seq (which would exercise the `seq` tiebreak) or identical minimizer-and-seq-but-different-qual. `test_clump_only_no_adapter_detection` produces identical seq+qual+different-id, which does exercise the id-tiebreak — good — but the mid-cascade fields aren't covered. Non-blocking; the byte-identity multiset test already guarantees no records go missing.

### N-6 (Nit) — Report path uses full input basename (`sample.fq.gz_clumping_report.txt`)

Matches existing `report_name` convention (`src/io.rs:342`) so this is consistency, not a bug. Flagged only because the resulting filename has two extension-like segments; if a future cleanup consolidates report naming, both should move together.

## 4. What the implementation does well

The rejection matrix in `Cli::validate` (src/cli.rs:678-820) is exhaustive and each rejection carries a distinct, user-oriented message (not a generic "conflicting flags"). The SE and PE code paths are parallel and well-commented, with a single sort primitive (`sort_paired_by_key`) guaranteeing pair lockstep by construction rather than by convention. The empty-input case explicitly emits a valid empty gzip member so `zcat`-based downstream tools don't choke — that's the kind of failure mode planners often miss. The plus-line + CRLF normalization contract-scope note is honestly surfaced in `--help`, the docs page, and unit-test comments. CI enforces byte-identity via `paste - - - -`-then-`sort`-then-`diff` (the load-bearing shape that catches line-mix-up bugs a naive `zcat|sort|diff` would miss), and cross-run determinism via `md5sum`. Integration tests spawn the built binary rather than exercising library code — closest thing to a user acceptance test. Deviations from PLAN.md are documented in `PROGRESS.md` with rationale rather than silently dropped. The three documented deviations (cores>=1 vs >=2, empty-input gzip member, sort determinism verify-only) are all reasonable and don't hide problems.

**File written:** `/Users/fkrueger/Github/TrimGalore/plans/clump-only/CODE_review_reviewer-A.md`
