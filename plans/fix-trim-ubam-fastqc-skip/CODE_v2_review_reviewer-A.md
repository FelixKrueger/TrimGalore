# Code Review (round 2) — `fix(ubam): address code-review feedback`

**Reviewer:** A (independent; Reviewer B's round-2 report not read)
**Commit:** `76384cc` on `fix/trim-ubam-fastqc-skip` (on top of `72540f0`)
**Scope:** incorporation audit of 7 round-1 recommendations + regression hunt on the follow-up itself
**Round-1 report:** `plans/fix-trim-ubam-fastqc-skip/CODE_review_reviewer-A.md`

---

## Summary

**Verdict: APPROVE WITH REVISIONS.**

All 7 recommendations landed, and 6 of them landed cleanly. Two land *better* than
round-1 asked for: the `fastqc.rs` `output_dir=None` docstring was corrected beyond my
suggested text (my round-1 replacement still carried the wrong "current dir if unset"
claim — the implementer checked upstream and fixed it properly), and the third CI step
closes a genuine driver-coverage hole.

One **MEDIUM regression is newly introduced by this follow-up**: adding `set -o pipefail`
to the SE step turns `unzip -l … | grep -q …` (`.github/workflows/ci.yml:863`) into a
**~1 % flaky red CI step**. I reproduced it: **5 failures in 500 runs, all exit 141
(SIGPIPE)** — on a build where the assertion actually *passes*. `grep -q` exits on first
match, `unzip` dies on the next write, and `pipefail` now surfaces that as step failure.
This did not exist before `76384cc` because the job's `shell: bash -l {0}` had no
`pipefail`. Two candidate fixes validated at 500 runs each, both 0 failures.

That is a one-line change and the only thing standing between this and a clean APPROVE.

---

## Incorporation audit — 7 items

| # | Round-1 item | Landed | Verified how |
|---|---|---|---|
| 1 | `set -euo pipefail` on 3 new blocks (A-2 / B-M-1) | ✅ | `ci.yml:856, 871, 888`; job shell confirmed `bash -l {0}` at `:678-679` (no `-e`, no `pipefail`) — the gap was real |
| 2 | 3rd CI step for the FASTQ-pair-in driver (B-M-2) | ✅ | `ci.yml:883-896`; **ran the exact command locally** — `BS-seq_10K_R1_val.bam`, `BS-seq_10K_R1_val_fastqc.zip`, `ZIP_COUNT=1` all as asserted |
| 3 | `cargo test` regression guard (A-1 / B-M-3) | ✅ | `tests/integration_ubam_out.rs:537-599`; both new tests pass; full file 15/15 pass |
| 4 | `src/fastqc.rs` docstring (A-3 / B-L-1) | ✅ **+improved** | Verified against `fastqc-rust-1.0.1/src/runner.rs:224-233` — fallback really is parent-of-input, not cwd |
| 5 | `.bam` extension-dependency comment ×3 (A-4) | ✅ | `main.rs:1638-1639, 1774-1775, 1896-1897`; upstream dispatch confirmed at `fastqc-rust-1.0.1/src/sequence/bam.rs:424-441` |
| 6 | Narrow the `--cores` NOTE (A-6 / B-L-3) | ✅ | `main.rs:220-226`; claim "FastQC still uses it" is true — all 3 blocks pass `cli.cores` |
| 7 | README PE mate pooling (B-L-2) | ✅ | `README.md:133-135`; matches observed behaviour (one zip per pair, confirmed locally) |

**Deferrals — both correctly discharged:**

- **A-5 (upload-artifact glob)** — correctly dropped as a round-1 misread on my part.
  `/tmp/ubam_out_*` was already at `ci.yml:861`, and `git show dev:.github/workflows/ci.yml`
  confirms it predates the branch. It covers all three new `_fqc` dirs. No change needed.
- **A-OSS (`--phred64` + uBAM quality offset)** — filed as **issue #358**, open, with the
  correct reproducer and the `src/bam.rs::BamWriter::write_record` citation. Properly
  out of scope for this PR.

---

## Issues by area

### Logic

No defects in the Rust changes. The three `main.rs` blocks are comment-only. The
`fastqc.rs` change is doc-only. Behaviour is unchanged from `72540f0`, which round-1
already approved.

### Errors / CI robustness

One new issue — **M-1** below. Everything else in the shell hardening is sound.

**The `|| true` on `ZIP_COUNT` is correct and cannot mask a failure.** I checked this
specifically because `|| true` is normally a smell. Empirically, in an empty directory:

| variant | behaviour |
|---|---|
| with `\|\| true` (as committed) | prints `Expected 1 fastqc zip, got 0`, then `exit 1` |
| without `\|\| true` | `exit 1`, **no diagnostic** — aborts at the assignment |

The only non-zero path through that pipeline is zero matches, which the following
`[ "$ZIP_COUNT" = "1" ]` rejects anyway. So `|| true` buys the diagnostic and gives up
nothing. Precedence is right too: `||` binds looser than `|`, so it applies to the whole
pipeline, not just `tr`.

### Efficiency

Two new integration tests add ~2 s wall-clock (file total 3.6 s). The third CI step
trims a 10 K-pair BS-seq fixture plus a FastQC pass — a few seconds in a job that already
installs samtools and builds release. Fine.

### Structure

Consistent with existing conventions: `.status()` (10 prior uses in the file),
`fresh_tmpdir` with unique slugs (verified — all 25 slugs across `tests/` are distinct),
`assert!` with actionable messages. The step rename to
`Validate PE uBAM-out (interleaved BAM in) — …` correctly disambiguates against the new
sibling step.

---

## Fixes applied

**None.** Same reasoning as round 1: Reviewer B is running concurrently with no shared
state, so a same-file edit risks a clobber; and the project's phase separation makes this
a review turn, not an implementation turn. Exact replacement text is supplied below so
M-1 can be applied in one paste.

---

## Recommendations

### MEDIUM — M-1. `set -o pipefail` makes the SE step ~1 % flaky red (new regression)

`.github/workflows/ci.yml:863`:

```yaml
          unzip -l /tmp/ubam_out_se_fqc/ubam_test_trimmed_fastqc.zip | grep -q fastqc_data.txt
```

`grep -q` exits the instant it matches. `fastqc_data.txt` is entry 21 of a 24-line
listing, so `unzip` is usually still writing when `grep` goes away — it takes SIGPIPE and
exits 141. Before this commit that was invisible (no `pipefail`, and the pipeline's status
was `grep`'s 0). With `pipefail` the step now fails **on a healthy build**.

Measured on this machine against a real fixture zip:

```
500 runs of: set -euo pipefail; unzip -l <zip> | grep -q fastqc_data.txt
  495  exit 0
    5  exit 141   (SIGPIPE)
```

This is the worst failure shape available: an intermittent red that points at a
nonexistent bug and burns whoever re-runs the job. ~1 % per run on this step.

Both candidate fixes validated at **500 runs, 0 failures**:

```bash
# Fix A — grep -c reads to EOF, so it structurally cannot SIGPIPE the writer.
unzip -l /tmp/ubam_out_se_fqc/ubam_test_trimmed_fastqc.zip | grep -c fastqc_data.txt >/dev/null
```

```bash
# Fix B — no pipe at all; the listing also lands in the failure artifact
# (it is under /tmp/ubam_out_*, already globbed at ci.yml:861).
unzip -l /tmp/ubam_out_se_fqc/ubam_test_trimmed_fastqc.zip > /tmp/ubam_out_se_fqc/ziplist.txt
grep -q fastqc_data.txt /tmp/ubam_out_se_fqc/ziplist.txt
```

**I recommend Fix B** — it keeps `grep -q`'s readable intent, is immune by construction,
and makes a red run diagnosable from the uploaded artifact. Negative control checked on
Fix A as well: `grep -c NO_SUCH_ENTRY >/dev/null` still exits 1, so gating is preserved
either way.

Keep `set -euo pipefail`. The hardening is right; only this one pipeline needs to stop
early-exiting.

**Forward guidance:** `ci.yml:487-490` in the `validation` job has four more
`unzip -l … | grep -q …` pipelines. They are safe *today* only because that job also
runs `bash -l {0}` without `pipefail`. If the job-wide hardening I floated in round 1
ever happens, those four inherit this exact flake — fix them in the same pass.

### LOW — L-2. New CI steps don't clear their output dir before counting

`mkdir -p` is idempotent, so a pre-existing `*_fastqc.zip` in `/tmp/ubam_out_pe_fqc`
would break the `ZIP_COUNT = 1` assertion. Harmless on a fresh GHA runner, and the Rust
tests are already stricter (`fresh_tmpdir` wipes first). Worth a `rm -rf <dir>` before
`mkdir -p` only if these steps ever get reused or run locally via `act`.

### LOW — L-3. `pipefail` hardening stops at the 3 new steps

The other seven `run:` blocks in `validation-ubam` remain non-gating, so the job is now
internally inconsistent — a reader can't tell whether the absence of `set -euo pipefail`
in a block is deliberate. Still out of scope for this PR (round 1 said the same), but it
is now a visible seam. Note M-1's guidance before doing that sweep.

### NIT — N-1. `run()` doc says "the input file's parent" for a param named `output_path`

`src/fastqc.rs:33-34`. Factually correct in context (fastqc-rust's fallback is the parent
of the file *handed to it*, i.e. the trimmed output), and in practice it resolves to the
same directory as trim_galore's input parent when `--output_dir` is unset. But
"input file" sitting three lines under `output_path` reads oddly. Suggest "…the analysed
file's parent directory". Purely cosmetic.

### NIT — N-2. PE assertions skip the HTML

SE checks `.zip` **and** `.html` (CI `:861-862`, test `:562-565`); both PE paths check
only `.zip`. Since HTML and zip are written by the same upstream call, this is only an
asymmetry, not a hole.

### NIT — N-3. Commit body says the phred64 issue "will be filed"

It has been — #358 — but the pushed commit body has no back-reference. Not worth an
amend; the CHANGELOG and PR thread carry the context.

---

## Verification performed

Everything below was executed, not inferred.

**Quality gates (both green):**

```
cargo fmt --all -- --check                              → FMT OK
cargo clippy --all-targets --release -- -D warnings      → clean
cargo test --test integration_ubam_out                  → 15 passed; 0 failed
```

**The two new tests pass and exercise the real path** — stderr confirms the call actually
happens rather than the assertions passing on stale files:

```
Running FastQC on /tmp/claude-501/tg_int_ubam_out_se_fastqc/ubam_test_trimmed.bam
Running FastQC on /tmp/claude-501/tg_int_ubam_out_pe_fastqc/ubam_paired_test_val.bam
test ubam_out_se_fastqc_produces_report ... ok
test ubam_out_pe_fastqc_produces_exactly_one_report ... ok
```

**New CI step #3 reproduced locally** with the release binary and the exact command line
from `ci.yml:889-895`:

```
BS-seq_10K_R1_val.bam          ✅ asserted path correct
BS-seq_10K_R1_val_fastqc.zip   ✅ asserted path correct
BS-seq_10K_R1_val_fastqc.html
ZIP_COUNT=1                    ✅ exactly-one guard holds
```

Naming traced independently through `src/io.rs:125-141` (`paired_bam_output_name` →
`strip_fastq_extensions(input_r1)` + `_val.bam`) and `fastqc-rust`'s
`strip_extensions` (`runner.rs:279-290`, which strips `.bam`) — so the CI assertion is
derived from the code, not just from one observed run.

**Upstream claims audited against `fastqc-rust-1.0.1` source:**

- Extension dispatch — `sequence/bam.rs:424-441`: `name_lower.ends_with(".bam")` →
  `BAMFile`, else fall through to `FastQFile`. The docstring and the three `main.rs`
  comments are accurate, including the "not content" emphasis.
- `output_dir=None` fallback — `runner.rs:224-233`: `config.output_dir` if set, else
  `group.files.first().parent()`, else `.`. Confirms the corrected docstring and refutes
  the pre-existing "current directory" claim.
- Related check: `apply_fastqc_args` (`src/fastqc.rs:73-106`) has no `-f`/`--format`
  mapping, so a user **cannot** override the format and make a `.bam` be parsed as FASTQ.
  The extension-dispatch dependency is therefore fully described by the new comments —
  no hidden second trigger.

**Shell semantics:** `validation-ubam` sets `shell: bash -l {0}` (`ci.yml:672-679`), which
GitHub uses verbatim — no injected `-eo pipefail`. Round-1 M-1's premise confirmed. YAML
re-parsed after the edit (`ruby -ryaml`): valid, all three fastqc steps present with the
renamed second step.

**`|| true` and SIGPIPE:** measured, results in the sections above.

**Call-site count:** `grep -rn "fastqc::run" src/` → 7 sites, all in `main.rs`. The
docstring's "Callers in `src/main.rs`" is accurate (round 1's "two callers" was the stale
claim).

---

## Conclusion

The follow-up is faithful to round-1 intent and, in two places, better than what was
asked for. The Rust side is comment/doc-only and verified against upstream source; the
test additions are portable, correctly-scoped, and follow the file's conventions; the CI
coverage gap is genuinely closed.

Ship it after the one-line M-1 change. `set -o pipefail` was the right call, but it needs
the `unzip … | grep -q` pipeline to stop early-exiting or the step will go red about once
per hundred runs on perfectly good builds — which would discredit exactly the regression
guard this PR exists to install.

**Verdict: APPROVE WITH REVISIONS** (M-1 before merge; L/NIT items optional).
