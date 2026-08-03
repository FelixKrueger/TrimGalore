# Plan Review: Trim Galore v2.0.0 Release Plan

**Reviewer:** B
**Plan file:** `plans/release-v2/plan.md`
**Date:** 2026-04-13

---

## Strengths

1. **Staged rollout is the right call.** The five-stage progression (crates.io -> author validation -> beta testers -> GitHub release -> bioconda) puts bioconda last, which correctly reflects its blast radius. This is a mature, conservative release strategy for a tool with thousands of downstream users.

2. **Byte-identical validation is central.** Making output identity the acceptance criterion — rather than "it runs" — is the gold standard for a drop-in rewrite. The existing CI workflow (`ci.yml`) already validates this across multiple modes (single-end, paired-end, hardtrim, clock, demux), giving confidence the infrastructure exists.

3. **Dependency elimination is well-articulated.** The plan clearly communicates the value proposition to both users and bioconda maintainers: removing the `cutadapt ==5.2` pin and the Python dependency eliminates a major class of solver conflicts. This is a genuine contribution to the bioconda ecosystem.

4. **Decisions are resolved and documented.** The plan doesn't leave key decisions hanging — Perl archival, SIMD strategy, Docker images, MSRV, and nf-core module timing are all decided with clear rationale.

5. **CI/CD workflows already exist.** The release workflow (`release.yml`) and CI workflow (`ci.yml`) are already implemented and cover the four target platforms. The plan is largely documenting what's built rather than speculating.

6. **The `--cores` / nf-core module note is important.** Identifying that the nf-core module's `task.cpus - 4` heuristic is obsolete for the Rust version is a subtle but high-value observation. The decision to defer this is also correct — the current heuristic still works, just suboptimally.

---

## Concerns

### High Severity

#### H1. Bioconda recipe uses wrong tarball target names (musl vs gnu mismatch)

The bioconda recipe in Section 4.1 references tarballs named `*-x86_64-unknown-linux-gnu.tar.gz` and `*-aarch64-unknown-linux-gnu.tar.gz`, but the release workflow (both the plan in Section 3.2 and the actual `.github/workflows/release.yml`) builds `x86_64-unknown-linux-musl` and `aarch64-unknown-linux-musl` targets. The recipe URLs will 404.

**Fix:** Change the bioconda recipe URLs from `linux-gnu` to `linux-musl` to match what the release workflow actually produces.

#### H2. Linux aarch64 cross-compilation from x86_64 runner is fragile

The existing `release.yml` builds `aarch64-unknown-linux-musl` on `ubuntu-latest` (x86_64) via cross-compilation with `gcc-aarch64-linux-gnu`. However, the plan in Section 3.2 says this should use `ubuntu-24.04-arm` (native ARM runner). These contradict each other.

Cross-compiling Rust with musl for aarch64 is notoriously tricky — `musl-tools` only provides the x86_64 musl libc, and `gcc-aarch64-linux-gnu` links against glibc by default. The current workflow sets `CC_aarch64_unknown_linux_musl=aarch64-linux-gnu-gcc`, which will attempt to use the GNU linker with the musl target. This can produce binaries that are not truly statically linked or that fail at runtime on Alpine/musl-based systems.

**Recommendation:** The plan's suggestion of using a native ARM runner (`ubuntu-24.04-arm`) is the safer approach. Update the release workflow to match the plan, or use the `cross` tool for cross-compilation. If native ARM runners are available on GitHub Actions (they are, since late 2024), that is by far the simplest path. Verify that the resulting aarch64-musl binary runs on both glibc and musl-based Linux (e.g., test in an Alpine Docker container).

#### H3. No explicit rollback plan for bioconda

The plan mentions keeping v0.6.11 available via `conda install trim-galore==0.6.11`, but there is no documented rollback procedure if v2.0.0 is broken post-bioconda-merge. Questions that need answers:

- Can bioconda "yank" a version, or does the old recipe need to be re-published?
- If the recipe switches from `noarch: generic` to platform-specific, does `conda install trim-galore==0.6.11` still work? (It should — conda keeps historical versions — but this should be explicitly tested before the v2.0.0 recipe is submitted.)
- What's the timeline for the rollback decision? Who monitors, what signals trigger a rollback?

**Recommendation:** Add a rollback section to the plan with: (a) explicit verification that v0.6.11 remains installable after v2.0.0 lands, (b) a monitoring period (e.g., 2 weeks post-bioconda-merge), and (c) a named person responsible for watching GitHub issues and bioconda feedback.

#### H4. No `cargo test` integration tests validating byte-identity in the test suite itself

The CI validation job compares outputs using shell-level md5sum comparisons, which is good. But `cargo test` (run in the `optimus-prime` job) only runs unit tests inside the Rust crate. When the repo is restructured and the `optimus_prime/` prefix is removed, the validation job needs to be rewritten.

More importantly, the plan's Phase 2 validation is described as manual ("Felix runs it..."). There's no mention of encoding the byte-identity checks as automated integration tests that run on every PR. If a future change breaks byte-identity, the current CI would not catch it until the manual validation step.

**Recommendation:** Convert the validation CI job's md5sum checks into either: (a) `cargo test` integration tests that shell out to the built binary, or (b) a standalone CI job that runs on every PR (not just as a post-merge validation). This should be part of Phase 1 restructuring.

### Medium Severity

#### M1. License field in Cargo.toml is technically incorrect

The current `Cargo.toml` uses `license = "GPL-3.0"`, and the plan proposes `license = "GPL-3.0-or-later"`. The repo's actual `LICENSE` file is GPL v3 only (the standard FSF text). If the intent is "GPL-3.0 only" (not "or later"), the correct SPDX identifier is `GPL-3.0-only`. If the intent truly is "or later," that's fine, but it should be a deliberate choice, not a typo.

**Recommendation:** Clarify whether the license is `GPL-3.0-only` or `GPL-3.0-or-later` and ensure the Cargo.toml, LICENSE file, and bioconda recipe are all consistent. crates.io validates SPDX identifiers strictly.

#### M2. Version string in reports says "(Optimus Prime)" but plan says "Oxidized Edition"

The report header (`report.rs`, line 124) currently writes:
```
Trim Galore version: 0.1.0 (Optimus Prime)
```

The plan says `--version` should report `Trim Galore - Oxidized Edition 2.0.0`. The branding is inconsistent between the `--version` output, the report header, and the internal crate name. Additionally, the main.rs startup banner (line 25) prints `Trim Galore v0.1.0` with no edition suffix.

This matters for MultiQC compatibility — MultiQC parses the version line from trimming reports. If the format deviates from what MultiQC expects, parsing could break.

**Recommendation:** (a) Decide on a single branding string and apply it consistently across `--version`, the report header, and the startup banner. (b) Verify what version-line regex MultiQC uses to parse Trim Galore reports and ensure the new format matches. The MultiQC Trim Galore module source is at `multiqc/modules/trim_galore/trim_galore.py` — the regex should be checked before release.

#### M3. The `exclude` field in Cargo.toml should be under `[package]`

The plan's proposed Cargo.toml (Section 1.2) places `exclude = [...]` as a standalone top-level key outside `[package]`. In Cargo.toml, the `exclude` field must be nested under `[package]`:

```toml
[package]
name = "trim-galore"
# ...
exclude = ["test_files/", "plans/", ...]
```

As written in the plan, Cargo will either ignore the field or error during `cargo package`.

**Recommendation:** Move the `exclude` array into the `[package]` section in the plan's proposed Cargo.toml.

#### M4. Missing SHA256 checksums workflow step

The release workflow generates tarballs but the plan doesn't mention generating `.sha256` checksum files, even though Section 3.2 says "Each produces: `trim_galore-{version}-{target}.tar.gz` + `.sha256`". The existing `release.yml` doesn't include a `sha256sum` step. Without checksums, the bioconda recipe has no way to verify download integrity (the `sha256: PLACEHOLDER` values need to come from somewhere).

**Recommendation:** Add a `sha256sum $ARCHIVE.tar.gz > $ARCHIVE.tar.gz.sha256` step to the release workflow's packaging step, and upload the `.sha256` files as release assets alongside the tarballs. Alternatively, document that SHA256 values are computed manually after the release and pasted into the bioconda recipe.

#### M5. Crate name availability is not verified

The plan notes "Need to verify `trim-galore` is available on crates.io" as a parenthetical, but this is a blocking prerequisite for Phase 1. If the name is taken, the entire crates.io strategy needs revision (e.g., `trim-galore-rs`, `trimgalore`, etc.).

**Recommendation:** Check `https://crates.io/crates/trim-galore` now, before any other work begins. If the name is taken, resolve the naming before finalizing the plan. This is a 30-second check that could invalidate Phase 1.

#### M6. The `--gzip` / `--no_gzip` flag naming inconsistency

The Perl Trim Galore uses `--gzip` and `--dont_gzip` as flags. The Rust CLI only has `--dont_gzip` (verified in cli.rs). The plan's feature completeness audit (Section 2.1) should explicitly call out which Perl flags are accepted as no-ops and which are dropped. Users who have `--gzip` in their scripts will get an error if the flag is silently removed.

Additionally, the Perl version accepts `--no_report_file`, `--suppress_warn`, `--path_to_cutadapt`, `--cores` (which maps to Cutadapt's `-j`), and other flags that may not have direct equivalents. Any flag that existed in the Perl version but is absent from the Rust version should either be accepted as a no-op (with a deprecation warning) or explicitly documented as removed.

**Recommendation:** Generate a complete diff of `trim_galore --help` (Perl) vs the Rust CLI's accepted flags. Document every discrepancy in the plan. For flags like `--path_to_cutadapt` or `--gzip` that are no longer meaningful, accept them as no-ops with a warning.

### Low Severity

#### L1. No Windows target

The plan explicitly covers Linux and macOS but does not mention Windows. While most bioinformatics work happens on Linux/macOS, WSL2 usage is growing, and some bioinformatics classrooms use Windows. The Rust codebase appears to have no Windows-specific blockers (no Unix-only syscalls visible in the dependencies).

**Recommendation:** Consider adding `x86_64-pc-windows-msvc` as a future release target (not blocking v2.0.0, but worth noting as a post-release enhancement). At minimum, document that Windows users should use WSL2.

#### L2. No smoke test for `cargo install trim-galore` from crates.io

The release workflow tests `cargo install --path optimus_prime` (local path), not `cargo install trim-galore` from the actual crates.io registry. After the initial publish, there's no automated verification that the published crate installs correctly.

**Recommendation:** Add a post-publish CI step (or manual checklist item) that runs `cargo install trim-galore` from crates.io and verifies `trim_galore --version` outputs the expected string.

#### L3. The `run_exports` pin in the bioconda recipe may be unnecessary

The proposed recipe includes `run_exports: {{ pin_subpackage("trim-galore", max_pin="x") }}`. `run_exports` is primarily used when a package provides shared libraries that other packages link against. Since `trim_galore` is a standalone binary (not a library), `run_exports` is unnecessary and may confuse bioconda reviewers.

**Recommendation:** Remove the `run_exports` block from the bioconda recipe unless there is a specific reason for it (e.g., bioconda policy requires it for all packages).

#### L4. Test matrix gaps

The validation CI job covers single-end, paired-end, hardtrim5, clock, and demux modes. It does not cover:
- `--rrbs` / `--non_directional`
- `--nextseq`
- `--poly_a`
- `--poly_g`
- `--small_rna` / `--nextera` adapter presets
- `--retain_unpaired`
- `--clip_R1` / `--clip_R2` / `--three_prime_clip_R1` / `--three_prime_clip_R2`

These are all modes with distinct code paths. Phase 2.2 mentions these in the manual validation, but they should be in CI.

**Recommendation:** Expand the validation CI job to cover RRBS, nextseq, poly-A/poly-G, clipping, and adapter preset modes. The test files already exist (`BS-seq_10K_R1/R2`, `10K_150bp`, `polyAT_R1/R2`, etc.).

#### L5. `Cargo.lock` should be committed for binary crates

The plan mentions moving `Cargo.lock` to the repo root, which is correct. But it's worth explicitly noting that for binary crates (as opposed to library crates), Cargo.lock should be committed to version control. The Rust community convention is: commit Cargo.lock for binaries, gitignore it for libraries. The `.gitignore` update in Phase 1.1 should not exclude Cargo.lock.

**Recommendation:** Explicitly note in Phase 1.1 that Cargo.lock must be committed (not gitignored) since this is a binary crate.

#### L6. The `[profile.release]` with `strip = true` removes debug symbols

Setting `strip = true` in the release profile removes all debug symbols, which makes crash reports less useful. If a user reports a segfault or panic, the backtrace will contain only addresses, not function names.

**Recommendation:** Consider using `strip = "debuginfo"` instead of `strip = true`. This strips DWARF debug info (saves most of the size) but preserves symbol names for backtraces. Alternatively, keep `strip = true` but build and archive a separate debug-symbols file (`.dwp` or `.dSYM`) for each release, available for download if crash analysis is needed.

---

## Recommendations Summary

**Before starting Phase 1:**
1. Fix the musl/gnu mismatch in the bioconda recipe URLs (H1)
2. Check crates.io for `trim-galore` name availability (M5)
3. Clarify GPL-3.0-only vs GPL-3.0-or-later (M1)
4. Move `exclude` inside `[package]` in the Cargo.toml template (M3)

**During Phase 1:**
5. Resolve the aarch64 build strategy — native ARM runner vs cross-compile (H2)
6. Unify branding across `--version`, report headers, and startup banner (M2)
7. Verify MultiQC's version-line regex before changing the format (M2)
8. Generate a complete Perl vs Rust CLI flag diff and handle missing flags (M6)
9. Add SHA256 generation to the release workflow (M4)

**During Phase 2-3:**
10. Add a rollback procedure document for bioconda (H3)
11. Convert md5sum validation checks into automated integration tests (H4)
12. Expand CI test matrix to cover RRBS, nextseq, poly-A/G, clipping modes (L4)

**Post-release:**
13. Consider Windows target for future release (L1)
14. Add post-publish crates.io install verification (L2)

---

## Verdict

The plan is well-structured, technically sound in its overall approach, and reflects careful thought about the staged rollout strategy. The most critical issues are the musl/gnu URL mismatch in the bioconda recipe (a copy-paste bug that would cause the recipe to fail) and the aarch64 cross-compilation ambiguity between the plan and the actual workflow. Both are straightforward to fix.

The absence of an explicit rollback plan for bioconda and the lack of automated byte-identity integration tests are the most consequential gaps — they represent risks that compound under pressure (i.e., exactly when you most need them). Adding these before the bioconda submission would significantly reduce the risk of the release.

Overall assessment: **Ready to proceed with the fixes above.** No fundamental design issues; the concerns are about precision and completeness rather than direction.
