# Review: Release Plan for Trim Galore v2.0.0 (Oxidized Edition)

**Reviewer:** A  
**Date:** 2026-04-13  
**Plan file:** `plans/release-v2/plan.md`

---

## Strengths

1. **Staged rollout is the right strategy.** The plan correctly sequences confidence-building (author validation, beta testing) before broad distribution via bioconda. This is exactly how you ship a full rewrite of a tool that thousands of pipelines depend on. The explicit "bioconda is last, not first" principle is spot on.

2. **Byte-identical validation as a first-class requirement.** Anchoring the release on decompressed output matching the Perl original gives a concrete, falsifiable definition of "drop-in replacement." The existing CI already validates this for several modes, which de-risks the transition.

3. **Dependency elimination is a genuine win for the ecosystem.** Removing the Python/Cutadapt/pigz dependency chain from bioconda's solver graph is a material improvement. The `cutadapt ==5.2` pin has been a known pain point, and this eliminates it entirely.

4. **Decisions section resolves open questions clearly.** The plan doesn't leave ambiguity about SIMD tiers, Docker images, or Perl script archival. Each resolved decision has a rationale.

5. **Existing CI and release workflows are already in place.** The workflows at `.github/workflows/ci.yml` and `.github/workflows/release.yml` already cover Rust builds, cross-compilation, byte-identical validation, and binary packaging. This is not a greenfield CI effort.

6. **Comprehensive feature parity.** The SUMMARY.md and cli.rs confirm that every Perl CLI flag has a Rust equivalent. The plan's Phase 2 audit step is appropriate nonetheless.

---

## Concerns

### HIGH Severity

**H1. Bioconda recipe uses prebuilt binaries but bioconda strongly prefers building from source.**

Bioconda's policy is that recipes should build from source whenever possible. The proposed recipe downloads prebuilt tarballs from GitHub Releases. While some packages do this (e.g., tools distributed as static binaries), it is not the preferred path and may face pushback from bioconda maintainers. Bioconda has Rust build infrastructure (`{{ compiler('rust') }}`) and many Rust bioinformatics tools (e.g., `minimap2`, `noodles`) build from source in bioconda.

**Recommendation:** Prepare two recipe variants: (a) a source-build recipe using `{{ compiler('rust') }}` with `cargo build --release` in `build.sh`, and (b) the prebuilt-binary recipe as a fallback. The source-build recipe eliminates the musl-vs-glibc question entirely since bioconda's build infrastructure handles linking. Lead with option (a) in the PR and have (b) ready if build times or cross-compilation are problematic.

---

**H2. Release workflow trigger mismatch between plan and reality.**

The plan (Phase 3.2) says the release workflow is triggered by `workflow_dispatch` (manual trigger, "following RustQC pattern"). The actual workflow at `.github/workflows/release.yml` triggers on `push: tags: ["v*"]`. These are different mechanisms. The plan's description of validation steps (read version from Cargo.toml, verify tag doesn't exist) implies pre-release checks that a tag-push trigger does not inherently provide -- once you push the tag, the build starts immediately.

**Recommendation:** Either update the plan to match the existing tag-push trigger (and note that version validation happens before tagging), or update the workflow to `workflow_dispatch` if you genuinely want a manual gate. Tag-push is simpler and more common; just be clear about the sequence: bump version in Cargo.toml, commit, tag, push tag.

---

**H3. Bioconda recipe URL uses `-gnu` but release workflow builds `-musl` for Linux.**

The plan's bioconda recipe (Section 4.1) references tarballs named `*-x86_64-unknown-linux-gnu.tar.gz` and `*-aarch64-unknown-linux-gnu.tar.gz`. However, the actual release workflow builds `x86_64-unknown-linux-musl` and `aarch64-unknown-linux-musl` targets. The tarball names will contain `musl`, not `gnu`. This mismatch means the bioconda recipe would 404 on installation.

**Recommendation:** Fix the recipe URLs to use `musl` targets, matching the release workflow. Or, if you switch to building from source in bioconda (see H1), this becomes moot since bioconda's build environment uses its own glibc.

---

**H4. aarch64 Linux cross-compilation from x86_64 runner is fragile.**

The release workflow cross-compiles `aarch64-unknown-linux-musl` on `ubuntu-latest` (x86_64) using `gcc-aarch64-linux-gnu`. The plan says this should use `ubuntu-24.04-arm` (native ARM runner). The actual workflow does NOT do this. Cross-compiling Rust for aarch64-musl with the system cross-linker works for simple binaries but can fail with C dependencies. The `flate2` crate (via `zlib-rs`) includes C/assembly code that may not cross-compile cleanly.

**Recommendation:** Align the workflow with the plan: use `ubuntu-24.04-arm` for native aarch64-linux builds. GitHub Actions now offers ARM runners. This eliminates cross-compilation issues entirely.

---

### MEDIUM Severity

**M1. Version string format mismatch between plan and code.**

The plan (Section 1.3) says `--version` should report `Trim Galore - Oxidized Edition 2.0.0`. The actual code in `main.rs:25` prints `Trim Galore v{CARGO_PKG_VERSION}` (which would be `Trim Galore v2.0.0`). The clap `version` attribute derives from Cargo.toml, which would give just `2.0.0`. Neither matches the planned "Oxidized Edition" branding.

The bioconda test (Section 4.1) checks `trim_galore --help | grep "Oxidized Edition"` -- this will fail unless the help text actually contains that string. The clap `about` in cli.rs currently says "A fast, single-pass NGS adapter and quality trimmer" with no "Oxidized Edition" mention.

**Recommendation:** Update the clap `#[clap(version = ...)]` to include "Oxidized Edition" in the version string, and ensure the `about` text also mentions it. Verify the bioconda test will pass before submitting the recipe.

---

**M2. No rollback plan documented.**

The plan has a checklist item "Keep Perl v0.6.11 available via `conda install trim-galore==0.6.11` for rollback" but does not detail what happens if v2.0.0 breaks a real-world use case after bioconda release. Questions not addressed:
- Does bioconda allow "yanking" a version?
- How quickly can a revert recipe PR be merged?
- Should there be a v2.0.0 -> v0.6.11 downgrade path documented for users?
- What monitoring (GitHub issues template, bioconda feedstock issues) triggers a rollback decision?

**Recommendation:** Add a brief rollback section. At minimum: (a) the bioconda recipe for v0.6.11 should remain in the feedstock history (it will by default), (b) document `conda install trim-galore==0.6.11` in the README upgrade section, (c) define a severity threshold (e.g., "any byte-non-identical output on a common use case triggers immediate revert").

---

**M3. No `cargo publish` step in the release workflow.**

The plan (Phase 3.2, Step 4) describes publishing to crates.io as part of the release workflow. The actual `release.yml` workflow has no `cargo publish` step -- it builds binaries, creates a GitHub Release, and tests `cargo install` from the local path. The crates.io publication is not automated.

**Recommendation:** Either add a `cargo publish` job to the release workflow (after the GitHub Release job succeeds) or explicitly document that crates.io publication is a manual step. For a first release, manual `cargo publish` is perfectly fine.

---

**M4. CI workflow needs updating for repo restructure.**

The plan moves Cargo.toml to the repo root, but the existing CI workflow (`ci.yml`) uses `working-directory: optimus_prime` for Rust builds and caches `optimus_prime/target`. After restructuring, all these paths break. The release workflow has the same issue (`working-directory: optimus_prime`, `optimus_prime/target/...` paths).

**Recommendation:** Phase 1 (repo restructure) should explicitly include updating both `.github/workflows/ci.yml` and `.github/workflows/release.yml` to remove `optimus_prime/` path prefixes. This is an obvious dependency but not called out in the plan.

---

**M5. MSRV check mentioned but not implemented in CI.**

The plan (Section 3.1) says CI should include an "MSRV check (Rust 1.74)." The actual CI workflow does not do this -- it uses `dtolnay/rust-toolchain@stable`. An MSRV check requires an additional job that pins the Rust version to 1.74.

**Recommendation:** Add a CI job that uses `dtolnay/rust-toolchain@1.74.0` and runs `cargo check` (not a full test suite -- just compilation). This catches accidental use of newer Rust features.

---

**M6. `exclude` list in Cargo.toml may be incomplete.**

The proposed `exclude` list omits `optimus_prime/` (which won't exist after restructure, fine), but also omits `*.fastq.gz`, `*.fastq`, and `*.fq.gz` files in `test_files/`. Test data files should not be in the crates.io package. The `exclude = ["test_files/"]` entry handles the directory, but double-check that no test files live outside it.

**Recommendation:** Run `cargo package --list` after restructuring (the plan already includes this in Phase 5.1) and verify no large binary files sneak in. Consider adding `"*.gz"` to the exclude list as a safety net.

---

**M7. The CI validation job does not test RRBS, NextSeq, poly-A, or polyG modes.**

The existing byte-identical validation in CI covers: single-end Illumina, paired-end BS-seq, hardtrim5, clock, and demux. It does not validate RRBS (`--rrbs`/`--non_directional`), NextSeq 2-colour (`--nextseq`), poly-A (`--poly_a`), poly-G, `--retain_unpaired`, `--trim-n`, `--max_n`, `--clip_R1/R2`, or `--three_prime_clip_R1/R2`.

**Recommendation:** Expand the CI validation matrix to cover at least RRBS (directional + non-directional), NextSeq, and poly-A modes before the v2.0.0 release. Test files for these modes exist in `test_files/` (PolyA.fastq.gz, polyAT_R1/R2.fastq.gz, etc.). This is especially important because these are the modes most likely to have subtle behavioral differences.

---

### LOW Severity

**L1. The plan does not mention Windows.**

This is a Rust project that could trivially cross-compile for Windows. Many bioinformaticians on Windows use WSL, but native Windows support would be a differentiator. The plan doesn't mention it at all.

**Recommendation:** Not a blocker for v2.0.0, but note it as a future consideration. If the code uses no Unix-specific APIs (likely, given it's a file processor), a `x86_64-pc-windows-msvc` target could be added to the release matrix with minimal effort.

---

**L2. No mention of `man` page or shell completions.**

Clap can auto-generate man pages and shell completions (bash, zsh, fish). The Perl version likely had no man page either, so this isn't a regression, but it's a nice-to-have for a v2.0.0 release.

**Recommendation:** Consider adding `clap_mangen` and `clap_complete` as build dependencies for a future release. Not a v2.0.0 concern.

---

**L3. crate name availability is listed as a to-do but is a hard blocker.**

Phase 5.1 says "Verify `trim-galore` crate name is available." If someone else has claimed this name on crates.io, the entire crates.io publication strategy needs to change (e.g., `trim-galore-oxidized` or negotiate a name transfer).

**Recommendation:** Check crate name availability NOW, before investing effort in Cargo.toml metadata. This takes 10 seconds on crates.io.

---

**L4. `clippy` and `fmt --check` not in CI.**

The plan (Section 3.1) says CI should run `cargo clippy -- -D warnings` and `cargo fmt --check`. The actual CI workflow does neither. These are good practices and trivial to add.

**Recommendation:** Add these as a separate CI job (or steps in the existing Rust job).

---

**L5. The `run_exports` in the bioconda recipe may be unnecessary.**

The `run_exports` section pins downstream packages that depend on `trim-galore`. Since `trim-galore` is a standalone CLI tool (not a library), no other package depends on it at build time. `run_exports` is typically for shared libraries.

**Recommendation:** Remove the `run_exports` section from the bioconda recipe unless there's a specific reason for it.

---

**L6. Bioconda recipe `source` section uses multiple URLs without `fn` (filename).**

Bioconda's conda-build may have trouble with multiple conditional source URLs. The standard pattern for platform-specific prebuilt binaries in bioconda is to use `fn:` to set explicit filenames, and ensure only one source section is active per platform via selectors.

**Recommendation:** If going the prebuilt-binary route, validate the recipe locally with `conda-build` before submitting. Better yet, use a single source pointing to the GitHub tarball of the source code and build from source (see H1).

---

## Recommendations Summary

### Before proceeding (do now)

1. Check `trim-galore` crate name availability on crates.io (L3).
2. Decide on bioconda strategy: build-from-source vs. prebuilt-binary (H1). This affects Phases 3 and 4 significantly.
3. Fix the `-gnu` vs `-musl` URL mismatch in the recipe (H3), or make it moot by building from source.

### During Phase 1 (repo restructure)

4. Update both CI workflows to remove `optimus_prime/` path prefix (M4).
5. Update clap version/about strings to include "Oxidized Edition" (M1).
6. Align release workflow trigger with plan (tag-push vs. workflow_dispatch) (H2).

### During Phase 2 (validation)

7. Expand CI byte-identical validation to cover RRBS, NextSeq, poly-A modes (M7).
8. Add MSRV check, clippy, and fmt to CI (M5, L4).

### During Phase 3 (CI/CD)

9. Use native ARM runner for aarch64-linux builds (H4).
10. Decide whether `cargo publish` is automated or manual (M3).

### During Phase 5 (bioconda)

11. Add rollback plan with specific criteria and documented downgrade path (M2).
12. Remove `run_exports` from recipe (L5).
13. Validate recipe with `conda-build` locally before PR (L6).

---

## Overall Assessment

This is a well-structured release plan for a significant migration. The staged rollout is appropriately conservative, the feature parity is thoroughly documented, and the performance gains are compelling. The main risks are in the bioconda recipe strategy (prebuilt vs. source-build) and a handful of mismatches between the plan text and the actual workflow files. None of the concerns are showstoppers -- they are all addressable within the existing plan structure. The biggest single recommendation is to seriously consider building from source in bioconda rather than distributing prebuilt binaries, as this aligns with bioconda policy and eliminates the musl/glibc and URL-naming concerns simultaneously.
