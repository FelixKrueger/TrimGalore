# Draft — doc sweep for the fastqc-rust 1.1.1 bump

**Status:** drafted, NOT applied. Needs an implementation trigger.
**Base:** `chore/fastqc-rust-1.1.1` @ `d5f41c5`
**Scope:** documentation and comments only. No behaviour change.

## Count correction

The review message said "nine places". That conflated two different staleness
classes. The accurate split:

- **Ten** sites assert `FastQC 0.12.1` (edits 1–10 below).
- **Three** sites list the `--fastqc_args` subset without `--png` (edits 11–13).
- **One** test comment describes `--svg` as still setting a flag (edit 14).
- **One** fixture note (`test_files/README.md:14`) describes an encoding
  heuristic that changed in 0.13.0 — see "Deferred" at the bottom.

`CHANGELOG.md:1398,1420` and `Docs/.../reference/changelog.md:654,676` also say
0.12.1 and are deliberately **not** touched: they are historical entries
describing what shipped at v2.1.0, and the docs changelog page is synced.

---

## Edits

### 1. `src/cli.rs:340` — `--fastqc` help text (user-visible)

```diff
-    /// Produces FastQC 0.12.1-compatible *_fastqc.html / *_fastqc.zip artifacts.
+    /// Produces FastQC 0.13.0-compatible *_fastqc.html / *_fastqc.zip artifacts.
```

### 2. `src/fastqc.rs:4-9` — module doc

The "same version we previously bundled in the Docker image" clause was true of
0.12.1 and is now false, so it goes rather than getting a new number.

```diff
 //! Replaces the v0.6.x-style shell-out to an external `fastqc` binary
-//! with an in-process call. Targets byte-equivalent output to Java
-//! FastQC 0.12.1 — the same version we previously bundled in the Docker
-//! image — see <https://github.com/ewels/FastQC-Rust> for the upstream
-//! crate. The bundled approach removes Java + the FastQC tarball as
-//! runtime dependencies, completing the "single static binary, zero
-//! deps" story for the v2.x rewrite.
+//! with an in-process call. Targets byte-equivalent output to Java
+//! FastQC 0.13.0 — see <https://github.com/ewels/FastQC-Rust> for the
+//! upstream crate. The bundled approach removes Java + the FastQC
+//! tarball as runtime dependencies, completing the "single static
+//! binary, zero deps" story for the v2.x rewrite.
```

### 3. `README.md:35`

```diff
-produces FastQC 0.12.1-compatible HTML + ZIP outputs without requiring Java
+produces FastQC 0.13.0-compatible HTML + ZIP outputs without requiring Java
```

### 4. `Docs/src/content/docs/quickstart.md:53`

```diff
-Outputs are FastQC 0.12.1-compatible HTML + ZIP files.
+Outputs are FastQC 0.13.0-compatible HTML + ZIP files.
```

### 5. `Docs/src/content/docs/guide/outputs.md:149`

```diff
-producing FastQC 0.12.1-compatible HTML + ZIP reports alongside the trimmed output.
+producing FastQC 0.13.0-compatible HTML + ZIP reports alongside the trimmed output.
```

### 6. `Docs/src/content/docs/guide/overview.md:10`

```diff
-library (FastQC 0.12.1-compatible HTML + ZIP, no Java required).
+library (FastQC 0.13.0-compatible HTML + ZIP, no Java required).
```

### 7. `Cargo.toml:36-42` — pin rationale

The "brand-new crate (v1.0.0 published 2026-04-26)" argument has expired, and
the 1.0.1-specific percentage-precision note is now two releases stale. Replaced
with the standing reason, trimmed to the house one-fact-per-line rule.

```diff
-# Bundled FastQC integration. Exact-pin: brand-new crate (v1.0.0 published
-# 2026-04-26) — treat upstream bumps as deliberate-test events. See
-# /Users/fkrueger/.claude/plans/bundled-fastqc-rust.md for the design notes.
-# v1.0.1 (2026-04-26) fixes percentage-precision in >>Overrepresented sequences
-# (ewels/FastQC-Rust#2); restored byte-identity to Java FastQC 0.12.1 for that
-# section.
+# Bundled FastQC integration, tracking Java FastQC 0.13.0. Exact-pinned: the
+# report is a tested output surface, so bumps are deliberate-test events.
+# See /Users/fkrueger/.claude/plans/bundled-fastqc-rust.md for the design notes.
 fastqc-rust = "=1.1.1"
```

### 8. `Cargo.toml:43-44` — noodles comment

```diff
 # uBAM input support (#316). Exact-pinned to the version already pulled
-# transitively via fastqc-rust 1.0.1 — using the `noodles` umbrella crate with
+# transitively via fastqc-rust 1.1.1 — using the `noodles` umbrella crate with
```

(Verified: 1.1.1 still resolves `noodles` to exactly 0.88.0, so the claim stays
true as well as current.)

### 9 + 10. `CLAUDE.md:67` and `AGENTS.md:67` — identical text, both need it

```diff
-`fastqc-rust` is exact-pinned (`=1.0.1`). The pin is deliberate: the crate is brand-new (v1.0.0 published 2026-04-26) and any bump should be taken as a deliberate-test event with output-byte-identity re-verified against Java FastQC 0.12.1. The CI `validation` job has a smoke test that asserts (a) `fastqc` is **not** on `$PATH` (so a green run proves the bundled library did the work) and (b) the produced `.zip` contains `summary.txt`, `fastqc_data.txt`, `Images/`, `Icons/`.
+`fastqc-rust` is exact-pinned (`=1.1.1`, tracking Java FastQC 0.13.0). The pin is deliberate: the bundled report is a tested output surface, so any bump is a deliberate-test event whose output should be re-verified against the Java FastQC release it tracks. The CI `validation` job has a smoke test that asserts (a) `fastqc` is **not** on `$PATH` (so a green run proves the bundled library did the work) and (b) the produced `.zip` contains `summary.txt`, `fastqc_data.txt`, `Images/`, `Icons/` — it does **not** check report content.
```

The trailing clause is added because the old text implied CI guards byte-identity.
It does not; it only checks four filenames are present.

### 11. `.github/dependabot.yml:50-56`

```diff
-    # fastqc-rust is exact-pinned (`=1.0.1`) because the bundled report must stay
-    # byte-identical to Java FastQC 0.12.1, which the validation job checks. The
-    # grouping above defeated that: PR #500 batched `=1.0.1` → `=1.1.0` into the
-    # weekly patch/minor PR, where it broke the build (`FastQCConfig` dropped
-    # `svg_output`, `threads` became `Option<usize>`) and would have retargeted
-    # the report at Java 0.13.0. That bump is a deliberate-test event with its own
-    # output re-verification, not a batched line item.
+    # fastqc-rust is exact-pinned (`=1.1.1`) because the bundled report is a tested
+    # output surface and a bump retargets it at a new Java FastQC release. The
+    # grouping above defeated that: PR #500 batched a bump into the weekly
+    # patch/minor PR, where it broke the build. Such a bump is a deliberate-test
+    # event with its own output re-verification, not a batched line item.
```

Drops the "which the validation job checks" overstatement and the now-historical
API-breakage detail (that is in `d5f41c5`'s commit message, where it belongs).

### 12. `src/cli.rs:348-351` — `--fastqc_args` help (user-visible)

```diff
     /// Additional arguments to pass to FastQC. Implies --fastqc.
     /// Common flags are translated to the bundled engine: --nogroup, --expgroup,
-    /// --quiet, --svg, --nano, --nofilter, --casava, -t/--threads, -o/--outdir.
+    /// --quiet, --png, --nano, --nofilter, --casava, -t/--threads, -o/--outdir.
+    /// --svg is accepted as a no-op; SVG graphs are the default.
     /// Unrecognised flags emit a warning and are ignored.
```

### 13. `src/fastqc.rs:87-89` — `Flags handled:` doc list

```diff
 /// Flags handled:
-///   `--nogroup`, `--expgroup`, `--quiet`, `--svg`, `--nano`,
-///   `--nofilter`, `--casava`, `-t`/`--threads`, `-o`/`--outdir`.
+///   `--nogroup`, `--expgroup`, `--quiet`, `--png`, `--svg` (no-op),
+///   `--nano`, `--nofilter`, `--casava`, `-t`/`--threads`, `-o`/`--outdir`.
```

### 14. `Docs/src/content/docs/guide/outputs.md:155`

```diff
-`--fastqc_args "..."` passes a subset of FastQC flags through. Currently supported: `--nogroup`, `--expgroup`, `--quiet`, `--svg`, `--nano`, `--nofilter`, `--casava`, `-t`/`--threads`, `-o`/`--outdir`. Other flags emit a warning and are ignored.
+`--fastqc_args "..."` passes a subset of FastQC flags through. Currently supported: `--nogroup`, `--expgroup`, `--quiet`, `--png`, `--nano`, `--nofilter`, `--casava`, `-t`/`--threads`, `-o`/`--outdir`. `--svg` is accepted as a no-op, since SVG graphs are the default. Other flags emit a warning and are ignored.
```

### 15. `tests/integration_preflight_tripwire.rs:29-32` — stale `--svg` claim

```diff
 //! Those two suffixes are the whole set: `--extract` is not in the curated `--fastqc_args`
-//! subset, so no `_fastqc/` directory appears, and `--svg` sets a rendering flag without
-//! writing an `.svg` file.
+//! subset, so no `_fastqc/` directory appears, and neither `--png` nor the no-op `--svg`
+//! writes an image file beside the primary artefacts.
```

---

## Proposed commit message

```
docs: retarget the FastQC references at 0.13.0 and document --png

The 1.1.1 pin tracks Java FastQC 0.13.0, but ten places still asserted
0.12.1 — including both --fastqc and --fastqc_args --help text and three
docs-site pages, so a GA would have shipped a binary naming the wrong
FastQC version. Three more places listed the --fastqc_args subset without
--png, and a test comment still described --svg as setting a rendering flag.

The Cargo.toml and dependabot.yml pin rationales argued from the crate
being brand-new, which stopped being true two releases ago; both now give
the standing reason instead. Both also claimed the validation job checks
report byte-identity against Java FastQC. It does not — it asserts four
filenames are present in the .zip — so that claim is gone and CLAUDE.md
and AGENTS.md say what the smoke test actually covers.

CHANGELOG.md and the synced changelog page keep their 0.12.1 references:
those describe what shipped at v2.1.0.

Documentation and comments only; no behaviour change.
```

---

## Deferred — needs a decision, not in this sweep

`test_files/README.md:14` explains why the phred64 fixture is not uniform-quality,
via a parenthetical describing FastQC's encoding heuristic as "upstream FastQC
0.12.1 behaviour: it keys off the lowest observed quality character".

**That heuristic changed in 0.13.0**, measured on this fixture:

| | Encoding reported | Per-base means |
|---|---|---|
| fastqc-rust 1.0.1 (FastQC 0.12.1) | `Illumina 1.5` | Q40 / Q2 |
| fastqc-rust 1.1.1 (FastQC 0.13.0) | `Sanger / Illumina 1.9` | Q71 / Q33 |

So the parenthetical's argument inverts under the new pin. How this line should
read depends on whether the `phred64` wiring below lands, so it is held back
rather than guessed at.
