# Plan: release v2.4.0

**Date:** 2026-11-02
**From:** `dev` @ `d04f06b`, 103 commits and 56 CHANGELOG entries ahead of `master` (v2.3.0, 7 July)
**Revised:** 2026-11-02, after #479 landed and #1988/#1989 clarified the nf-core side

## Version

**v2.4.0**, a minor bump. The Unreleased section carries behaviour changes, not just fixes:
`--clumpify` bin sizing (output bytes move again), the poly-G/poly-A stage order, `--library`
presets replacing four clip flags, paired runs writing every output to one directory, uBAM output
dropping FASTQ header descriptions, several new refusals that previously succeeded, and a `gzp`
floor raised to 2.0.3. None break a documented interface, so not 3.0.0; too much behaviour change
for 2.3.1.

## Decisions needed before step 1

1. **Confirm v2.4.0.**
2. **Ship with `Security audit` red?** The three `zip` advisories (#501–#503) are unreachable from
   the binary — extraction is gated on `config.do_unzip`, which `src/fastqc.rs` never sets, and
   `--extract` is not in the `--fastqc_args` allowlist. No fixed version is reachable from our
   lockfile: `fastqc-rust 1.0.1` requires `zip ^2`, the fix is `>= 9.0.0`, and the bump is waiting
   upstream in ewels/FastQC-Rust#14. Recommendation: **ship**. The release workflow does not gate on
   audit, holding 102 commits for a third-party crate bump costs more than it buys, and the
   reasoning is already recorded on all three issues.
3. **Which open PRs land first** (section below).
4. **#479** — rebase and merge, or close. See below.

## Open PRs

| PR | Base | Disposition |
|---|---|---|
| #509 wiki + template 0.10 | dev | Merge before the FF, so the wiki ships with the release |
| #507 gh-actions patch/minor | dev | Merge |
| #500, #486 | dev | See below — both now need a rebase onto `d04f06b` |
| #500 cargo patch/minor (6) | dev | Merge — then re-run `cargo audit` locally, in case a bump clears something |
| #486 npm patch/minor (7) | dev | Merge |
| #508 satori 0.26.0 → 0.33.5 | master | **Seven minors on a docs dep. Review, do not wave through.** Retarget to dev |
| #506, #499, #496, #491, #490, #489, #487, #484, #483 | master | Security bumps that target the default branch. Most resolve or close themselves once `master` catches up. Re-triage *after* the FF, not before |
| ~~#479~~ gzp `mem::forget` | dev | ✅ **Done.** Rebased, comments trimmed, rebase-merged as `d04f06b` via #526 with @sujeito-operator's authorship preserved. Raises the `gzp` floor to 2.0.3 |

## Version-bearing files

| File | Change |
|---|---|
| `Cargo.toml` + `Cargo.lock` | `2.3.0` → `2.4.0` |
| `CITATION.cff` | `version`, and `date-released` to the release date |
| `README.md:90` | `:latest` "currently `v2.3.0`" and the `:v2.3.0` example |
| `Docs/.../install.md:63,64,85` | tag table rows and the `--version` output sample |
| `Docs/.../performance/threading.md:65` | "**The next release uses roughly twice this on the plain path.**" — this *is* that release. Rewrite to present tense. CHANGELOG already has an entry anticipating it |
| `Docs/.../reference/credits.md:35` | Add the v2.4.0 Zenodo DOI — **after** Zenodo mints it, so post-release |
| `Docs/.../reference/changelog.md` | Synced page. Do not hand-edit |
| `.github/workflows/release.yml:13` | Historical comment about v2.3.0. Leave |

## Sequence

1. Land the dev PRs above. Confirm `ci.yml` green on `dev` afterwards (audit excepted).
2. On `dev`: bump `Cargo.toml`, run `cargo build` to move `Cargo.lock`, update `CITATION.cff`,
   `README.md`, `install.md`, `threading.md`. Convert `### Unreleased` to
   `### Version 2.4.0 (Released on <date>)` and open a fresh empty `### Unreleased`.
3. Commit, push `dev`, let CI run. **Do not tag.**
4. Fast-forward `master` from `dev` and push `master`.
5. Dispatch `release.yml` with `--ref master` and **`dry_run: true`** first. The workflow supports
   it and suppresses every destructive effect. Cheap insurance given the v2.3.0 partial ship.
6. Dispatch for real with `--ref master`.
7. **Verify all four channels** — `gh run watch` returning is not success:
   tag, GitHub Release, crates.io, Docker `:latest`. Check `conclusion`, not just completion.
8. Post-release: add the Zenodo DOI to `credits.md`; re-triage the nine master-targeted dependabot
   PRs; delete the merged `wiki/init`.

## Hazards, from the v2.3.0 incident

- **Never pre-tag** between the version push and the dispatch. The workflow creates the tag. A tag
  at a different commit is a hard error; at the same commit it is tolerated but not the intent.
- **Push-only silently skips the release jobs.** The dispatch is required.
- **`gh run watch` returning ≠ succeeded.** Verify each channel independently.
- The release does not gate on lint or audit, so `dev` must be green *before* the FF, not after.

## Open questions

- **Q1.** Does the `fastqc-rust` pin move in this release? Only if Phil merges ewels/FastQC-Rust#14
  and cuts a release first. Recommendation: do not block; the pin moves in a follow-up, which is
  also the deliberate-test event CLAUDE.md asks for.
- **Q2.** Is `master` a clean fast-forward from `dev`? `master`'s HEAD is `c6528e5`, a merge of
  master into dev, which suggests yes — verify with `git merge-base --is-ancestor` before step 4.

---

## Revision, 2026-11-02

### Landed since the plan was written

`d04f06b` (#479 via #526) — drops the local `mem::forget` in `Sink::finish`'s `--cores N` arm and
raises the `gzp` floor to 2.0.3. Rebase-merged so @sujeito-operator keeps authorship. That makes
**104 commits / 56 entries** and changes `Cargo.toml`, so the version bump in step 2 now sits on
top of a dependency change rather than a clean tree.

**Landed 2026-11-03:** #509 (rebase-merged, so the build and update stay two commits), #507 and
#486 (squashed). `dev` is now `da8bbd0`. All three were audit-red only.

**#500 (cargo-patch-minor) is held, not merged.** Dependabot batched `fastqc-rust` `=1.0.1` →
`=1.1.0` into the weekly patch/minor group. It fails *every* job, not just the audit: `FastQCConfig`
dropped `svg_output` and `threads` became `Option<usize>`, so `src/fastqc.rs` fails with `E0308`,
`E0609` and `E0277`. Commented asking for the bump to be split out.

### Two work items to fold into the release

**1. Make a repeated boolean flag a no-op. Decided: booleans only.**

`--poly_a --poly_a` and `--nextseq 20 --nextseq 20` both hard-error with
`cannot be used multiple times`. Perl v0.6.x did not — verified against `Getopt::Long`, where a
repeated boolean simply sets again and a repeated value flag is last-wins. So this is a v2
regression, and it is what turned an nf-core config bug into a dead pipeline for a user
(nf-core/rnaseq#1989).

**Scope: booleans only.** A repeated boolean carries no information to lose, so accepting it is
free. A repeated value flag with two *different* values is genuinely ambiguous and the refusal is
worth keeping — someone passing `--length 20 --length 50` has a config bug they should hear about.

**Mechanism, probed in-tree and reverted:** add `overrides_with = "<field>"` pointing at itself on
each of the 28 `bool` fields in `src/cli.rs`. Verified on `poly_a`:

| Case | Result |
|---|---|
| `--poly_a --poly_a` | accepted, flag set (1,229 reads trimmed) |
| `--polyA --poly_a` (alias + primary) | accepted, flag set — aliases share the arg id |
| `--nextseq 20 --nextseq 20` | still refused, as intended |

Rejected alternative: `Command::args_override_self(true)`, which exists in clap 4 and would be a
one-liner, but applies to every argument including value flags — the scope we just ruled out.

28 mechanical attribute edits, no behaviour change beyond idempotence. Needs a test per class:
one repeated boolean accepted, one repeated value flag still refused, so a later `args_override_self`
cannot silently widen it.

**2. Document the non-restartable-input refusal (#379).** The fix landed in
`format::non_restartable_kind`, with the refusal at `cli.rs` and an integration test, but it
appears on **zero** docs-site pages — while being on four wiki pages. A user hitting
`trim_galore <(zcat reads.fq.gz)` gets a correct error and no documentation behind it. A few lines
on the input page: the input must be a regular file, why, and the `zcat … > tmp.fq` workaround.

### Also worth a line in the release notes

trim_galore 2.1.0 shipped in nf-core/rnaseq 3.26.0, and `--poly_g` was broken for Read 2 until
v2.3.0 (#321) — it looked for poly-C at R2's 5' end rather than poly-G at the 3'. Anyone still
pinned to 2.1.x on paired data has mis-trimmed R2. 3.27.0 ships 2.3.0. The v2.4.0 notes should say
this plainly rather than leaving it in a closed issue, because the failure is silent.

### fastqc-rust 1.1.0 — deliberately not in this release

Checked after #500 failed. 1.1.0 (9 Oct) is a Java-version tracking release, not a maintenance bump:

- **It does not fix the advisories.** Still `zip ^2`, so the pin move buys nothing on #501–#503.
  ewels/FastQC-Rust#14 remains the only route.
- **It tracks Java FastQC v0.13.0**, where `=1.0.1` tracks 0.12.1. Phred+33 becomes the default and
  encoding autodetection is gone; Basic Statistics gains Mean/Median Length and drops "Sequences
  flagged as poor quality"; `--min_length` filters reads rather than padding plots; SVG replaces PNG
  as the default graph format. The validation job's smoke test checks that the `.zip` *contains*
  `summary.txt`, `fastqc_data.txt`, `Images/` and `Icons/` — it would pass while the contents moved.
- **The API break is wider than `svg_output`:** `threads` → `Option<usize>`, new `max_length`,
  `phred64`, `decompress_threads`, `png_output`, and `SequenceFileGroup::new` now takes file openers.
- `--threads` is now a whole-run budget with parallel gzip via rapidgzip, and `--template` no longer
  uses `-t`. Our `--fastqc_args` allowlist forwards `-t`/`--threads` straight through, so what we
  forward has changed meaning underneath us.
- One upside: the `native-zlib` feature is gone, which is the fix for ewels/FastQC-Rust#3 — the
  default features that silently overrode our `zlib-rs`-only flate2 config. Binaries are pure Rust
  and static now.

**Decision: stay on `=1.0.1` for v2.4.0.** The migration needs output re-verification against
0.12.1 → 0.13.0, an API port in `src/fastqc.rs`, and a review of what `--fastqc_args` forwards. That
is its own release, and it does not solve the problem we would reach for it to solve.

### Dependabot now enforces the exact pins

`chore/dependabot-ignore-exact-pins` adds `fastqc-rust` and `mimalloc` to the cargo `ignore` list,
beside the existing `noodles` entry. All three are exact-pinned in `Cargo.toml` with a stated reason,
and the patch/minor grouping was overriding that — #500 is what it looks like when it does.
`mimalloc` had the same exposure and is listed before it bites, since allocator changes move the peak
RSS that the clumpify-memory job gates on.

### Still open from the original plan

Q1 (does the `fastqc-rust` pin move?) — unchanged, and now clearer: ewels/FastQC-Rust#14 is
rebased, green and waiting on Phil. Do not block the release on it.

The three `zip` advisories (#501–#503) now carry a comment each explaining why they are
unreachable, so shipping with `Security audit` red is documented rather than silent.
