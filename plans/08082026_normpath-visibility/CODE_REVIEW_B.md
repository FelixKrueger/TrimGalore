# Code Review B — #399 `norm_path` visibility demotion + doc re-homing

**Reviewer:** B (independent; A reviewing in parallel)
**Branch:** `fix/399-normpath-private` @ `315557d` — base `dev` @ `ac08a97`
**Diff:** `git diff ac08a97..315557d` — one file (`src/io.rs`), 3 hunks, +10/−12
**Plan:** `plans/08082026_normpath-visibility/PLAN.md` (r2)

## Verdict

**Approve.** No Critical, no High. The demotion is correct, complete, and verified by
real compiles including a positive control that proves the safety net fires. The doc
re-homing is a genuine improvement, not a shuffle: the trade-off paragraph now sits on
the item that actually makes the trade-off, on the page rendered docs readers land on,
and two existing pointers in `main.rs` that already said "see `io::collision_key`" now
resolve to it. One **Medium** on the accuracy of the relocated prose (pre-existing
wording, but promoted to a public API page by this change, so this is the cheap moment
to fix it) and three **Low** prose nits.

---

## ⚠️ Tree hazard discovered mid-review — affects anyone re-running these numbers

**The shared working tree moved during this review.** It is now on branch
`fix/400-clump-help-text` @ `0676fdf`, which does **not** contain #399 —
`src/io.rs` in the working directory is the **pre-change** version (`pub fn norm_path`
with the 7-line doc, at line 45).

Consequences worth flagging to the caller:

- Any verification run against the *working tree* after that switch measures the wrong
  branch. A "green" result there says nothing about #399.
- **This is silent.** `cargo test`/`clippy` pass on both versions and the test count is
  identical (561 either way — the change adds and removes no tests), so a reviewer who
  ran the suite after the switch would see exactly the numbers they expected and draw a
  conclusion about the wrong source. **Reviewer A's numbers are worth checking for this
  same exposure.**
- I therefore re-ran the whole suite against a copy pinned with
  `git archive 315557d`, which is immune to the shared tree. All numbers in the
  Verification section below are from the pinned copy unless stated otherwise.

I did not switch branches or touch the tree, per instructions.

---

## Verification performed

Pinned copy: `git archive 315557d -- src tests test_files Cargo.toml Cargo.lock build.rs`
extracted to `<scratchpad>/probe399b`; confirmed `src/io.rs:39` reads `fn norm_path`
(private) before running anything.

| # | Check | Result |
|---|-------|--------|
| 1 | `cargo fmt --all -- --check` (pinned) | **exit 0**, no diff |
| 2 | `cargo clippy --all-targets -- -D warnings` (pinned) | **exit 0**, zero warning/error lines |
| 3 | Full `cargo test --release` (pinned) — CI's exact invocation | **exit 0, 561 passed / 0 failed.** Matches the commit message's "561 tests green" claim exactly. |
| 4 | `cargo test test_norm_path_case_folds` | **pass** — `io::tests::test_norm_path_case_folds ... ok`. Plan assumption A2 (child module reaches a private parent item) verified by *running* it, not by citing the language rule. |
| 5 | **Positive control — external caller** | Planted `tests/probe_e0603.rs` calling `trim_galore::io::norm_path` in the pinned copy. `cargo check` → **exit 0, zero E0603** (compiles lib+bin only, never sees `tests/`). `cargo check --all-targets` → **exit 101**, ``error[E0603]: function `norm_path` is private``. Probe removed → `--all-targets` back to **exit 0, 0 warnings**. The harness demonstrably both passes and fails. |
| 6 | **Is the net actually in CI?** | Yes — `ci.yml:157` runs `cargo clippy --all-targets --release -- -D warnings`, and `ci.yml:76` runs `cargo test --release`. The `--all-targets` that control 5 proves is *required* is the flag CI uses. |
| 7 | `--all-targets` scope | `cargo metadata`: 1 lib, 1 bin, 11 integration tests, 1 example. So a missed external caller anywhere in `tests/` is caught. |
| 8 | `cargo doc --no-deps` (real tree, pre-switch) | **exit 0.** `target/doc/trim_galore/io/` contains `fn.collision_key.html` and `fn.path_identity_key.html` but **no `fn.norm_path.html`** — which is itself proof this run saw the private version. |
| 9 | Caller search | `norm_path` appears **only** in `src/io.rs` (definition, the one call in `collision_key`, and the test) plus `src/cli.rs:1903`, which is a `//` line comment naming the *test* — unaffected. Searched `src/ tests/ examples/ benches/ build.rs Cargo.toml Docs/`; also confirmed `lib.rs` is 17 bare `pub mod` with **no `pub use`**, and checked the `use crate::io as naming` aliases in `specialty.rs`/`clump_only.rs` (no `naming::norm_path` anywhere). |
| 10 | Diff scope | `git diff --stat` = `src/io.rs` only. No stray files, no CHANGELOG entry — consistent with plan step 7. |
| 11 | Commit message | Accurate on every factual claim I checked, including the `pub(crate)`-rejected rationale the plan asked for and the E0603/`cargo check` asymmetry. |

### A flake I hit, chased down, and cleared — recorded so nobody re-derives it

My **first** pinned run was `cargo test` in **debug**, and it came back
**559 passed / 2 failed**:

```
test ubam_out_se_fastqc_produces_report ... FAILED
test ubam_out_pe_fastqc_produces_exactly_one_report ... FAILED
integration_ubam_out.rs:561: FastQC zip missing — --fastqc silently skipped on the uBAM-output path
```

**This is not caused by #399, and I did not take it on trust in either direction.** What
I established:

1. The failing assertion is `io.rs`-independent — it is `--fastqc` on the uBAM-output
   path. The assertion immediately *before* it (output BAM exists) **passed**, and the
   process exited zero, so trimming worked and only the FastQC artifact was absent.
2. Reproduced manually with the probe's own debug binary → **the zip and html were
   produced correctly.** So the binary is fine.
3. Same pinned source, **release**, same test target → **23/23 pass**.
4. Same pinned source, **debug, `--test-threads=1`** → **23/23 pass**.
5. Same pinned source, **debug, parallel, repeat run** → **23/23 pass**.
6. Full **`cargo test --release`** (CI's invocation) → **561/0**.

Four subsequent test runs green on identical source, plus a successful manual
reproduction attempt ⇒ a non-deterministic failure in my sandboxed probe under first-run
load, not a source defect. There is no mechanism by which a visibility keyword on a
path-lowercasing helper reaches FastQC zip creation.

**Residual observation, pre-existing and out of scope for #399** (flagged only because I
have the data): these two tests appear non-deterministic under `cargo test` in **debug**
with default parallelism. CI is unaffected — `ci.yml:76` runs `cargo test --release`. But
CLAUDE.md tells contributors to run plain `cargo test`, so a contributor could see these
two fail spuriously and go hunting. Worth its own issue if it recurs; **not** something
#399 should absorb, and I could not reproduce it a second time to characterise it
properly.

---

## Area 1 — Logic / correctness

**Clean.** Findings:

- **The demotion compiles and nothing reachable breaks.** Verified by a real
  `clippy --all-targets` on pinned source, plus the E0603 positive control proving that
  result is meaningful rather than vacuous.
- **The fold is still pinned.** `test_norm_path_case_folds` survives, calls the private
  fn from the child `mod tests`, and still asserts all five properties (lowercase
  identity, uppercase fold, mixed fold, absolute-path fold, alias equality).
- **No `dead_code` exposure.** `norm_path` retains a production caller
  (`collision_key`, `io.rs:81`), so privatization cannot orphan it. Clippy confirms.
- **The `path_identity_key` cross-reference is accurate.** `collision_key`'s new
  sentence says input identity "uses the case-preserving `path_identity_key`". Verified:
  `path_identity_key` is called at `cli.rs:537` (duplicate R1/R2 within a pair),
  `cli.rs:549-550` (duplicate pair detection), `cli.rs:635-638` (duplicate input), and
  `cli.rs:934/941` (passthrough input matching) — all input-identity questions, none of
  them output naming. The claim is true as written.
- **The trade-off paragraph is more accurate on `collision_key` than it was on
  `norm_path`.** Agreed with the plan's reasoning, and it survives scrutiny: a function
  returning `String` cannot "false-positive", but `collision_key`'s doc block is framed
  as a *predicate* ("Would these two paths be the same *output* file?"), and under that
  framing "may false-positive" means "may answer yes when the answer is no" — which is
  exactly right. On `norm_path` (a bare `to_ascii_lowercase`) the sentence was a
  category error.
- **One nuance the plan did not claim but which supports the move:** `main.rs:790` and
  `main.rs:2541` both carry the comment *"Pre-flight across pairs before any I/O; see
  `io::collision_key` for the key."* Before this change a reader following that pointer
  arrived at a two-line doc with no trade-off in it. Now the pointer lands on the
  trade-off. The move repairs two pre-existing pointers rather than creating any.

## Area 2 — Errors, rustdoc, and whether CI would catch anything

**No new rustdoc breakage, and I can show why plus show that nothing would have caught
it if there were.**

- **No bracketed intra-doc link references `norm_path` anywhere.** `src/io.rs` contains
  **zero** bracketed intra-doc links in the entire file — it uses plain backticks
  throughout. So privatization cannot have created a `private_intra_doc_links` warning.
  `cargo doc --no-deps` exits 0 and emits nothing about `io`.
- **The hazard class is real in this crate, and it is completely ungated.** `cargo doc`
  already emits 12 warnings on `dev`, **including exactly the failure mode in question**:

  ```
  warning: public documentation for `alignment` links to private item `myers_proves_no_match`
           = note: `#[warn(rustdoc::private_intra_doc_links)]` on by default
  ```

  That warning is sitting in the tree unnoticed. It confirms both plan claims at once:
  the hazard is live, and nothing catches it —
  `grep -rniE "cargo doc|rustdoc|RUSTDOCFLAGS|document-private"` over `.github/` returns
  **nothing**. `ci.yml:160`'s "Docs build (Astro)" is the documentation *website*. So
  this change is safe because of a style choice, not because of a gate.
- **The rest of the crate does use bracketed links** (`fastq.rs`, `trimmer.rs`,
  `specialty.rs`, `adapter.rs`, `bam.rs`, `alignment.rs`, `fastqc.rs`, `report.rs`,
  `main.rs`). So the plain-backtick choice matches `io.rs` locally while diverging from
  the crate. See Low-2.
- **Doc-tests:** no doc example anywhere references the symbol, so the `--all-targets`
  blind spot the plan flagged (A5) stays theoretical. Full `cargo test` green anyway.
- **Semver:** removing a `pub fn` from a crate published at 2.3.0 is strictly a breaking
  API removal. The plan owns this (A3) and the risk is negligible — the library surface
  is uncurated (17 bare `pub mod`, no `pub use`, no `publish = false`), and
  `release.yml:442` publishes with `--no-verify`. No `cargo-semver-checks` in CI, so
  nothing will object. No action for this PR; the plan's filed follow-up (a `//!` note
  declaring the surface an implementation detail) is the right durable fix.

## Area 3 — Structure / style

**Comment budget: the change is net-negative in comment lines** (`norm_path` −6,
`collision_key` +5, test −1). Worth stating plainly because CLAUDE.md's convention
exists to stop comment growth, and this diff shrinks it while relocating.

On the convention itself ("one line default, two max, state the fact not the
evidence"): `collision_key`'s doc block is now 6 lines, which a literal reading would
flag. I do **not** think that is the right reading here. That rule governs inline `//`
comments in code bodies; `io.rs`'s established convention for `///` API blocks is
substantially longer — `is_gzipped`, the function **immediately above** `norm_path`,
carries a 21-line doc block (`io.rs:12-32`), and several test items carry 8–12-line
blocks. The change also moves an existing paragraph rather than authoring new prose.
Consistent with the file.

The two rewritten comments:

- `norm_path`'s reduced one-liner is **compliant and clear**, and correctly omits the
  "(private)" framing the plan dropped — the `fn` keyword on the next line says it.
- The test comment is now 2 lines (down from 3), within budget. But see Low-1 and Low-3.

**Did anything become less clear?** One thing, and it is the only place I would say yes:
the old `norm_path` doc contained *"`collision_key`, **which absolutises first**"*. That
clause is gone. It was load-bearing for a specific reader — the one who asks "why is
there a separate test for the fold instead of just testing `collision_key`?" The answer
(the plan's own surviving rationale, B §1.4: `collision_key` absolutises, so re-expressing
the fold test through it would make the relative-path assertions cwd-dependent) is no
longer stated anywhere near either the fold or its test. It is recoverable —
`lexical_normalise`'s doc four lines below says "absolutised" — so this is Low, not
High. Noted under Low-3.

## Area 4 — Dangling references

**None.** Checked exhaustively:

- `grep "Pragmatic trade-off"`, `grep "silent data loss"`, `grep "loud early error"`,
  `grep -i "opt-in case"` over `src/` each return **exactly one hit**, all at the new
  location (`io.rs:77-79`). The paragraph was moved, not copied — no second stale copy
  survives. (I chased this specifically: an earlier review round
  described this paragraph as "duplicated verbatim", so a leftover copy was a live
  possibility. There isn't one.)
- No comment anywhere points a reader at `norm_path` for the trade-off.
- `cli.rs:1903` names `io::tests::test_norm_path_case_folds` — the test still exists
  under that exact name, so the pointer is intact.
- `io.rs:1210`'s test doc says outputs are "compared on `collision_key` (the pre-flight's
  own metric)" — still true.
- **Issue-reference coverage survives the move.** `#389` was dropped from `norm_path`'s
  doc, but it lives on at `io.rs:74` on `path_identity_key` — which is the *correct*
  home, since #389 was an input-identity confusion. `#216`/`#383` remain on
  `collision_key`. Nothing was orphaned.

## Area 5 — Efficiency

**Nil, as expected.** Visibility is an annotation with no codegen consequence; the sole
call site was already intra-crate. If anything privatization marginally *helps* codegen
(no external linkage to preserve on a one-line function), which is unmeasurable and not
a reason for the change.

---

## Recommendations by priority

### Critical — none
### High — none

### Medium

**M-1. The relocated trade-off names the rarest instance of its own condition, and now
does so on a public page.**

`io.rs:77` reads *"on opt-in case-sensitive APFS volumes this may false-positive"*. The
condition for a false positive is **any case-sensitive filesystem** — which includes
**every Linux run**, i.e. most TrimGalore invocations and all of CI. Naming only "opt-in
case-sensitive APFS volumes" makes a routine, deliberately-accepted behaviour sound like
an exotic edge case.

This is not speculation — the repo's own test says so. `io.rs:1170-1174`:

> REGRESSION: re-keying `Cli::validate`'s input-identity checks on the case-folded key
> made the #216 CI guard fail — it feeds four genuinely distinct files (`Sample_R1` /
> `SAMPLE_R1`) **on a case-sensitive filesystem** and asserts the *output* pre-flight
> refuses them.

So CI exercises this false positive on Linux on every run, by design, as the safe
default. The wording is **pre-existing and unchanged by this diff** — but this change is
what promotes it from a private helper's doc to `fn.collision_key.html`, a rendered
public API page, which makes now the cheapest possible moment to correct it.

*Trivial fix* — one word of scope:

```rust
/// Pragmatic trade-off: on case-sensitive filesystems (Linux, or opt-in
/// case-sensitive APFS) this may false-positive, but the penalty is a loud
/// early error rather than silent data loss.
```

Same line count. I did not apply it — it is prose the author may want to word
differently, and the plan explicitly left doc copy reviewable.

### Low

**L-1. The test comment's first sentence is a fragment and reads awkwardly.** *(trivial
fix)*

```rust
// The fold `collision_key` is built on (callers reach it through that key,
// never directly). Plain lowercase is identity; upper/mixed folds down.
```

"The fold `collision_key` is built on" is a noun phrase with a trailing preposition and
no main verb. The old version was grammatical ("norm_path is the case-folded
normalisation that …"). Terser is good; this crossed into needing a second read.
Suggestion, same 2 lines:

```rust
// The fold that `collision_key` is built on; production callers only ever
// reach it through that key. Lowercase is identity; upper/mixed fold down.
```

**L-2. Consider making the `path_identity_key` reference a real link.** *(optional —
author's call, genuine trade-off both ways)*

`collision_key`'s new sentence exists to send a reader to `path_identity_key`. Both
items are `pub` and rustdoc renders both pages (I verified `fn.path_identity_key.html`
exists), so ``[`path_identity_key`]`` would resolve to a working hyperlink — exactly
where the plan wants public readers to be able to go. The plan's stated reason for plain
backticks (nothing in CI would catch a broken or private link) is a sound default, but
does not apply to this particular link, whose target I confirmed is public and rendered.

Counter-argument, which is why this is Low and not a recommendation: `io.rs` has zero
bracketed links today, so this would be the first, and 8 other modules using the style
does not oblige this one. Leaving it is defensible; I'd lean to leaving it for local
consistency.

**L-3. The dropped "which absolutises first" clause is the one real clarity loss.**
*(optional)*

See Area 3. If it matters to the author, the cheapest home is the test comment, since
the fact's main use is explaining why the test targets the fold directly rather than
going through `collision_key`. Not worth a line on its own if L-1 is taken instead.

**L-4. "callers reach it through that key, never directly" is contradicted one line
later.** *(trivial fix; subsumed by L-1)*

The very next statement in the same test calls `norm_path` directly. "Callers" plainly
means production callers, and L-1's rewording ("production callers") resolves it.

---

## Comparison notes for the caller

Things I would specifically check against Reviewer A's report:

1. **Whether A's verification ran against the switched tree.** See the tree-hazard
   section — this failure is silent and produces the expected numbers.
2. **Whether A also concluded the trade-off sentence is accurately homed on
   `collision_key`.** I think it is, on the predicate-framing argument; a reviewer could
   reasonably argue it belongs on `preflight_output_collisions`, where the bail actually
   happens. I considered and rejected that: `collision_key` also serves `cli.rs:935`'s
   passthrough check, so the key's doc covers both consumers while the pre-flight's would
   cover one.
3. **M-1 (scope of the false-positive condition)** — this is the only substantive
   finding I have, and it is easy to miss because the wording is unchanged by the diff.
4. **Whether A saw the two `integration_ubam_out` `--fastqc` failures.** If A ran
   `cargo test` in debug they may have, and they are alarming at first glance
   ("--fastqc silently skipped"). They are a flake, not #399 — see the section above for
   the five runs that clear it. If A reported them as a defect, that is a false positive;
   if A ran only `--release`, they will not have seen them at all.
