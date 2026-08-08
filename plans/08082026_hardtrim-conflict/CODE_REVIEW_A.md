# Code Review A — Reject `--hardtrim5` + `--hardtrim3` together (#386)

**Target:** `fix/386-hardtrim-conflict` off `dev` @ `c4599f5`, uncommitted (`src/cli.rs`, `CHANGELOG.md`).
**Reviewer:** single reviewer. The usual dual-reviewer default was scaled down for diff size on the
orchestrator's call (one clap attribute, two tests, one CHANGELOG entry).
**Nothing was changed.** Recommendations only — the tree is uncommitted.

## Summary

Correct, minimal, and in the right layer. The exclusion is real (verified against the binary in both
flag orders), no other code path depends on both flags being set, nothing shipped in the repo passes
both, and the CHANGELOG's behavioural claim checks out against `main.rs`. No Critical or High findings.
Three Low and two Medium items, all optional polish; the Mediums are a missing remediation hint in the
error text and the docs page's `## Compatibility` section, which now omits the only hard
incompatibility the mode has.

## Verification performed

**Parse-time behaviour** (release binary, `c4599f5`, rebuilt 07:51 after the 07:50 `cli.rs` edit; run
from a scratch cwd, never the repo root):

| Invocation | Result |
| --- | --- |
| `--hardtrim5 20 --hardtrim3 15 F` | exit 2, stderr: `the argument '--hardtrim5 <HARDTRIM5>' cannot be used with '--hardtrim3 <HARDTRIM3>'` |
| `--hardtrim3 15 --hardtrim5 20 F` | exit 2, same message with the names swapped — **symmetry confirmed empirically**, not assumed |
| `--hardtrim3=15 --hardtrim5=20 F` | exit 2, rejected in `=` form too |
| `--hardtrim5 20 --hardtrim3 0 F` | exit 2 — the conflict wins over the `1..999` range check in `validate()`; still a usage error, message names both flags |
| `--hardtrim5 20 F` | exit 0, writes `…​.20bp_5prime.fq.gz` |
| `--hardtrim3 15 F` | exit 0, writes `…​.15bp_3prime.fq.gz` |
| `--help` | exit 0, both entries render, no conflict noise |

**Layer question (brief item 1).** Parse-time is the whole surface for the binary: `main.rs:175` is the
only `Cli::try_parse_from` call, `Cli` is never constructed programmatically anywhere in `src/` or
`tests/`, there are no `env = …` attributes, neither flag has a `default_value`, and there is no
config-file or `@argfile` input path. `rewrite_perl_short_flags` only maps `-r1`/`-r2`/`-a2`, so it
cannot manufacture or hide either flag. `main.rs`'s parse-error handler special-cases only
`DisplayHelp`/`DisplayHelpOnMissingArgumentOrSubcommand`/`DisplayVersion` and sends everything else to
`err.exit()` — so `ArgConflict` correctly lands on stderr with exit 2, as measured.

**Sweep of every `hardtrim3` reader (brief item 2).** `main.rs:323-324` (paired-mode exemption,
`hardtrim5.is_none() && hardtrim3.is_none()`), `main.rs:394`/`:425` (dispatch), `cli.rs:718-722`
(`--clumpify`), `:851-855` (`--clump_only`), `:921-925` (`--passthrough`), `:952-960` (range check),
`:986-987` (the `-a2`-unusable reason, `hardtrim5.is_some() || hardtrim3.is_some()`). None requires
both to be set, and none changes meaning when only one can be: the `&&`/`||` forms stay correct with a
strictly narrower input domain. No flag implies or defaults `hardtrim3`.

**Nothing shipped breaks.** No test, CI step, README, or docs page passes both flags.
`tests/integration_adapter2.rs:375-382` has cases *named* `hardtrim5`/`hardtrim3` but each passes one
flag. `ci.yml:658-768` (collision pre-flight, Perl md5 parity) uses `--hardtrim5` alone.
`tests/integration_output_collision.rs:329-431` uses one flag per case.

**CHANGELOG accuracy (brief item 5).** Verified independently of the filing: `main.rs:394-422`'s
`if let Some(n) = cli.hardtrim5` block ends in an unconditional `return Ok(())`, so the 3' request was
dropped *regardless of CLI order* — the entry's "only the 5' trim ran" is exact, not order-dependent.
The Perl claim is taken as given per the brief.

**Gates.** `cargo fmt --all -- --check` clean (re-run, not taken on trust). `cargo test --lib hardtrim`
→ 6/6 pass, including both new tests and the four pre-existing `hardtrim` guards. Tests run in the
debug profile, so clap 4.6.4's arg-graph `debug_assert` is live — a typo'd conflict ID would panic
rather than silently no-op; the ID `hardtrim3` is the derive field name and is correct.

## Findings

### Medium

**M1 — the error carries no remediation, unlike every comparable refusal in this codebase.**
clap's generic "cannot be used with" tells the user what is forbidden but not what to do. The
CHANGELOG entry does say it ("Run the two trims as separate invocations") — the diagnostic doesn't.
That is out of step with the surrounding work: `#379`'s non-restartable message ships a remediation
command, `#383`/`#385`'s collision messages name `--output_dir`, and the *other* hardtrim exclusions
bail from `validate()` with bespoke text (`--clump_only and --hardtrim5 are mutually exclusive`,
`cli.rs:851-855`). Moving the check into `Cli::validate()` §3.4a alongside those would let the message
carry the one-sentence fix.
Trade-off, and I'd call it genuinely close: clap gives correct symmetric wording for free (measured
above) and matches the 7 existing `conflicts_with` uses; `validate()` costs a hand-written both-order
check plus tests for both orders, and loses the free symmetry. If the remediation sentence is not
worth that, keeping clap is defensible — but then the CHANGELOG is the only place the fix is written
down, and users hit the error before they read the CHANGELOG. **Recommend, user's call.**

**M2 — `docs/src/content/docs/modes/hardtrim.md` documents both modes and, at line 36-38, has a
`## Compatibility` section that now omits the only hard incompatibility.** Its current content
("Both modes accept multiple input files … operates per-file independently") reads as permissive; a
user could reasonably infer the two flags compose. One sentence at line 38 closes it. Plan §2 scoped
docs out, so this is a follow-up rather than a defect in the diff — but the site auto-deploys from
`dev`, so it ships as soon as this merges.
Separately: `docs/src/content/docs/reference/changelog.md` is a hand-synced copy and already lags
(`#383`/`#384` are absent, file last touched Aug 3). **Pre-existing, not this change's regression** —
noted only so it isn't mistaken for one.

### Low

**L1 — `test_each_hardtrim_alone_still_accepted` uses `Cli::parse_from`, which calls
`std::process::exit(2)` on a parse error rather than returning.** If a future over-broad conflict
breaks the alone case, this test kills the whole test binary instead of failing one test — a
confusing failure mode for the exact regression the test exists to catch.
`try_parse_from(..).expect("…")` fails cleanly and matches the sibling test's style. Counter-argument
that keeps this at Low: `parse_from` is the established idiom in this module (10+ call sites,
`cli.rs:1148` onward), so the new test is consistent with the file. The `validate()` call itself is
right — it does run, and the fixture path resolves because `cargo test` runs from the crate root.

**L2 — `--help` doesn't state the exclusion.** The newer flags spell theirs out (`--clump_only` and
`--passthrough` both enumerate incompatibilities in their doc comments); the five adapter presets set
the opposite precedent (mutually conflicting, none documented). A clause on `--hardtrim5`'s doc
comment would settle it in the newer direction. Optional.

**L3 — the exclusion is not enforced for a programmatically-built `Cli`.** `pub struct Cli` with
public fields in a library crate means a consumer (or future in-repo caller) can set both and get the
old silent drop back. No live path today — `main.rs` is the only dispatcher — so this is only worth
acting on if M1 moves the check into `validate()`, which closes it for free.

## Recommendation

Ship as-is if M1 is judged not worth the churn; the change is correct either way. The two items I'd
actually action are **M2** (one sentence in `modes/hardtrim.md`) and, if the check moves to
`validate()`, **M1 + L3** together. Nothing here blocks the commit.
