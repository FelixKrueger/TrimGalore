# Code Review — Build Provenance Core (Reviewer A)

**Scope**: `build.rs`, `src/cli.rs` (`#[clap(...)]` on `Cli` only), `src/main.rs` (startup banner only).
**Plan**: `plans/build-provenance/PLAN.md` §build.rs, §clap version, §banner.
**Verdict**: Core logic is correct and meets the plan. Two minor observations, no blockers.

---

## Summary

All six mandated checks pass. `SOURCE_DATE_EPOCH` hard-fail is implemented with a clear `panic!` message. Git fallback safely degrades to `"unknown"` across every failure mode (missing git, missing `.git`, non-zero exit, non-UTF-8, empty stdout). All three rerun directives are present. `-V` remains single-line; `--version` carries the provenance body. The banner is emitted exactly once, to stderr, before any trim work.

---

## Checklist results

### 1. `SOURCE_DATE_EPOCH` hard-fail — PASS
- `build.rs:18-23`: `env::var("SOURCE_DATE_EPOCH")` → `Ok(s).parse::<u64>().unwrap_or_else(|_| panic!(...))`. Message is descriptive and includes the offending value via `{s:?}` (quoted debug form). Negative inputs (e.g. `-1`) fail `u64::from_str` and trigger the panic, as required.
- Unset env var cleanly falls through to `SystemTime::now()`. No silent fallback.

### 2. Git fallback → `"unknown"` — PASS
- `build.rs:5-15`: the chain
  `.ok()` (Command launch failure, e.g. git not on PATH)
  → `.filter(|o| o.status.success())` (non-zero exit, e.g. not a repo / shallow with no HEAD)
  → `.and_then(|o| String::from_utf8(o.stdout).ok())` (non-UTF-8, unlikely but defensive)
  → `.map(|s| s.trim().to_string())`
  → `.filter(|s| !s.is_empty())` (empty stdout)
  → `.unwrap_or_else(|| "unknown".to_string())`
  covers every documented failure mode in PLAN §6.1 (missing `.git`, git not on PATH, `cargo install` tarball path). Good layered defence.

### 3. Rerun directives — PASS
- `build.rs:70-72`: `.git/HEAD`, `.git/index`, and `SOURCE_DATE_EPOCH` are all emitted verbatim as the plan specifies.
- Note (non-blocking): `.git/HEAD`'s *content* only changes on branch switch or detached-head moves, not on commits-to-current-branch. The `.git/index` trigger covers the staging side, and the final commit updates `.git/refs/heads/<branch>` — which this file does NOT watch. Practically this still works because `.git/index` is rewritten on `git commit` (index is refreshed), so rebuild triggers. The plan's rationale holds, but if a user commits via a tool that doesn't refresh index (rare), hash staleness is theoretically possible. Matches ruSTAR's approach; acceptable.

### 4. `-V` short flag stays terse — PASS
- `src/cli.rs:13`: `version = concat!(env!("CARGO_PKG_VERSION"), " (Oxidized Edition)")` — single line, no provenance. Clap maps `-V` to `version`, so the short flag is unchanged from the pre-feature behaviour. Plan goal §2 satisfied.

### 5. `--version` long flag carries provenance — PASS
- `src/cli.rs:14-17`: `long_version = concat!(env!("CARGO_PKG_VERSION"), " (Oxidized Edition)\n", env!("VERSION_BODY"))`. Clap maps `--version` to `long_version` when present. Output will be the required two-line form. `env!("VERSION_BODY")` resolves at compile time via the `cargo:rustc-env=` directive in `build.rs:68`.

### 6. Banner emits exactly once, before trim work, to stderr — PASS
- `src/main.rs:29-34`: three `eprintln!` calls immediately after `cli.validate()?` and before `FastqReader::sanity_check`. All specialty-mode dispatches (`hardtrim5`/`hardtrim3`/`clock`/`implicon`) and the main trim path begin *after* line 34, so the banner precedes every trim path including early-exit specialty modes.
- `eprintln!` writes to stderr — correct stream.
- No loop or conditional around the banner block — emits exactly once per process invocation.

---

## Findings

### Low — `.git/HEAD` vs `.git/refs/heads/*` rebuild-trigger gap (informational)
`cargo:rerun-if-changed=.git/HEAD` tracks the HEAD ref-pointer, not the commit that ref resolves to. On a normal `git commit` on the current branch, `.git/HEAD`'s content is unchanged; the update lands in `.git/refs/heads/<branch>`. The `.git/index` trigger is what actually catches fresh commits (index is rewritten on commit). This is the same pattern ruSTAR uses and is documented to work in practice — just flagging that the plan's rationale ("a fresh commit triggers a rebuild") is achieved via the index trigger, not the HEAD trigger. No action needed.

### Low — Panic message wording nit
`build.rs:21-22`: the panic message says "non-negative decimal seconds-since-epoch integer" which is accurate, but the PLAN §6.6 draft wording ("must be a u64 seconds-since-epoch") differs slightly. Both convey the same constraint; the implementation's wording is actually clearer. No change needed.

### Informational — ISO 8601 formatter is dependency-free and correct
`build.rs:31-56`: `format_iso8601_utc` + `civil_from_days` (Howard Hinnant's algorithm) correctly converts u64 epoch seconds to `YYYY-MM-DDTHH:MM:SSZ` without pulling in `chrono`. Spot-checked algorithm: for `epoch=1700000000` (2023-11-14T22:13:20Z), the math yields `days=19675`, `secs_of_day=80000 → 22:13:20`, and `civil_from_days(19675)` produces `(2023, 11, 14)`. Matches expected output. No dependency bloat, aligns with the plan's "no `chrono`" constraint.

### Informational — Banner emits before sanity check — by design
`src/main.rs:29-36`: the banner precedes `FastqReader::sanity_check(&cli.input[0])`. If sanity fails, the banner is still printed. This is the right behaviour for user-facing identification (version/provenance is useful context in bug reports even when input is malformed).

---

## Fixes applied

None — no defects at a severity warranting direct fix.

## Recommendations

None blocking. Optional future work (not for this PR):
- If the `.git/HEAD`-vs-refs gap ever surfaces in practice, add `cargo:rerun-if-changed=.git/refs/heads` (directory watch) for completeness.
