# Build Provenance — Scoped Coverage Audit

**Plan:** `/Users/fkrueger/Github/TrimGalore/plans/build-provenance/PLAN.md`
**Branch:** `optimus_prime` (uncommitted working-tree diff)
**Date:** 2026-04-19

| Item | Verdict | Evidence |
|---|---|---|
| 1. `build.rs` exists at repo root | DONE | `build.rs:1` present |
| 2. emits `cargo:rustc-env=GIT_SHORT_HASH` | DONE | `build.rs:66` |
| 3. emits `cargo:rustc-env=BUILD_TIMESTAMP` | DONE | `build.rs:67` |
| 4. emits `cargo:rustc-env=VERSION_BODY` | DONE | `build.rs:68` |
| 5. panics on malformed `SOURCE_DATE_EPOCH` | DONE | `build.rs:19-23` `panic!` on parse err |
| 6. git-failure fallback to `"unknown"` | DONE | `build.rs:14` `.unwrap_or_else("unknown")` |
| 7. rerun triggers: HEAD, index, SOURCE_DATE_EPOCH | DONE | `build.rs:70-72` |
| 8. `src/cli.rs` `version = ...` terse | DONE | `src/cli.rs:13` |
| 9. `src/cli.rs` `long_version = ...` with VERSION_BODY | DONE | `src/cli.rs:14-17` |
| 10. banner prints to stderr before trim work | DONE | `src/main.rs:29-34` before `sanity_check` |
| 11. `Cargo.toml` — no version bump | DONE | `Cargo.toml:3` still `2.1.0-beta.1`; no `build =` key |
| 12. `lint` job (fmt-check + clippy -D warnings) | DONE | `ci.yml:84-104` |
| 13. `validation.needs` includes `lint` | DONE | `ci.yml:122` `needs: [rust-tests, lint]` |
| 14. `reproducibility` job with dedicated cache | DONE | `ci.yml:47-82` key `...-cargo-repro-...` |
| 15. `audit` job using `rustsec/audit-check@v2` | DONE | `ci.yml:106-117` |
| 16. `--version` content-regex verification | DONE | `ci.yml:37-45` regex `^[0-9a-f]{7,40} — ...` |
| 17. dependabot.yml — cargo + gh-actions, weekly Mon, limit 5, FelixKrueger | DONE | `.github/dependabot.yml:14-32` |
| 18. `release.yml` header comment: intentional non-gating | DONE | `release.yml:15-22` "Lint-gating note" block |

## Summary

- Total: 18
- DONE: 18
- PARTIAL/MISSING/DEVIATED: 0

VERDICT: COMPLETE
