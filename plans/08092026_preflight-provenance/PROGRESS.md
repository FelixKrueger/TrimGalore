# Progress: the collision pre-flight names the inputs at fault (#397)

**Last updated:** 2026-08-09

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | PLAN.md — recommends `planned: &[(PathBuf, OutputSource)]` (tuple slice + `PlannedOutput` alias, no builder) and three duplicate-output message shapes; 11 physical call sites / 13 logical dispatch paths + 4 helpers enumerated; full test/CI pin inventory (every pinned substring preserved); 4 Critical questions incl. **C4 — live data loss on the SE uBAM arm found by Probe 3**; awaiting Felix's answers |
| Plan Review | 📋 Planned | Dual independent reviewers once Felix has done the manual pass |
| Impl Plan | 📋 Planned | Expect to fold into PLAN.md per pipeline convention rather than a separate IMPL.md |
| Implementation | 📋 Planned | Blocked on C1-C4 |
| Code Review | 📋 Planned | |
| Coverage | 📋 Planned | |

## Open decisions for Felix

- **C1** — is message shape 2c (a positional reused across pairs) in #397's scope? Recommend yes — the case is already reachable and pinned (`tests/integration_output_collision.rs:1215`), it just gets unfollowable advice today.
- **C2** — paired primaries carry `Pair(r1, r2)` uniformly, or the more pointed `Input(r1)` for `_val_1`? Recommend uniform `Pair`; this is the single design point a #398 answer would change.
- **C3** — the `--clump_only --paired` uBAM per-pair report keeps `Input(chunk[0])` (truthful) despite the asymmetry with its sibling `Pair` candidate. Needs blessing.
- **C4** — Probe 3 found live data loss at `90dda11`: `--output-format ubam sample.fastq sample.fastq_trimming_report.txt` overwrites the second input before failing to read it. Separate issue first, or fold into #397's PR?

## Related

- **#398** (undecided) — paired-arm report/primary directory asymmetry. PLAN.md §Integration traces all four contact points; only **C2** would change if #398 is answered first.

## History

- 2026-08-09: C1-C4 answered (2c in, uniform Pair, Input(chunk[0]) blessed, C4 split to #409/PR #410); plan unblocked, must rebase onto merged dev

- 2026-08-09: Plan → ✅ Complete (PLAN.md written from a live investigation of `dev` @ `90dda11`, five probes against the current debug binary; presented for manual review)
