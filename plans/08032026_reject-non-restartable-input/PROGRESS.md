# Progress: Reject non-restartable input with an accurate diagnostic

**Last updated:** 2026-08-05

## Status

| Step | Status | Notes |
|------|--------|-------|
| Plan | ✅ Complete | `PLAN.md`, revised twice — v1 → v2 → v3 |
| Plan Review | ✅ Complete | Manual review, then dual agents: `PLAN_REVIEW_A.md`, `PLAN_REVIEW_B.md` |
| Impl Plan | ✅ n/a | The design plan carries its own implementation outline (§5) |
| Implementation | ✅ Complete | §12 implementation notes; 470 tests pass, fmt + clippy clean |
| Code Review | ✅ Complete | Dual agents: `CODE_REVIEW_A.md`, `CODE_REVIEW_B.md`. All findings taken — §13 |
| Coverage | ✅ Complete | `COVERAGE.md` — **COMPLETE**, 62 items, 0 missing, 0 partial |

**Open before merge:** V3 — byte-identity against Perl Trim Galore 0.6.11 needs the CI `validation` / `validation-ubam` jobs. Cannot be run locally.

## Issue

[#379](https://github.com/FelixKrueger/TrimGalore/issues/379) — `<(zcat reads.fq.gz)`, `/dev/stdin` and FIFOs could not be read. Scope: option **A** of three — reject accurately, do **not** add pipe support (§10 costs B and C).

## What shipped

Three layers, each covering what the others cannot:

1. **`Cli::validate`** — a path `stat` (replacing `exists()`, so zero extra syscalls) names a pipe, FIFO or socket before anything opens it. This is the only layer that can diagnose a **writer-less FIFO**, because `fs::metadata` returns where `File::open` blocks forever.
2. **`detect_input_format`** — an `fstat` on the handle it already holds, before the first `read`.
3. **The second-open check** — not restartable if `stream_position() != 0` **or** a length-aware leading-byte comparison fails.

Plus `--passthrough` routed through the guard, a uBAM `--passthrough` rejected up front, and the `--demux` barcode file guarded against the blocking-open case.

## The findings that shaped it

Each was a claim verified against a narrower set than it named — the recurring failure of this plan.

| Stage | Finding |
|---|---|
| Before v1 | A **seekability probe** — what #379 proposes — accepts `/dev/stdin` redirected from a regular file, which is seekable and still does not restart (`/dev/fd/N` is `dup` on Darwin). Killed by measurement |
| Manual review of v1 | The byte comparison **alone hangs**: its own second open blocks on a FIFO whose writer has exited, the commonest FIFO shape. A design defect, recorded in v1 only as a test risk. Fixed by the fstat pre-filter |
| Dual plan review of v2 | Both reviewers independently found the guard **missed `--passthrough` entirely**, and that v2 stamped "Verified" on a claim checked against `cli.input` only |
| Dual code review of v3 | Both reviewers independently found the layer-0 change **regressed `--passthrough <uBAM>`** from a pre-output error to a deep failure with three empty files on disk. §3.2's justification was false and §2.4 already said why |

## Why the tests are shaped the way they are

`File::open` on a writer-less FIFO blocks, so a careless test hangs CI rather than failing it — which reads as a slow job, not a red X. Mitigations, all in §5.2 and §5.3:

- every FIFO test bounds the call under test in-process (`recv_timeout`), because `timeout(1)` exists on neither Darwin nor `macos-latest`;
- the integration tests use `spawn` + `kill`, never `Command::output()`, which cannot be interrupted;
- `timeout-minutes: 45` on the job and a nextest `terminate-after` — the structural backstop that survives an implementer skipping every rule.

Proven, not assumed: forcing the type check to `false` makes the FIFO test **fail in 10 s** with `blocked; it must fail, not hang`.

## History

- 2026-08-03: Plan → ✅. Scope A confirmed after three options were costed
- 2026-08-04: v2 after manual review (fstat pre-filter); dual plan review; v3 (`--passthrough`, offset predicate, CI backstop); implementation; dual code review + coverage audit; fix round — all findings taken
- 2026-08-05: Committed on `fix/reject-non-restartable-input`
