# EPIC — Pluggable I/O formats for TrimGalore + Bismark

**Tracking issues:** [TrimGalore#315](https://github.com/FelixKrueger/TrimGalore/issues/315), [Bismark#1025](https://github.com/FelixKrueger/Bismark/issues/1025). Spun out of [TrimGalore#316](https://github.com/FelixKrueger/TrimGalore/issues/316) (uBAM input, now landed in TrimGalore#317) once the four-format conversation (uBAM, BINSEQ, mim, paraseq) made clear that the shape under all of them is a single pluggable I/O abstraction, not four bespoke features.

## 1 · Goal

Commit both TrimGalore and Bismark to a **symmetric pluggable I/O architecture** — one input front-end with N reader backends, one output back-end with N writer backends — so that uBAM, BINSEQ, and mim each become an opt-in backend behind the same trait surface rather than a bespoke feature flag.

Concrete user-facing outcomes:

- **TrimGalore** reads uBAM (✅ already landed via #317), BINSEQ, and mim; emits FASTQ today plus opt-in BINSEQ, mim-sidecar, and uBAM outputs.
- **Bismark** reads uBAM, BINSEQ, and mim; existing alignment-output behaviour is unchanged (byte-identity to v0.25.1 preserved).
- **TrimGalore → Bismark** pipeline can run with any format pair without a FASTQ round-trip in between.

## 2 · Scope

### In scope

- Trait-based pluggable input + output abstractions in both tools.
- uBAM, BINSEQ, and mim as concrete backends.
- Cross-repo coordination on the trait shape (shared crate vs. mirrored definitions — to be decided in Phase 3).
- Cargo-feature gating so slim builds remain lean.
- CI matrix covering each backend per tool plus the end-to-end TrimGalore-emits / Bismark-reads handshake for each format.

### Out of scope (this epic)

- paraseq integration (parked per [TrimGalore SPIKE_paraseq.md](../paraseq-spike/spikes/SPIKE_paraseq.md) findings; revisit if the architectural rewrite of `parallel.rs` ever lands).
- CRAM input (different beast — would need reference fetching, breaks the single-static-binary invariant).
- Bismark's alignment-output formats (separate question; keep byte-identity to v0.25.1).
- BAM output as a primary trim_galore mode (only as an opt-in carrier-of-provenance writer; FASTQ stays the default).

## 3 · Phase breakdown

Three phases. Phase 1 lands first (foundation + the natural extension of work already done in #317). Phase 2 can begin in parallel once the trait shape is stable. Phase 3 is the integration moment.

### Phase 1 — TrimGalore uBAM output (re-scoped in PLAN v2)

[`phase1-trimgalore-formats/PLAN.md`](phase1-trimgalore-formats/PLAN.md)

**Re-scoped 2026-06-25** after dual plan-review surfaced critical issues with the original BINSEQ/mim + `RecordSink`-trait design (see `phase1-trimgalore-formats/PLAN_REVIEW_{A,B}.md`).

- **uBAM input:** ✅ done in #317.
- **uBAM output:** opt-in via `--output-format ubam`, single-threaded, **one interleaved BAM** for paired output (matching samtools/Picard/fgbio/CellRanger convention). Aux-tag propagation via the existing `--preserve-tags` flag.
- **NO trait abstraction.** FASTQ output stays on the existing `parallel.rs` per-worker-compress + concatenate model unchanged. uBAM output dispatches to a separate serial code path with a concrete `BamWriter`.

**Deferred from Phase 1** (rolled into a future "Phase 4 — additional formats" once crate ecosystems mature):

- BINSEQ in/out — `binseq` crate is currently unstable (15 versions in 9 months, two breaking format revs) and pulls `zstd-sys` with a C `build.rs` that breaks the reproducibility invariant.
- mim sidecar output + mim opportunistic input — `mim-index 0.1.1` requires Rust 1.91; project MSRV floor is 1.88. Hard showstopper until either side moves.

### Phase 2 — Bismark uBAM input

**Status: ✅ DONE in Bismark (merged 2026-06-25).** No `phase2-bismark-formats/PLAN.md` was written in this repo — work was specified + executed entirely in the Bismark repo under tracking issue [Bismark#1025](https://github.com/FelixKrueger/Bismark/issues/1025). PRs:

- [Bismark#1026](https://github.com/FelixKrueger/Bismark/pull/1026) — `feat(aligner): single-end unaligned-BAM (uBAM) input` (merged 2026-06-25T20:36:48Z)
- [Bismark#1027](https://github.com/FelixKrueger/Bismark/pull/1027) — `feat(aligner): paired-end unaligned-BAM (uBAM) input` (merged 2026-06-25T21:59:09Z)
- [Bismark#1028](https://github.com/FelixKrueger/Bismark/pull/1028) — docs: uBAM input + Milestones line (merged 2026-06-25T22:20:02Z)

Bismark inherited the same `RecordSource` trait shape that has been stable in TrimGalore since [#317](https://github.com/FelixKrueger/TrimGalore/pull/317). Byte-identity invariant on post-alignment output preserved (input front-end is fully additive — no Perl equivalent to match). The shared-vs-mirrored-trait decision was deferred to Phase 3 (so far it's mirrored).

**Bismark phase numbering ≠ TrimGalore EPIC phase numbering.** Bismark calls this work "#1025 Phase 1" in its own PR titles; in this EPIC it's Phase 2. Future-stale-read guard.

### Phase 3 — Cross-tool integration + shared-trait decision

`phase3-cross-tool-integration/PLAN.md` _(to be written)_

- **Shared-trait decision** ratified (shared crate vs. mirrored definitions). Whichever lands, Phase 1 + Phase 2 are updated to match.
- **Pipeline integration tests:** end-to-end matrix where TrimGalore emits each output format → Bismark reads it; assert record-count parity and downstream alignment correctness on a small fixture.
- **Documentation pass across both repos:** usage examples, format-compatibility matrix, performance notes from the cross-tool benchmark.
- **Release coordination:** when the two tools cut their next minor versions (TrimGalore v2.3, Bismark v0.26 or equivalent), the feature lands paired so users get a working ecosystem in one upgrade, not two staggered ones.

## 4 · Sub-plan table

| # | Phase | Plan file | Status | Depends on |
|---|-------|-----------|--------|------------|
| 1 | TrimGalore uBAM output (re-scoped) | [`phase1-trimgalore-formats/PLAN.md`](phase1-trimgalore-formats/PLAN.md) — v2.1 | ✅ **Done** — TrimGalore [#320](https://github.com/FelixKrueger/TrimGalore/pull/320) + [#322](https://github.com/FelixKrueger/TrimGalore/pull/322) fix, on `dev` 2026-06-26 | — |
| 2 | Bismark uBAM input | _(no in-repo plan; spec + work both lived in Bismark)_ | ✅ **Done** — Bismark [#1026](https://github.com/FelixKrueger/Bismark/pull/1026) + [#1027](https://github.com/FelixKrueger/Bismark/pull/1027) + [#1028](https://github.com/FelixKrueger/Bismark/pull/1028), merged 2026-06-25 | #1 (dispatch pattern) — but actually independent in practice; landed same week |
| 3 | Cross-tool integration | _(to be written)_ | 🟡 **Unblocked** — both halves now exist | #1 and #2 |
| 4 | BINSEQ + mim formats (deferred from Phase 1) | _(to be written; TrimGalore-side gated — see §8)_ | ⏸ **TrimGalore-side deferred** per §8 gates. **🟡 Bismark side: BINSEQ INPUT already merged** via [Bismark#1029](https://github.com/FelixKrueger/Bismark/pull/1029) 2026-06-26 (Bismark accepted zstd-sys trade-off; TrimGalore's single-static-binary invariant has a tighter bar) | independent — slots in when TrimGalore-side gate conditions in §8 are met |

## 5 · Shared assumptions

- **Single static binary invariant** holds across both tools — every new format backend must be pure Rust, no external runtime deps. Rules out `rust-htslib` (libhts), Java-based readers, anything requiring system libraries. The existing `noodles` umbrella (used in TrimGalore via #317, pulled transitively into Bismark already for post-alignment tools) covers uBAM. BINSEQ and mim have pure-Rust crates from their respective Arc/COMBINE-lab authors.
- **Byte-identity baselines stay where they are.** TrimGalore's FASTQ output must continue to match Perl 0.6.11 for the existing flag matrix. Bismark's alignment output must continue to match v0.25.1. New format backends are **additive** — they introduce new output paths, never alter the FASTQ/SAM/BAM baselines.
- **Cargo features for slim builds.** Each non-default backend gets its own feature flag. Prebuilt-release binaries enable all; source builds for size-constrained environments can opt out.
- **Cost floor doesn't scale with LLM speed.** The "agentic AI lowers the addition cost" framing (per Phil's argument in the issues) shifts the *implementation* cost downward but not the *correctness* cost — byte-identity invariants, CI coverage, fixture maintenance, edge-case handling. The right architectural response is one abstraction with N backends (amortising the floor) rather than N feature flags (multiplying it). This epic is structured around that response.

## 6 · Integration points

- **`RecordSource` trait** (input) and **`RecordSink` trait** (output) are the load-bearing seams. Phase 1 defines them in TrimGalore (RecordSource already exists post-#317); Phase 2 adopts the same shape in Bismark; Phase 3 decides whether to formalise as a shared crate.
- **Cargo feature flags** are the user-facing dial. Naming convention: `--features bam` (already pulled via noodles), `--features binseq`, `--features mim`. Default-on for prebuilt binaries.
- **Format-detection layer.** TrimGalore's `src/format.rs` already does content-based detection (peek magic byte → decompress + probe for `BAM\1`). Phase 1 extends it for BINSEQ magic + mim sidecar presence. Phase 2 ports the same detection logic to Bismark (or imports it from a shared crate).
- **CI matrix.** Each repo's `validation` job grows a per-format axis. Cross-tool handshake validation lives in a new dedicated job (likely in TrimGalore's CI, since it's the producer and has the matrix already) or in nf-core/methylseq's integration tests if Phil wants to host them externally.
- **Documentation.** Each tool's README gets a "supported input/output formats" table. CLAUDE.md in each repo gets a "format backends are pluggable via Cargo features" note. The CHANGELOG entries reference each other so users tracking one tool see the cross-tool angle.

## 7 · Sequencing notes

- **Phase 1 already has its first input backend done.** uBAM input landed in TrimGalore#317; the rest of Phase 1 (BINSEQ input, mim input, output-side abstractions) builds on the `RecordSource` trait that PR introduced.
- **Phase 2 is gated on Phase 1's trait stability**, not on every Phase-1 backend landing. As soon as the `RecordSource`/`RecordSink` shape is locked, Bismark can begin consuming it (or mirroring it).
- **Phase 3 is the integration moment.** It's deliberately small — it consolidates rather than introducing new functionality. Its main deliverable is the cross-tool handshake test matrix + the shared-crate decision.
- **Realistic timeline:** Phase 1 (BINSEQ in / mim in / outputs) ~2 weeks of focused work; Phase 2 (Bismark inputs) ~2 weeks; Phase 3 (integration + docs) ~1 week. Coordinated release window mid-July if started now.

## 8 · Phase 4 gating conditions (objective, TrimGalore-side only)

Phase 4 (BINSEQ + mim) on the **TrimGalore side** starts ONLY when ALL of the following are true. Until then, Phase 4 stays ⏸️ Deferred on TG. Re-evaluation cadence: every 6 months from the Phase 1 ship date.

> **Note on cross-repo asymmetry.** Bismark made its own call on BINSEQ — it adopted BINSEQ input via [Bismark#1029](https://github.com/FelixKrueger/Bismark/pull/1029) (2026-06-26). The gates below are derived from invariants specific to TrimGalore: the single-static-binary commitment, the reproducibility CI job (`SOURCE_DATE_EPOCH` double-build bit-identity), and the no-system-library policy. Bismark doesn't carry those invariants, so the same upstream dep can pass Bismark's bar and fail TrimGalore's. Both decisions are correct for their respective repos. See also the [`feedback_cross_repo_epic_hygiene`](../../../.claude) memory for the general pattern.
>
> **Revision note (v5, 2026-06-30).** Audited by @ewels in [#315](https://github.com/FelixKrueger/TrimGalore/issues/315). Two original gates were dropped as proxies-not-substance:
> - "binseq must be 1.0+" — replaced. Rust crates routinely ship production-quality at 0.x; what matters is API + format stability over a window, not a version-string milestone.
> - "mim-index MSRV ≤ TG's 1.88 MSRV" — dropped entirely. TrimGalore is consumed as a binary (`cargo install` / bioconda prebuilt / Docker `:latest` / GitHub Release binaries); MSRV bumps don't break user installations, only source-builders, and source-builders can `rustup update`. MSRV is acknowledged but not blocking.

| Condition | Specific criterion | Current state (2026-06-30) |
|---|---|---|
| BINSEQ format + API stability | `binseq`'s on-disk format AND public API have been stable for 90+ days from the candidate pin point: no format revs, no breaking API changes during the window | 15 versions in 9 months including 2 breaking format revs since 2025-Q4; the stability window has not yet closed. **NOT MET.** |
| BINSEQ reproducibility | `binseq` either drops the `zstd-sys` C build dependency OR adopts a pure-Rust zstd backend by default — required by TG's single-static-binary invariant + reproducibility CI | `zstd-sys` is mandatory in the current dep tree. **NOT MET.** *Bismark accepted this trade-off via Bismark#1029 because Bismark doesn't have the same reproducibility-CI invariant.* |
| mim format spec | Public format spec exists answering at minimum: are byte offsets in compressed-stream or decompressed-stream coordinates? Per-record or per-block? What is the chunk-boundary contract? | No public spec; only Rust implementation behaviour. **NOT MET.** |

If after 12 months from Phase 1 ship date no progress is observed on the remaining gates, the epic explicitly closes Phase 4 as "won't do" and the BINSEQ + mim conversation is reset.

## 9 · Open questions (deliberately left for the phase plans)

- **Shared crate vs. mirrored trait?** Phase 3 spike decides. Shared crate is cleaner long-term but adds release-coordination overhead between the two tools. The `RecordSource` trait is the load-bearing seam (used by both repos for uBAM input). Output-side abstractions remain per-tool — TrimGalore's uBAM-out + Bismark's alignment-out are different shapes and don't share a trait.

---

## Revision history

- **v1, 2026-06-25** — initial epic, three phases. Spawned from cross-repo discussion in TrimGalore#315 + Bismark#1025. uBAM input for TrimGalore (#317) already landed and counts as Phase 1's first deliverable.
- **v2, 2026-06-25** — Phase 1 re-scoped (BINSEQ + mim deferred to a new Phase 4) after dual plan-review on Phase 1 PLAN v1 surfaced 5 critical findings. Phase 4 added.
- **v3, 2026-06-26** — Phase 2 description corrected: it inherits `RecordSource` from #317 (NOT Phase 1's "trait stability" — Phase 1 v2.1 has no new trait). §8 Phase 4 gating conditions made objective (specific version + spec criteria, 12-month sunset). Original §8 "Open questions" renumbered to §9.
- **v4, 2026-06-27** — Phase 2 marked DONE: Bismark uBAM input shipped via Bismark#1026 + #1027 + #1028 (merged 2026-06-25). Phase 1 also marked DONE: TrimGalore#320 + #322 fix merged into `dev` 2026-06-26. Phase 4 row updated to surface that the Bismark side has BINSEQ INPUT merged via Bismark#1029 (2026-06-26), even though TrimGalore's §8 gates still hold for the TrimGalore-side adoption — each repo has its own decision-making bar. **Phase 3 (cross-tool integration) is now unblocked** since both producer (#320) and consumer (Bismark#1027) halves exist. Bismark uses its own phase numbering ("#1025 Phase 1" = uBAM, "#1025 Phase 2" = BINSEQ) which doesn't match this EPIC's numbering — noted inline in Phase 2.
- **v5, 2026-06-30** — §8 gates audited by @ewels via [#315](https://github.com/FelixKrueger/TrimGalore/issues/315) comment. Two of the original v3 gates were proxies, not substance, and have been removed: the "binseq must be 1.0+" criterion (Rust crates routinely ship production-quality at 0.x; restated as a 90-day format+API stability window) and the "mim-index MSRV ≤ TG's 1.88 MSRV" criterion (dropped entirely — TG is a binary-distributed tool, MSRV bumps don't affect user installations). Two load-bearing gates remain: zstd-sys reproducibility (TG-specific single-static-binary invariant — Bismark accepted this trade-off via #1029, TG can't) and mim format spec absence. Asymmetric-cross-repo-invariants framing added inline to §8 head note.
