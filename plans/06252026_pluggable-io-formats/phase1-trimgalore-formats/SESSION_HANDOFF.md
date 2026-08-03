# Session handoff — Phase 1 implementation, Step 2 entry

This is the prompt to paste into a fresh `/clear`-ed Claude Code session to resume implementation of the uBAM-output feature. It assumes the session is launched in the TrimGalore repo at `/Users/fkrueger/Github/TrimGalore`.

---

## Paste this into the new session

I'm continuing implementation of the uBAM-output feature for TrimGalore (Phase 1 of the pluggable-I/O-formats epic). The prior session paused mid-Step-2 after hitting a cluster of noodles 0.88 API specifics; you're picking up from a clean checkpoint with Steps 0-1 done.

**MANDATORY reading order before touching code:**

1. `/Users/fkrueger/.claude/CLAUDE.md` — global rules. Especially: the exact-trigger requirement for implementation (already satisfied; the prior session received "implement"), the dual-reviewer convention, "make thinking visible / no silent skips" discipline.
2. `/Users/fkrueger/Github/TrimGalore/CLAUDE.md` — project conventions. Especially: "no external runtime deps" invariant, byte-identity-to-Perl-0.6.11 for FASTQ output, `cargo test --release` from crate root.
3. `~/.claude/skills/code-implementation/SKILL.md` — implementation discipline (you ARE the code-implementation agent for this session).
4. `plans/06252026_pluggable-io-formats/EPIC.md` v3 — epic + Phase 4 deferred-formats gate criteria.
5. **`plans/06252026_pluggable-io-formats/phase1-trimgalore-formats/PLAN.md` v2.1 — the spec you implement against.** Read §0 (revision history at a glance), §3 (behaviour), §4 (signatures), §5 (implementation outline — your task list), §9 (validation).
6. `plans/06252026_pluggable-io-formats/phase1-trimgalore-formats/PLAN_REVIEW_{A,B,C,D}.md` — the four prior reviewer reports, for context on which decisions are locked vs. open.

**What's already done (do NOT re-implement):**

- **Step 0** — `write_paired_reports` helper extracted from `run_paired` + `run_paired_ubam_single_file`. Both call sites refactored to call the helper. Pure refactor, 296 → 296 tests pass.
- **Step 1** — `--output-format <fastq|ubam>` CLI flag in `src/cli.rs`:
  - `OutputFormat` enum (with `clap::ValueEnum` derive) declared ABOVE the `Cli` struct (gotcha: original draft put it between `Cli`'s `#[derive(Parser,Debug)]` line and the struct, which hijacked the derives — that's now fixed).
  - `Cli::output_format` field with `default_value_t = OutputFormat::Fastq`.
  - 5 PLAN §3.4a CLI-level rejection rules in `Cli::validate()`: ubam + clumpify / passthrough / clock / implicon / demux. Each has a unit test.
  - PLAN §3.4b format-detection-time rule in `src/main.rs`: `--preserve-tags` + all-FASTQ-inputs + `--output-format ubam` → hard error (was just a warning when output is FASTQ).
- 305 tests pass (286 lib + 9 new cli + 2 + 8 integration). fmt clean, clippy clean.

**What's pending (your task):**

Steps 2-7 per PLAN §5. The natural next move is Step 2 (`BamWriter`).

**The 5 noodles 0.88 API specifics that blocked the first-pass BamWriter draft** (see TODO comment in `src/bam.rs` where the stub is; the draft was rolled back to keep the checkpoint clean):

1. **`sam::Header::builder().set_header(...)` signature.** The first-pass tried `set_header(Header::default())` (where `Header` is `noodles::sam::header::record::value::map::Header`) but the actual signature wants `Map<Header>`. Investigate: `~/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/noodles-sam-0.69.0/src/header/builder.rs`.
2. **`sam::header::record::value::map::tag::program` path** doesn't exist. The first-pass tried `tag::program::VERSION` and `tag::program::COMMAND_LINE` constants. Find the actual paths — likely something like `map::program::tag::VERSION` or similar. Source under `~/.cargo/registry/src/.../noodles-sam-0.69.0/src/header/record/value/map/program/`.
3. **`record_buf::Name`** doesn't exist at the path `noodles::sam::alignment::record_buf::Name`. The first-pass tried `Name::from(name_bytes.into())`. Find the actual name type — may be `BString`, or `&[u8]` directly assigned, or a `Name` re-exported from a different module.
4. **`write_alignment_record` is a trait method.** Needs `use noodles::sam::alignment::io::Write;` to be in scope. Source: `~/.cargo/registry/src/.../noodles-sam-0.69.0/src/alignment/io/write.rs:13`.
5. **One unresolved type mismatch** in record-construction — likely related to one of the above once they're fixed. Run `cargo build --release` after each fix to see the next error.

**Approach for Step 2:**

Don't fight noodles' API in the dark. Before writing more code, briefly probe the actual paths via `grep -rn "pub struct\|pub fn\|pub use" ~/.cargo/registry/src/index.crates.io-*/noodles-sam-0.69.0/src/alignment/record_buf/ ~/.cargo/registry/src/index.crates.io-*/noodles-bam-0.73.0/src/` for the specific items. Then write `BamWriter::create` + `write_record` + `finish` per PLAN §4. Then write the `split_name_and_tag_tail` + `parse_tag_tail` helpers per PLAN §3.6 (with A/Z/i/f only; reject B/H arrays). Then unit tests per §5 step 2.5.

**Per-step verification gate is non-negotiable:**

After each PLAN §5 step, run:

```sh
cargo build --release
cargo test --release          # must be ≥ 305 passing (305 was the Step 1 baseline)
cargo fmt --all -- --check
cargo clippy --release --all-targets -- -D warnings
```

If any fails, fix before moving on. Do NOT silently skip a step or move on with a failing gate. The prior session's Step 0 + Step 1 are committable mini-milestones; same discipline applies to Steps 2-7.

**Locked decisions — do NOT relitigate:**

- uBAM output is opt-in via `--output-format ubam` (resolved in PLAN v1).
- Paired uBAM output is ONE interleaved BAM (PLAN v2 — matches samtools/Picard/fgbio).
- No `RecordSink` trait (PLAN v2 — FASTQ and uBAM live on different code paths).
- Aux tag round-trip: parse from `FastqRecord.id` tail per PLAN §3.6 (PLAN v2.1; option (a)).
- BINSEQ + mim deferred to Phase 4 (PLAN v2 + EPIC §8).
- Per-record `is_unmapped()` check on EVERY BAM record (post-#317).

**Critical reminders:**

- `--output-format ubam` is currently parseable + validated but `main.rs` dispatch on it is NOT yet wired. A user invoking it today gets FASTQ output silently. Step 4 wires the dispatch (calling Step 3's trimmer entry points which take a `BamWriter` from Step 2). DON'T add the dispatch before Step 2 + Step 3 exist or the binary will fail to build.
- Don't auto-commit. Wait for explicit user direction.
- After all 7 steps land, the workflow says to run `/code-reviewer` (dual reviewers, fresh contexts via Agent tool) and `/plan-manager` to audit coverage. Don't skip this — it's caught real bugs at every prior milestone in this feature.

Begin by reading the artifacts above, then resume Step 2.

---

## After-pasting

Save this prompt as your first message in the `/clear`-ed session. The new session has zero context from this one — everything it needs is in the prompt + the files it references.
