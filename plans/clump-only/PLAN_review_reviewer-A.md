# Plan review — `--clump_only` (Reviewer A)

## 1. Verdict

**APPROVE WITH REVISIONS.** Direction is sound and the parallel-clump machinery is reused correctly. Blocking issues are (a) the byte-identity contract as written is unenforceable against the current reader/writer, and (b) the rejection matrix has gaps that let contract-violating flags slip through.

## 2. Critical findings

### C1 — Plus-line is discarded on read, emitted as bare `+` on write. Contract in §Behavior is unenforceable.
The plan promises "Plus-line separator preserved as in input" (§Contract) and "byte-identical" record output. Reality:
- `src/fastq.rs:407-411` — `next_record()` explicitly reads and *discards* line 3.
- `src/fastq.rs:66` — `write_to()` hard-codes `b"+\n"`.
- `FastqRecord` (`src/fastq.rs:32-39`) has no field for the plus-line.

Any input whose plus line carries a repeat-header (`+SRR…`) will emit `+` on output — headers, seq, qual are byte-exact but plus-line is not.
**Fix:** either (a) narrow the contract to header+seq+qual byte-identity and explicitly document plus-line normalisation to bare `+`, or (b) add a `plus_tail: Option<String>` to `FastqRecord` and honour it in `write_to`. The plan's §Assumptions bullet "verify" is not sufficient — this is code that must change to match the stated contract.

### C2 — CI byte-identity test cannot catch C1 with the chosen fixture.
`test_files/BS-seq_10K_R1.fastq.gz` uses bare `+` on line 3 (verified by decompressing). The `zcat|sort|diff` CI step (§Implementation outline §9) will therefore pass even if the plus-line is silently mangled.
**Fix:** either add a fixture with `+<repeat-header>` plus-lines to the CI byte-identity step, or (if C1 is resolved by narrowing the contract) explicitly state the CI test's semantics as "header+seq+qual line preservation".

### C3 — Rejection matrix omits three flags that violate byte-identity.
The matrix (§Rejection matrix) misses:
- **`--rename`** (`src/cli.rs:257`) — mutates record ids (`FastqRecord::append_to_id`, `src/fastq.rs:148`). Trivially violates header byte-identity if allowed through.
- **`--discard_untrimmed`** (`src/cli.rs:335`) — filter that drops records. Depending on where it's honoured, this could produce silent record loss (violates the multiset-preservation invariant).
- **`--consider_already_trimmed`** (`src/cli.rs:262`) — trim-shape flag, no meaning in clump-only.

**Fix:** add all three to the rejection matrix and to the CLI-validate unit tests (§Validation #5).

### C4 — uBAM input is silently accepted.
Format detection in `main.rs:159-167` runs unconditionally before the specialty dispatch. If a user passes `--clump_only x.bam`, `open_sync_reader` returns a `BamReader`, records flow through with no plus-line concept, and aux-tag propagation depends on `--preserve-tags` — none of which was reviewed in the plan. §Rejection matrix does not name uBAM input; §CHANGELOG says "FASTQ in/out only in v1" but there is no code path enforcing that.
**Fix:** in `Cli::validate` (or the post-detect check in `main.rs:159-198`), reject `--clump_only` when any `input_formats[i] == InputFormat::UnalignedBam`.

### C5 — §Determinism analysis is factually wrong; the "add tiebreaker" alternative is redundant.
The plan (§Determinism, §Implementation outline §4) proposes converting `sort_single_by_key`/`sort_paired_by_key` to stable sort or adding a `(key, input_index)` secondary. Reality (`src/clump.rs:271-313`):
- Sort is *already* stable (`Vec::sort_by`).
- Tiebreaker is *already* content-based: `key → seq → qual → id`.

Cross-run determinism on identical input therefore holds today for `--clumpify` (already), and will hold for `--clump_only` by inheritance. Nothing in `src/clump.rs` needs to change.
**Fix:** delete §Implementation outline step 4 (or reduce it to a two-line "verify existing behaviour" note) and delete the (a)/(b)/perf-caveat prose in §Determinism. The doc-comment on `sort_single_by_key` (lines 267-270) already promises the property the plan is asking for.

### C6 — `--fastqc` on specialty modes is a new integration point, not a copy-paste.
The plan (§SE steps 8, §Integration/Writes) says "same integration point as the trim path". Grep shows `fastqc::run` is called only in `main.rs:882` and `main.rs:1188` — inside the normal-trim paths, after those functions return. Existing specialty modes (`hardtrim5/3`, `clock`, `implicon`, `demux`) have no FastQC hook at all. Adding it for clump-only is a small but real new integration.
**Fix:** explicitly own this as a new hook in the implementation outline (§Implementation outline §2 or a new §12), and confirm the intended UX ("if `--fastqc` is set alongside `--clump_only`, run FastQC on the reordered outputs after clump completes"). Otherwise silently drop the claim.

## 3. Notable but non-blocking

- **N1** `src/io.rs` helper name. Plan cites `trimmed_output_name` at "line 65-95"; actual is `single_end_output_name` at `src/io.rs:72`. Cosmetic drift.
- **N2** `--basename` interaction. Plan's `clump_only_output_name` signature has no `basename` param; the existing PE naming honours it. Decide whether `--basename foo --clump_only` yields `foo_clumped_{1,2}.fq.gz` or is rejected.
- **N3** `--phred64`/`--phred33` are missing from the rejection matrix. Neither mutates bytes today, but they signal an interpretation that clump-only doesn't apply — reject for clarity.
- **N4** Report filename collision claim (§Report shape, "Why not `_trimming_report.txt`"). The scan-glob argument is correct but weak — it hinges on downstream tooling matching `*_trimming_report.*` exactly and not `*_report.*`. Fine to accept, but the rationale is "we chose a distinct suffix", not "we proved no scanner globs `*clumping*`". The claim `Verify no existing pipeline scans *_clumping* at implementation time` (Resolved decision #3) should be a soft yellow flag, not a green tick.
- **N5** Bin flush order at EOF. The plan's "Write bins in bin-index order" (§SE step 6) is imprecise. Actual behaviour (`src/parallel.rs:1163-1167`) is: bins are flushed as they hit budget during streaming (dispatched round-robin to workers) and the *remaining* non-empty bins are flushed in bin-index order at EOF. Cross-run determinism still holds because record→bin dispatch is deterministic. Reword to avoid future confusion.
- **N6** The `read_id_prefix` machinery (`src/fastq.rs:183`) that enforces R1/R2 header sync in the trim path does *not* run in the clumpy readers today (see `src/parallel.rs:802-836`). Clump-only inherits this "no header cross-check" behaviour — worth calling out in §Assumptions.

## 4. What the plan does well

Clear reuse of the existing bin/sort primitives without threading a skip-trim boolean; correct choice to own a sibling `src/clump_only.rs` module rather than bolt a branch onto `read_single_clumpy`. §Validation cleanly maps every claimed invariant to a concrete test, and the multiset/lockstep/determinism triad is the right decomposition. Reorder-report scope is deliberately narrower than `TrimStats` — good taste. Report-filename divergence from `*_trimming_report.*` is defensible and consciously argued.
