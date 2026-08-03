# Plan review — `--clump_only` (Reviewer B)

## 1. Verdict

**APPROVE WITH REVISIONS** — plan is well-structured and reuses primitives correctly, but the load-bearing "byte-identity" invariant contradicts existing reader/writer semantics and must be either softened in wording or backed by a `FastqRecord` change before implementation.

## 2. Critical findings (ranked)

### C1. Plus-line preservation is UNENFORCEABLE with the current reader/writer

- **What's wrong.** The Contract at PLAN §Behavior line 57 asserts *"Plus-line separator preserved as in input"*. This is false by construction. `FastqReader::next_record` at `src/fastq.rs:373-377` and again at `407-411` **reads the plus line into a buffer and discards its contents** — only its presence is verified. `FastqRecord` (`src/fastq.rs:32-39`) has no plus-line field, and `FastqRecord::write_to` at `src/fastq.rs:66` **always** emits the literal bytes `b"+\n"`. Any input with `+<header-repeat>` or comments on line 3 will be normalized to `+` in the output.
- **Evidence.** `src/fastq.rs:373-377` (reader discards) + `src/fastq.rs:66` (writer emits bare `+`). PLAN Assumption 4 (line 345) already flags this as "verify" but the plan's Self-Review (line 392) labels it non-blocking — that's wrong for an invariant sold as the mode's entire point.
- **Fix.** Pick one: (a) soften Contract line 57 to *"header, sequence, and quality bytes byte-exact; plus-line normalized to `+`"* and update issue #353 / CHANGELOG / docs to say so; or (b) add `pub plus_line: String` to `FastqRecord`, populate it in the reader (both `Direct` and `Threaded` paths), and use it in `write_to`. Option (b) is a much bigger change with cross-cutting perf implications for the trim path. Option (a) is the pragmatic call, but requires an explicit decision, not silence.

### C2. CRLF → LF collapse also breaks byte-identity

- **What's wrong.** `src/fastq.rs:364, 371, 384, 398, 405, 418` all call `trim_end_matches(['\n', '\r'])` on the four record lines; the writer then writes only `\n`. Windows-line-ending FASTQ input silently becomes LF-only output. Not byte-identical.
- **Fix.** Either doc it as a limitation in the same wording as C1, or normalize on the way in *and* record the original terminator per record. Practical answer: doc-only, but say it.

### C3. CI byte-identity check is line-multiset, not record-multiset

- **What's wrong.** PLAN §9 step's check `zcat X | sort > a; zcat Y | sort > b; diff -q a b` compares the **multiset of lines**, not of records. A bug that mixed a header line into another record's slot would still pass, because both files would contain the same 40 000 lines (10 000 headers, 10 000 seqs, 10 000 `+`, 10 000 quals). The test cannot catch the failure modes it's supposed to guard.
- **Fix.** `zcat X | paste - - - - | sort > a; zcat Y | paste - - - - | sort > b; diff -q a b`. Groups the 4 lines per record before sorting.

### C4. Rejection matrix cannot reliably reject flags with clap defaults

- **What's wrong.** PLAN §Rejection matrix lists `-q/--quality`, `--stringency`, `-e/--error_rate`. All three are `#[clap(default_value = "...")]` (see `src/cli.rs:48, 92, 96`) with concrete types, not `Option<T>`. `Cli::validate()` cannot tell "user passed `-q 20`" from "user passed nothing". Checking `self.quality != 20` (the existing `--nextseq` conflict pattern at `src/cli.rs:581`) silently accepts explicit `-q 20 --clump_only`, and any check that fires on default 20 rejects every `--clump_only` invocation.
- **Fix.** Either (a) accept these knobs silently on the `--clump_only` path (they're ignored because the trim path never runs — arguably fine, since it can't do harm), or (b) use clap's `ArgMatches::value_source()` to detect explicit user-set values. Pick one and say so.

### C5. `--basename` / `--no_report_file` handling missing

- **What's wrong.** Neither the signature nor the rejection matrix mentions `--basename`. Existing PE naming (`src/io.rs:166-172`) branches on it. If `--clump_only --basename foo` is silently accepted, output naming is undefined. Same for `--no_report_file` — should it suppress `*_clumping_report.txt`?
- **Fix.** Pick a rule per flag and document in the matrix. Cheapest: accept both, thread through the naming helper and the report writer.

### C6. §4 is looking for a problem that doesn't exist

- **What's wrong.** PLAN §Determinism (line 80) tells the implementer to "verify" that the sort primitives are stable, with two fallbacks. `src/clump.rs:275` and `src/clump.rs:302` already use `Vec::sort_by`, which is stable, and the tiebreaker cascades to `(key, seq, qual, id)`. No verification needed, no Option A/B/C branch, no perf risk from switching to stable.
- **Fix.** Delete §4's branch structure; leave a one-liner referencing the existing sort behaviour.

### C7. Report line "Compression ratio" is meaningless with `--dont_gzip`

- **What's wrong.** PLAN §Report shape prints `<input_bytes>/<output_bytes>`. With `--dont_gzip` (compatible per line 114) and a gzip input, ratio is inverted (input compressed, output uncompressed). Confusing at best, misleading in reports.
- **Fix.** Only emit the ratio line when both sides are gzip.

## 3. Notable but non-blocking

- Report byte-count acquisition method not specified (streaming counter vs `fs::metadata` after close). Prefer post-close metadata to avoid double-counting.
- `STATIC_OVERHEAD_BYTES = 512 MiB` at `src/clump.rs:168` was calibrated *including* FastQC + adapter-align memory. Clump-only doesn't use those; the reservation is over-generous, not incorrect. Fine for v1.
- §9 sanity test greps stderr for `"Mode: --clump_only"` — no such line exists today; the plan should specify **which module** prints it and **when** (before or after layout resolution).
- `sort_paired_by_key` tiebreaker only compares R1 fields (`src/clump.rs:302-307`). Stable sort makes this fine, but worth noting explicitly in Assumptions.

## 4. What the plan does well

Reuses the right primitives (`canonical_minimizer`, `bin_for`, `sort_*_by_key`, `resolve_layout`) with the right justification for owning a separate dispatcher rather than threading a "skip trim" boolean through `read_single_clumpy`. Rejection matrix is comprehensive on the multi-flag combinatorics. The `*_clumping_report.txt` filename decision is well-argued against nf-core's `*_trimming_report.*` scan glob (verified — nothing in `src/` builds that filename). CI job structure is right, only the diff mechanic is wrong. Signatures are minimal and typed correctly (`ClumpOnlyStats`, path/gzip/cores/memory/compression). Efficiency claim (skip-trim is a strict win) is sound.
