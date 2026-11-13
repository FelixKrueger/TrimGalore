# PLAN_REVIEW_B — Replace `bgzf::MultithreadedWriter` on the parallel uBAM output path (#542)

**Reviewer:** B (independent; no contact with Reviewer A)
**Plan reviewed:** `plans/11132026_issue-542-bgzf-writer/PLAN.md`
**Tree:** `dev` @ `d862f7d`, clean (`git diff --stat` empty before and after this review)
**Method:** noodles-bgzf-0.35.0 and noodles-bam-0.89.0 read from
`~/.cargo/registry/src/index.crates.io-*/`; nine behavioural probes built against the pinned
`noodles-bgzf = "=0.35.0"` in a throwaway crate under the scratchpad; real BAMs reframed and run
through `samtools` 1.x and `bgzip`. No in-tree file was touched. `cargo test` was **not** run (the
other reviewer holds the `target/` lock; the baseline is CI-green `dev`).

**Verdict: safe to implement, but not as written.** The design is sound and I could not fault the
mechanism. Two things must be settled before coding, one of which changes the §4 API shape.

---

## 1. §3 re-derived from source — the four facts

I re-derived all four rather than accepting them. **All four hold.** Probe transcript at the end of
this section.

### 3.1 The issue's sketch appends an EOF block to every batch — **HOLDS**

`noodles-bgzf-0.35.0/src/writer.rs:238-246`:

```rust
impl<W: Write> Drop for Writer<W> {
    fn drop(&mut self) {
        if self.inner.is_some() { let _ = self.try_finish(); }
    }
}
```

`try_finish` (`:189-198`) does `self.flush()?` then `inner.write_all(BGZF_EOF)`. `into_inner`
(`:118-120`) is `self.inner.take().unwrap()` — it leaves `None`, so the `Drop` that runs at the end
of `into_inner`'s own body is a no-op. The plan's reading is exactly right, and
**`into_inner()` is load-bearing, not stylistic** is the correct characterisation.

So the comment in #542's sketch (`// buf holds complete BGZF blocks, no EOF marker`) is false. My
probe, with a 26-byte payload: scoped `write_all`+`flush` then dropped → 82 bytes ending in the EOF
block; `flush`+`into_inner` → 54 bytes, no EOF. (The plan's 66/38 are the same shapes with a
shorter payload — the 28-byte delta is identical, which is the number that matters.)

### 3.2 The terminator needs no private constant — **HOLDS**

A fresh `bgzf::Writer` over an empty `Vec`, `try_finish()`, `into_inner()` → **exactly 28 bytes,
byte-equal to `BGZF_EOF`**. `flush()` short-circuits on an empty staging buffer
(`writer.rs:263-269`), so nothing but the terminator is emitted. The plan's choice over a copied
constant is correct. (See I4 — the plan then contradicts this principle in its own test.)

### 3.3 Compression level is already equal on both sides — **HOLDS, nothing to measure**

`writer/builder.rs:5-9` and `multithreaded_writer/builder.rs:77-83` both default
`compression_level: CompressionLevel::default()`, and `compression_level.rs:93-97` is `Self(6)`.
Stronger than the plan says: both paths funnel into the *same* `deflate::encode` with the same
feature-selected `CompressionLevelImpl` (`libdeflater::CompressionLvl` under the `libdeflate`
feature, else `flate2::Compression`), so the backend cannot diverge between them within one build
either. #529 sets only `set_worker_count`. Level 6 before and after on all three paths. The issue's
"compare rather than assume" is answered by reading; no measurement needed.

### 3.4 No byte-identity invariant on the `.bam` container — **HOLDS, and understates the guard**

Verified: `assert_ubam_eq` (`tests/integration_ubam_out.rs:125-153`) compares header-minus-`@PG` and
decoded record tuples; `ubam_out_cores_gt_1_matches_single_core_across_batches` (`:327-362`) compares
records; the CI `--clump_only` determinism steps compare `samtools view | sort | md5sum`
(`ci.yml:1074-1100`); the Perl matrix is FASTQ-only.

**What the plan misses:** the two unit parity tests it treats as needing an "extension" already pin a
*stronger* invariant than §3.4 credits. `src/parallel.rs:2880-2890` and `:2946-2954`:

```rust
// the BGZF framing may differ; the decompressed BAM stream must not
let (a, b) = (inflate(&serial_out)?, inflate(&par_out)?);
assert_eq!(a.len(), b.len(), "decompressed sizes differ");
assert!(a == b, "decompressed BAM streams differ");
```

That comment was written for exactly this change, and it is the real guard against data loss. It also
tells you which failure modes are and are not covered: a lost block changes `a.len()`, so it fails
loudly; a stray EOF block decompresses to zero bytes, so it passes. Naming this in §3.4 would have
corrected the risk ranking in §9 (see C2).

### Probe transcript

| Shape | Result |
| --- | --- |
| scoped `write_all` + `flush`, dropped | 82 bytes, **ends with EOF** |
| `write_all` + `flush` + `into_inner` | 54 bytes, no EOF |
| fresh writer, `try_finish` + `into_inner` | **28 bytes, == `BGZF_EOF`** |
| `try_finish` **without** `into_inner` | **two** EOF blocks |
| `into_inner` **without** `flush` (26-byte payload) | **0 bytes — whole payload lost** |
| empty encoder, `flush` + `into_inner` | 0 bytes (§5 edge case confirmed) |
| empty encoder, dropped | 28 bytes — a stray EOF (why `flush`+`into_inner` beats `try_finish`) |
| header member + 2 batch members + terminator | round-trips through `bgzf::Reader` |
| 200 KB batch, `flush` + `into_inner` | several members, no EOF (§5 edge case confirmed) |

---

## 2. Load-bearing claims in §§4-8 and §10 that §3 never checked

I checked these; the plan did not.

**§5.1 "the header is its own BGZF member" — verified, untested by the plan.** I took a real
1000-record uBAM, stripped its framing, and rebuilt it in the plan's exact shape: header bytes
(175 B) as their own member, records as further members, one terminator. `samtools view -c` → 1000;
`samtools view -H` → 4 lines; `samtools view | md5` identical to the samtools-written original. The
shape works. Nothing in §9 verifies it (see I1).

**§8 "records and order unchanged, so … `samtools` … unaffected" — true, and unguarded.** See I1.

**§7 "CPU: unchanged in total, redistributed" — the weakest claim in the plan.** Total CPU is
unchanged; *available parallelism for the compression stage* is not. See I2 with measurements.

**§7 "Memory: lower" — plausible, unverified, and the stated figure is off.** `into_bgzf` holds the
raw and the compressed buffer at once; the sink's result channel now carries ~15% as many bytes.
Lower is the right direction. But "each worker gains one `Vec` roughly half the encoded batch size"
comes from the issue's `with_capacity(encoded.len() / 2)`; BAM records at level 6 compress to
**0.146-0.150** in my measurement, so that reserve is 3.3× too large (O3).

**§6.3 "Keep the field-order comment — `inner` before `pending` still matters for drop order" —
half right.** The *order* still matters: on NTFS, `remove_file` on an open handle fails, so a
`pending` dropped first would leave a stale `.partial`. The *comment text* ("so the BGZF stream is
closed before the temporary is removed") becomes false, because the field is no longer a BGZF
stream. Keeping it verbatim violates the repo's own "state the current fact" rule (I3).

**§10 "compression level … deliberately left alone" — dismissed without engaging the repo's own
precedent.** `src/fastq.rs:29` is `DEFAULT_GZIP_LEVEL: u32 = 1`, with a comment explaining that
level 6 is what users opt into via `--compression`. The uBAM path sits at 6 and `--compression` is
documented FASTQ-only (`cli.rs:205, 231`). That asymmetry is pre-existing and fine — but it is
exactly what makes "the FASTQ path already does this" a weaker precedent than §2 claims, because
this change puts a 4.5×-more-expensive compressor on the worker critical path (I2).

**§10 "`--cores 1` never reaches this path" — verified.** `cli.cores > 1` gates all three uBAM
dispatch branches (`main.rs:2379, 2553, 2699`). Note the FASTQ branches use
`cli.cores > 1 || cli.clumpify` (`:1536, 1763, 2059`); the uBAM ones do not, so `--clumpify` cannot
sneak in. `clump_only.rs` imports only `BamWriter` (`:44, 698, 943`) — no `ParallelBamSink` — so §2's
out-of-scope claim holds even at the `--cores 2` the CI uses there.

---

## 3. Action items

### Critical

**C1 — §9.3 is not implementable as written, so the defect the issue exists to fix ends up
untested.** Three compounding problems:

1. `unwritable_file` is at `src/fastq.rs:814`, not `:836`.
2. It is defined *inside* `mod tests` (`src/fastq.rs:808`), so neither `src/bam.rs`'s test module nor
   `tests/` can call it. §9.3 says "Integration test" while pointing at a private unit-test helper.
3. Fatally: it returns a read-only `File`, and `ParallelBamSink` has nowhere to put one.
   `create` takes a `Path` and makes its own handle via `PendingOutput::create` (`io.rs:890`). §4's
   API has **no writer-injection point at all**, so "sink over an unwritable target" can only mean a
   read-only *directory* — which makes `create` fail and never reaches `write_compressed`. That tests
   `PendingOutput`, not C-1.

This matters now rather than later because the fix changes §4. Split the type:

```rust
struct BgzfSink<W: Write> { inner: W }          // compress-append-terminate, testable over any W
pub struct ParallelBamSink {                     // = BgzfSink<BufWriter<File>> + PendingOutput
    inner: BgzfSink<BufWriter<File>>,
    pending: Option<crate::io::PendingOutput>,
    header: Header,
}
```

§9.3 then drives `BgzfSink` over a writer that fails after N bytes and asserts `Err` rather than a
panic. The repo already has the harder version of this trick — `std::io::pipe()` with the read end
closed, at `src/fastq.rs:849-872`, written for #434 precisely because read-only files cannot produce
a *post-creation* write failure. Reuse that reasoning. Keep the existing
`ubam_out_cores_gt_1_failure_leaves_existing_output_untouched` (`tests/integration_ubam_out.rs:366+`)
as the publication half; it covers encoder failure, not write failure.

If the split is refused, say so in §9.3 and state plainly that C-1's retirement rests on code reading
alone.

**C2 — §9.1's assertion cannot fail on the mistake the plan calls the easiest to make.**
`slice::ends_with` returns `false` when the needle is longer than the haystack. So "assert the bytes
do **not** end with the 28-byte EOF block" is satisfied by an **empty** buffer — which is exactly what
`into_inner()` without `flush()` produces for any batch under ~64 KiB (probe: 26-byte payload → 0
bytes out). The second clause ("exactly one appears in a finished file") also passes, because
`finish` writes its one terminator regardless. A §9.1 written to the letter would go green on a
`bgzf_members` that discards every record.

The assertion points the same way as the mutation — the pattern Felix's own harness notes flag. Add a
positive round-trip:

```rust
let blocks = enc.into_bgzf()?;               // non-empty fixture
assert!(!blocks.0.ends_with(&eof_bytes()));  // no mid-stream terminator
let mut v = Vec::new();
bgzf::Reader::new(&blocks.0[..]).read_to_end(&mut v)?;
assert_eq!(v, expected_raw);                 // and nothing was dropped
```

Related recalibration the plan should absorb: the loss *is* caught today, by the `inflate` parity
tests (§3.4 above) and by `ubam_out_cores_gt_1_matches_single_core_across_batches` — just not by the
test §9.1 nominates for it. Conversely I measured the stray-EOF defect's actual cost: `bad.bam`
22,556 B vs `good.bam` 22,444 B over five batches — **28 bytes per batch**, ~0.01% on a realistic
1.17 MB batch, and `samtools view -c`, `samtools quickcheck`, `bgzip -dc`, `gzip -dc` and
`bgzf::Reader` all read it correctly (1000/1000 records). So "would have shipped as a silently bloated
file" is right, but the bloat is negligible; §12's framing of it as the headline hazard, ahead of
C-1's untested write path, is inverted relative to impact. Keep §9.1 — but for the right reason, which
is that readers outside htslib and noodles are unverified (I had Java but no htsjdk to test with).

### Important

**I1 — nothing reads the new framing with htslib, and `ParallelBamSink` has no CI coverage at all.**
After this change TrimGalore owns BGZF framing instead of delegating it to a well-exercised upstream
writer, so this is the moment the guard is needed. Today:

- every automated check on the parallel path reads with **noodles** (`assert_ubam_eq`, `bam_tuples`,
  the `inflate` parity tests);
- every `--cores` + `--output-format ubam` line in `ci.yml` (`:1038, 1058, 1078, 1080, 1095`) is
  `--clump_only`, i.e. the serial `BamWriter`;
- every trim-path uBAM step (`:1266, 1284, 1309, 1328, 1355, 1373`) omits `--cores`, i.e. also serial.

So `ParallelBamSink` is exercised only by `cargo test`, and never by htslib. I verified by hand that
the plan's framing is fine (§2 above), so this is an unguarded invariant rather than a predicted
failure — but it is cheap to close: one step in the job that already installs samtools at
`ci.yml:509-513`, running `--output-format ubam --cores 4` and asserting `samtools view -c` plus
`samtools quickcheck`. Five lines for the only structural check on a surface the project now owns.

**I2 — §7's parallelism claim is wrong in direction, and §9.5 has no decision rule.** Measured on
this machine, a 1.17 MB batch (the shape of 4096×150 bp) through `bgzf::Writer`:

| level | ms/batch | MB/s/core | compressed ratio |
| --- | --- | --- | --- |
| 1 | 3.1 | 384 | 0.150 |
| 6 | 14.1 | 85 | 0.146 |
| 9 | 20.6 | 58 | 0.145 |

(Input tiled from a real uBAM, so the ratios are optimistic in absolute terms; the *relative* figures
are what matter.)

14 ms of compression per batch is the same order as the trimming work it is being folded into. Before,
2N threads shared trim and compress; after, N threads do both. "CPU unchanged in total, redistributed"
is true of total CPU and false of available parallelism — on a host with ≥2N idle cores this is a real
throughput loss, plausibly approaching 1.5×. The saving grace, which the plan should state because it
is the actual argument for the change: when `--cores` matches the cgroup allocation — the nf-core
`--cores task.cpus` contract the issue itself cites — the old design was *oversubscribing*, so there
is nothing to lose and the thread-budget fix is a pure win. The plan's inference from "the speed-up
plateaued at `--cores 4`" is not evidence either way; a plateau is equally consistent with the writer
side binding.

Two concrete asks. Give §9.5 a threshold and an abort condition, not just the expectation "within
noise, or faster" — the plan currently has no stated response to a 30% regression, and §9.4/§9.5 are
described in §12 as things that "would be a follow-up, not a reason to hold this". And name the
fallback: `writer::Builder::set_compression_level(CompressionLevel::FAST)` recovers 4.5× of
compression CPU for ~3% more bytes, and would align the uBAM path with the repo's own
`DEFAULT_GZIP_LEVEL = 1`. I am **not** asking for the level change in this PR — it alters output size
and is a product decision — only that §10 record it as the named mitigation instead of dismissing the
level question as "no case raised".

**I3 — four stale doc comments the plan does not list, plus one it says to keep that becomes false.**
§2's placement table names only `main.rs`'s thread comment. Also stale after the change:

| Location | Text that stops being true |
| --- | --- |
| `src/bam.rs:1109-1110` | "compressed to BGZF by a pool of threads (`noodles` multithreaded writer)" |
| `src/bam.rs:1112-1113` | "so the BGZF stream is closed before the temporary is removed" (field order still matters on NTFS; reword to the drop-order fact, don't keep the text) |
| `src/bam.rs:1120` | "`workers` is the number of compression threads" |
| `src/parallel.rs:352-355` | "the sink compresses them to BGZF on its own pool of `cores` threads" |
| `src/parallel.rs:459-461` | "the sink compresses on `cores` threads" |

A doc comment describing a thread pool that no longer exists is the precise failure the repo's
comment rule targets. §6.3's "Keep the field-order comment" should read "reword it".

Judged against that same rule, §6.1's own guidance ("Comment, one line, on the invariant:
`into_inner` disarms the EOF-appending `Drop`") is **correct** — one line, current fact, no
narration of the old code, and sited on the helper it constrains. No change needed there.

**I4 — §9.1 needs a hardcoded 28-byte constant, which §3.2 argues against.** `BGZF_EOF` is
`pub(crate)` *in noodles*, so unreachable from TrimGalore in test code too. Copying the 28 bytes into
a test creates the second source of truth §3.2 refused. Derive it instead, which my probe confirms
works: a fresh `bgzf::Writer` over an empty `Vec`, `try_finish()`, `into_inner()` → exactly those 28
bytes. One helper in the test module, and noodles stays authoritative.

### Optional

- **O1** — §4 declares `pub fn finish(self)`; it needs `mut self` (it takes `pending` and flushes).
- **O2** — line-number drift that will cost the implementer a grep: the stale comment is
  `main.rs:2700` (§2 says `:2703` — that is the `create` call; §6.5 says `:2702`); `unwritable_file`
  is `fastq.rs:814`, not `:836`; the `clump_only` uBAM CI block starts at `ci.yml:1031`, not `:1075`.
- **O3** — drop the `with_capacity(encoded.len() / 2)` from the issue's sketch. Measured ratio is
  0.146-0.150, so `/ 2` over-reserves 3.3×. `/ 8`, or plain `Vec::new()` and let it grow.
- **O4** — §4 gives `BgzfBlocks` a `len`/`is_empty` pair with no stated caller. The §9.1 test can
  reach `.0` from the same module. Either name the caller or drop both (the `pub` newtype with no
  public constructor is the good part of that design, and it survives without them).
- **O5** — the whole point of the issue is that an `ENOSPC` finally gets reported, so report it with
  the path. `write_encoded`'s current context is `"failed to write BAM records"` and `finish`'s is
  `"failed to finalise uBAM output (BGZF EOF marker)"`; neither names the file.
  `PendingOutput::tmp_path()` (`io.rs:921`) exists for exactly this.
- **O6** — §2 says "Three files, plus tests"; §6.7 also touches `CHANGELOG.md`.
- **O7** — §8 calls the signature change a breaking library change but does not say what version bump
  it implies for a crates.io-published crate. One clause.
- **O8** — pre-existing oddity worth knowing before writing §9.1's "exactly one" assertion: the
  **serial** `BamWriter` emits **two** EOF blocks. `bam::io::Writer::try_finish` forwards to
  `bgzf::Writer::try_finish` (`noodles-bam-0.89.0/src/io/writer.rs:161-163`) without taking `inner`,
  so the subsequent `drop(self)` appends a second one. Confirmed on disk:
  `test_files/ubam_out_se_REFERENCE.bam` contains 2, samtools' own `test_files/ubam_test.bam`
  contains 1. Harmless — samtools reads both — but it means the new parallel output (one terminator)
  will differ from the serial reference in block count, and that a stray EOF block has been shipping
  from the serial path, unnoticed by every test and tool in this repo. Assert "exactly one" only
  against the parallel output, and consider a one-line follow-up issue for the serial double.

---

## 4. Alternatives considered

- **Keep a small dedicated compressor pool (e.g. 2 threads, independent of `--cores`).** Retires
  C-1 and H-2 the same way — the writes become fallible on a thread we own — while keeping compression
  off the trim workers' critical path. Costs a channel and an ordering guarantee the plan currently
  gets for free from the reorder buffer. Worth naming in §9.5's fallback alongside the level change,
  not worth preferring up front.
- **Lower the BGZF level rather than restructure the threading.** Would recover throughput without
  touching the thread model, but retires none of the three defects. Not an alternative to this change,
  only to its performance cost.
- **A deprecated `write_encoded` shim delegating to `write_compressed`.** Would soften O7 at the cost
  of keeping the `&[u8]` boundary the `BgzfBlocks` newtype exists to close. The plan is right to
  rename outright; the newtype is the best idea in §12 and should not be diluted.

---

## 5. Bottom line

The mechanism is right and I could not break it: §3's four facts all hold under independent
re-derivation, the framing reads correctly in samtools, and `BgzfBlocks` closes the silent-corruption
boundary properly. Settle **C1** before writing code, because it changes §4's shape, and fix **C2**
before writing the test it describes. **I1** and **I2** are what I would insist on at review time:
one is five lines of CI for a surface the project is about to start owning, the other is a measured
risk the plan currently argues away rather than gates.
