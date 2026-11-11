# Draft issue — validate the `--phred64` claim against the input

**Status:** drafted for review, NOT filed.
**Suggested labels:** check they exist on the repo before using any.

---

## Title

`--phred64` on Phred+33 input is only caught by FastQC, and only after output is written

## Body

`--phred64` is a claim about the input that nothing checks. When the claim is wrong the
consequences differ by path, and none of them is a clean refusal:

- **Trimming, no `--fastqc`.** Every score drops by 31 (floored at 0), turning typical Q2–Q41 into
  Q0–Q10. At the default `-q 20` this usually empties the file; at lower cutoffs it produces
  plausible-looking poor-quality data. The startup NOTE says this is "easy to miss downstream",
  which is accurate, and a NOTE is all there is.
- **`--clump_only`, no `--fastqc`.** No quality arithmetic happens, so the output is byte-identical
  to the input and the flag is genuinely inert.
- **Anything with `--fastqc`.** Since v2.4.0 the bundled engine is told the encoding, so FastQC
  refuses a file holding any quality character below ASCII 64 — correct, but it fails *after* the
  trimmed or clumped output and the trimming reports are already on disk, and it exits 1. A pipeline
  gating on exit status reads a successful trim as a failed step.

The check is cheap and unambiguous: a Phred+64 file cannot contain a quality character below ASCII
64, because Q0 *is* ASCII 64. So the input alone decides it, and Trim Galore already streams the
input for adapter auto-detection.

### Proposal

Track the lowest quality character during the scan Trim Galore already performs, and if `--phred64`
is set and that character is below ASCII 64, refuse before writing anything — naming the flag, the
offending character, and the likely cause.

```
Error: --phred64 declares Phred+64 input, but this file contains quality
character '?' (ASCII 63). A Phred+64 file cannot hold anything below
ASCII 64, so the input is most likely Phred+33. Drop --phred64.
```

### Why this is not a v2.4.0 change

It changes behaviour on a documented path that has nothing to do with FastQC: `--phred64` without
`--fastqc` currently runs to completion with degraded output, and this would make it a hard refusal.
That is a better default, but it is a behaviour change needing its own tests and release note rather
than something to slip in beside a dependency bump.

### Open questions

1. **Refuse or warn?** Refusing is right for the trimming path. `--clump_only` without `--fastqc` is
   genuinely lossless under the flag, so a refusal there would reject a working invocation — the
   same tension `src/main.rs:332` already records for the uniform BAM-input rejection, resolved the
   same way (uniform beats mode-dependent) unless there is reason to differ.
2. **How much input to scan?** Adapter detection reads the first 1M reads. A low character later in
   the file would still reach FastQC and fail there, so this reduces the window rather than closing
   it. Worth deciding whether that is acceptable or whether the check belongs in the per-record path.
3. **Should the FastQC step still abort?** Independent of this, there is a case that a post-trimming
   QC failure should not fail a run whose trimming succeeded. Changing it would affect every FastQC
   failure mode, not just this one.

### References

- `src/main.rs:332` — the uniform BAM-input rejection and its mode-dependent-rule reasoning.
- `src/main.rs:~487` — the `--phred64` startup NOTE, hoisted in v2.4.0 so it fires on every path.
- `src/fastqc.rs` — where `--phred64` reaches `FastQCConfig`.
- `fastqc-rust` `src/utils/phred.rs:56` — `detect()`, which returns `Err` below ASCII 64.
- `tests/integration_fastqc_phred64.rs` — the current coverage.
