# SPIKE — observed write set per dispatch arm

**Date:** 2026-08-10 · **Binary:** release, rebuilt from `dev` @ `3008b37`
**Purpose:** derive #414's expectation by observation rather than hand-enumeration.

## Method

`probe_arms.sh` and `probe_flag_variants.sh` each stage fixtures in a fresh
per-case directory, snapshot it recursively, run one arm with `cwd` = that
directory, snapshot again, and print the created set. 40 arm/flag combinations.

Reproduce:

```bash
cargo build --release        # target/ was stamped 38deb1e, NOT HEAD — rebuild
zsh probe_arms.sh          ./target/release/trim_galore /tmp/probe1 "$PWD"
zsh probe_flag_variants.sh ./target/release/trim_galore /tmp/probe2 "$PWD"
```

## Findings

1. **Baseline green.** Across all 11 pre-flighted call sites, every created file
   was already in that arm's candidate list. `--demux`'s
   `_demultiplexing_summary.txt` — the likeliest sixth instance — *is* planned
   (`demux.rs:156`).
2. **`--fastqc` is the whole unplanned set,** and it is a rule: per primary `P`,
   `<P-stem>_fastqc.{html,zip}` beside `P`, including for the `--passthrough`
   carrier and for a `.bam` primary. `--fastqc_args "-o DIR"` relocates both.
3. **`--clumpify` leaves no temp files.**
4. **Two arms have no pre-flight** — single-file interleaved paired, FASTQ-out and
   uBAM-out. Both hand-roll `{stem}_R{N}_trimming_report.{txt,json}`.
5. **Those two cannot alias their input**, proven adversarially: a BAM named
   `x_val_1.fq` yields `x_val_1_val_1.fq`; one named `y_R1_trimming_report.txt`
   yields `y_R1_trimming_report_R1_trimming_report.txt`. Input survived both runs.
   Every namer appends to a stem taken from the input's own filename, so a fixed
   point would need an infinitely long name.
6. **`--retain_unpaired` creates both `_unpaired_N` files even when empty** — so
   `planned ⊆ created` must not be asserted.

## Probe bugs worth inheriting as case requirements

- Paired mates must carry **identical read IDs** or `--passthrough` aborts with a
  sync error after writing three files.
- The `--demux` barcode file needs **real tabs**; `print -r --` does not expand
  `\t`.
- `--clumpify` requires `--cores >= 2`.
- `--basename` with multiple single-end inputs is refused at the CLI.

## Raw output

`OBSERVED_write_sets_by_arm.txt`, `OBSERVED_write_sets_by_flag.txt`.
