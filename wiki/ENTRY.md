# Trim Galore wiki

Trim Galore is fast adapter and quality trimming for NGS data, with specialty
modes, built-in QC, and a compression and archival mode. It ships as one
binary with no external runtime dependencies, and is used mostly from
pipelines, plus directly at the command line.

Start here. Find the row that matches the task, read those pages, and stop.
`index.md` lists every page.

| The task mentions | Read these pages, and stop |
|---|---|
| Whether to use this tool at all, a comparison with Cutadapt, what it is for | `overview.md` |
| Installing, bioconda, conda, Docker, a container, crates.io, Homebrew | `usage/install.md` |
| A first run, trimming a sample, single-end, paired-end, multiple pairs, auto-detection | `usage/trim-reads.md`, `reference/cli-flags.md` |
| An adapter sequence, `-a`, a preset, Nextera, smallRNA, a FASTA of adapters | `usage/trim-reads.md`, `reference/cli-flags.md` |
| RRBS, MspI, non-directional, end-repair, methylation bias | `usage/rrbs.md` |
| clumpify, compression, archival, output size, `--memory`, a memory floor | `usage/clumpify.md`, `internals/parallel-and-clumping.md` |
| uBAM, a BAM input or output, aux tags, `--preserve-tags`, an interleaved file | `usage/ubam.md`, `reference/inputs.md` |
| What files it accepts, a format, `.bgz`, a pipe, process substitution, stdin, a samplesheet | `reference/inputs.md`, `troubleshooting.md` |
| An output filename, `_val_1`, `_trimmed`, a trimming report, JSON, MultiQC | `reference/outputs.md` |
| A specific flag, a default value, specialty modes, hardtrim, clock, implicon, demux, passthrough | `reference/cli-flags.md` |
| A flag that does nothing, `--gzip`, `--path_to_cutadapt`, an old v0.6 command line | `reference/legacy-flags.md` |
| The order of trimming steps, poly-G versus poly-A, where filtering happens | `internals/trimming-pipeline.md` |
| `--cores`, parallelism, speed, why output bytes differ between runs | `internals/parallel-and-clumping.md`, `usage/trim-reads.md` |
| An error, a refusal, a warning, a crash, missing output | `troubleshooting.md` |
| Byte-identity with the Perl version, the v0.6.x relationship | `overview.md` |
