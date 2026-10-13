# Page index

Every page in this wiki, grouped by part.

## Root

| Page | Description |
|---|---|
| `overview.md` | The purpose, who runs it, and when it is and is not the right tool. |
| `troubleshooting.md` | What Trim Galore refuses, what it only warns about, and the remediation each message points at. |

## usage/

| Page | Description |
|---|---|
| `usage/install.md` | The distribution channels, which one to use in a pipeline, and what the binary needs at runtime. |
| `usage/trim-reads.md` | The ordinary single-end and paired-end runs, what happens without any flags, and the multi-pair case. |
| `usage/rrbs.md` | Running MspI-digested RRBS data, what `--rrbs` changes, and the non-directional variant. |
| `usage/clumpify.md` | Using `--clumpify` while trimming, `--clump_only` for reorder-only archival, and how to pick `--memory`. |
| `usage/ubam.md` | Running unaligned BAM input, writing uBAM output, carrying aux tags, and the v1 limits of the BAM path. |

## reference/

| Page | Description |
|---|---|
| `reference/inputs.md` | How Trim Galore counts and pairs input files, how it decides what each file is, and the auxiliary input files the flags take. |
| `reference/outputs.md` | Where outputs land, how each one is named, whether they are gzipped, and what the trimming reports contain. |
| `reference/cli-flags.md` | Every flag Trim Galore accepts, grouped by what it affects, with defaults. |
| `reference/legacy-flags.md` | The v0.6.x flags Trim Galore still parses, and why each one no longer does anything. |

## internals/

| Page | Description |
|---|---|
| `internals/trimming-pipeline.md` | The per-read step order in `trim_read`, which module owns each step, and where filtering and auto-detection attach. |
| `internals/parallel-and-clumping.md` | How `--cores` distributes trimming and compression, how `--clumpify` reorders reads into bins, and how the memory budget is resolved. |
