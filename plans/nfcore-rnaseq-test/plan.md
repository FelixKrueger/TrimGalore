# Plan: Draft PR to nf-core/rnaseq for TrimGalore v2.0 Testing

## Goal

Create a draft PR to nf-core/rnaseq that swaps the TrimGalore container from
`biocontainers/trim-galore:0.6.10` to `ghcr.io/felixkrueger/trimgalore:dev`,
applies necessary tweaks for v2.0 compatibility, and validates via nf-test.

---

## Compatibility Analysis

### What works out of the box

| Feature | Status | Notes |
|---------|--------|-------|
| CLI flags (`--paired`, `--cores`, `--gzip`) | OK | v2.0 accepts all these; `--gzip` is hidden/deprecated (on by default) |
| Output file patterns (`*_trimmed.fq.gz`, `*_val_{1,2}.fq.gz`) | OK | Identical naming |
| Report file naming (`*_trimming_report.txt`) | OK | Same pattern |
| Report parsing regex: `sequences processed in total` | OK | Exact string present in v2.0 reports |
| Report parsing regex: `shorter than the length cutoff` | OK | Exact string present in v2.0 reports |
| Version extraction: `trim_galore --version \| grep -Eo "[0-9]+(\\.[0-9]+)+"` | OK | Outputs `2.0.0` |
| `--cores N` parallel compression | OK | Supported in v2.0 |

### What needs changes

| Issue | Severity | Fix |
|-------|----------|-----|
| `--fastqc_args` in default ext.args | **Breaking** | Container has no FastQC. Remove from preset args; run FastQC as separate process if needed. |
| Snapshot mismatches in nf-tests | **Expected** | Reports have new format (Cutadapt-compat section). Snapshots need updating. |
| trim_html / trim_zip outputs | **Missing** | With `--fastqc_args` removed, no `.html`/`.zip` from TRIMGALORE process. These channels become empty. |

### The `--fastqc_args` problem in detail

**Current default** (subworkflow config line 15):
```groovy
def preset_args_map = ["--fastqc_args": "'-t ${task.cpus}'"]
```

This passes `--fastqc_args '-t N'` to trim_galore, which implies `--fastqc` and
runs FastQC inside the trim_galore process. The old biocontainers image bundled
FastQC + Java + Python. Our container is just the Rust binary.

**Fix**: Remove `--fastqc_args` from the preset args map. The subworkflow
already runs FASTQC as a separate process on raw reads. If trimmed-read FastQC
is needed, it should be a separate FASTQC process call — which is more
nf-core-native anyway (one tool per container).

---

## Deliverables

| # | File | Change |
|---|------|--------|
| 1 | `modules/nf-core/trimgalore/main.nf` | Swap container to `ghcr.io/felixkrueger/trimgalore:dev` |
| 2 | `subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/nextflow.config` | Remove `--fastqc_args` from preset args |
| 3 | nf-test snapshots | Update after running tests |

---

## Detailed Changes

### 1. Container swap (`modules/nf-core/trimgalore/main.nf`)

Replace lines 6-8:
```groovy
// Before:
container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/trim-galore:0.6.10--hdfd78af_2' :
    'biocontainers/trim-galore:0.6.10--hdfd78af_2'}"

// After:
container "ghcr.io/felixkrueger/trimgalore:dev"
```

Single container reference works for both Docker and Singularity (Singularity
pulls Docker images from GHCR via `docker://` protocol automatically).

### 2. Remove `--fastqc_args` preset (`subworkflows/.../nextflow.config`)

Change the TRIMGALORE ext.args closure to remove the fastqc_args preset:
```groovy
// Before:
def preset_args_map = ["--fastqc_args": "'-t ${task.cpus}'"]

// After:
def preset_args_map = [:]
```

This means:
- `params.extra_trimgalore_args` still works for user-supplied args
- No FastQC runs inside the TRIMGALORE process
- Trimmed-read FastQC is deferred to a separate step (future enhancement)

### 3. Update nf-test snapshots

After running tests, update snapshots with `--update-snapshot` flag.
Snapshot changes are expected for:
- Report files (new v2.0 format with Cutadapt-compat section)
- Version strings (0.6.10 → 2.0.0)
- Possibly missing html/zip outputs if snapshots include them

---

## Prerequisites

### Branch setup

User has write access to nf-core/rnaseq (origin).

1. Create branch: `git checkout -b test/trimgalore-v2`
2. After changes: `git push -u origin test/trimgalore-v2`
3. Draft PR against `nf-core/rnaseq:dev`

---

## Test Strategy

### Phase 1: Module-level nf-tests

```bash
cd ~/Github/rnaseq
nf-test test modules/nf-core/trimgalore/tests/main.nf.test --update-snapshot
```

Validates: container pulls, trim_galore runs, output file patterns, version extraction.

### Phase 2: Subworkflow-level nf-tests

```bash
nf-test test subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/tests/main.nf.test --update-snapshot
```

Validates: full FASTQC → UMI → TRIMGALORE pipeline, report parsing, read count filtering.

### Phase 3: Full pipeline comparison

```bash
# Run with new container
nextflow run ~/Github/rnaseq -profile test,docker --outdir results_v2

# Compare key outputs against dev branch reference
# (alignment rates, gene counts, MultiQC report)
```

---

## Implementation Order

1. Create `test/trimgalore-v2` branch from `dev`
3. Swap container in module main.nf
4. Remove `--fastqc_args` from subworkflow config
5. Run module nf-tests → update snapshots
6. Run subworkflow nf-tests → update snapshots
7. Commit all changes
8. Create draft PR against `nf-core/rnaseq:dev`
9. Run full pipeline test and compare outputs
