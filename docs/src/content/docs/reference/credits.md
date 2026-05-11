---
title: Credits & license
description: Project history, maintenance, and license.
---

## History

Trim Galore was developed at [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/) by [Felix Krueger](https://github.com/FelixKrueger) (now part of [Altos Labs](https://altoslabs.com/)). For roughly a decade it shipped as a Perl wrapper around [Cutadapt](https://cutadapt.readthedocs.io/) and was the de-facto standard for adapter trimming in many bisulfite sequencing pipelines.

In 2026, Trim Galore was rewritten in Rust as a single static binary (v2.x). It is a drop-in replacement for v0.6.x with extra capabilities: poly-G auto-detection, a generic poly-A trimmer, per-pair adapter detection, cleaner multi-adapter syntax, and a structured JSON report. See [Migrating from v0.6.x](/reference/migration/) for the changes.

Current development is at [github.com/FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore).

## Contributors

Trim Galore is maintained by [Felix Krueger](https://github.com/FelixKrueger). The Rust rewrite has benefitted from contributions and design input from a number of people; the [GitHub contributors page](https://github.com/FelixKrueger/TrimGalore/graphs/contributors) tracks the full list. Notable recent contributions:

- **[Phil Ewels](https://github.com/ewels)** — bundled FastQC integration via [`fastqc-rust`](https://crates.io/crates/fastqc-rust), docs site infrastructure (Astro/Starlight), GitHub Actions hardening, hero animation polish.
- **[Dongze He](https://github.com/an-altosian) (`@an-altosian`)** — Phase-1B Perl-parity hunt and Buckberry-scale performance audit. Reported and prototyped the lowercase clip-flag fix (#242), \--max\_n fraction-mode UX, the gzip-extension and \--retain\_unpaired Perl-parity fixes (#245), the \--clock/\--implicon convenience widening (#245-D), the test-coverage gap inventory (#246), CI hardening recommendations (#247), and the Buckberry-scale performance audit (#248) including the gzip-level and single-buffered-write wins landed in v2.1.0-beta.6.

If you've contributed and would like a line here, please open a PR — the list is intentionally curated rather than auto-generated.

## Citing Trim Galore

If you use Trim Galore in published work, cite the project repository:

> Krueger, F. *Trim Galore!* Consistent adapter and quality trimming for FASTQ files. <https://github.com/FelixKrueger/TrimGalore>

Trim Galore has a Zenodo DOI for archival citation. Pick the version-specific DOI from the [release page](https://github.com/FelixKrueger/TrimGalore/releases) for the version you used.

## License

Trim Galore is released under the **GNU General Public License v3.0** (GPL-3.0). The full text is in the [LICENSE file](https://github.com/FelixKrueger/TrimGalore/blob/master/LICENSE).

## Reporting issues

- Bugs and feature requests: [GitHub Issues](https://github.com/FelixKrueger/TrimGalore/issues).
- Beta feedback: open an issue with the `beta-feedback` label.
- Security: contact the maintainer privately rather than filing a public issue.

## Related tools

- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Quality reports for raw and trimmed FASTQ.
- [MultiQC](https://multiqc.info). Aggregates Trim Galore reports (text and JSON) across samples.
- [Bismark](https://github.com/FelixKrueger/Bismark). Bisulfite read aligner. The standard downstream step for RRBS workflows.
- [UmiBam](https://github.com/FelixKrueger/Umi-Grinder). UMI-aware deduplication for Bismark BAM files. Used with `--clock` and `--implicon`.
- [Cutadapt](https://cutadapt.readthedocs.io/). The adapter-trimming engine that the Perl Trim Galore wrapped. v2 has its own built-in adapter trimmer.
