# Trim Galore Changelog

## Changes since Trim Galore Release v0.4.2

### Trim Galore

Changed the option `--rrbs` for paired-end libraries from removing 2 additional base pairs from the 3' end of both reads to trim 2 bp from the 3' end only for Read 1 and set `--clip_r2 2` for Read 2 instead. This is because Read 2 does not technically need 3' trimming since the end of Read 2 is not affected by the artificial methylation states introduced by the [end-repair] fill-in reaction. Instead, the first couple of positions of Read 2 suffer from the same fill-in problems as standard [paired-end libraries](https://sequencing.qcfail.com/articles/library-end-repair-reaction-introduces-methylation-biases-in-paired-end-pe-bisulfite-seq-applications/). Also see [this issue](https://github.com/FelixKrueger/TrimGalore/issues/3).

Added a closing statement for the REPORT filehandle since it occasionally swallowed the last line...

Setting `--length` now takes priority over the smallRNA adapter (which would set the length cutoff to 18 bp).
