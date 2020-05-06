# Trim Galore Changelog

### Version 0.6.5dev

- Changed the way in which we test for the version of Cutadapt, more here: https://github.com/FelixKrueger/TrimGalore/issues/85

- Allowed specifying of multiple adapters for special cases. Works either via the command line, e.g.: `-a  " AGCTCCCG -a TTTCATTATAT -a TTTATTCGGATTTAT"` or via a FastA file, like so: `-a "file:multiple_adapters.fa"`  More info here: https://github.com/FelixKrueger/TrimGalore/issues/86.

- Added new special trimming mode for UMIs for the IMPLICON method ([`--implicon`](https://github.com/FelixKrueger/TrimGalore/issues/90)). In this mode, an 8bp UMI (unique molecular identifier) sequence is transferred from the start of Read 2 to the readID of both reads to allow UMI-aware deduplication (e.g. with `deduplicate_bismark --barcode` or [UmiBam](https://github.com/FelixKrueger/Umi-Grinder). Following this, Trim Galore will exit.

### Version 0.6.5 (Release on 19 Nov 2019)

- Added checks for whitespace(s) within input filenames, or a potential output folder name (supplied with `-o`). `[FATAL ERROR]` messages will advise users to use `_` instead.

- In a `--paired --basename BASE` scenario, the output files will now be called `BASE_val_1.fq.gz BASE_val_2.fq.gz` as described in the documentation (we previously also added `_R1` and `_R2`). This had to be addressed twice (0f631e5f979281fd4f18faef39818399a068a4b3 and 9ad019635a8a7f1aebb56f309889a7841a0ae42e) as the first approach was generating the Read 1 twice.

- removed a superflous warning statement for directional RRBS mode


### Version 0.6.4 (Release on 01 Aug 2019)

- Changed the adapter auto-detection procedure so that inconclusive detection always defaults to `--illumina`, unless none of the top 2, equal contaminants was 'Illumina', in which case it now defaults to `--nextera`. A warning message about this is now printed to the screen as well as to the trimming report.

- In addition to that, added the option `--consider_already_trimmed INT`. If no specific adapter exceeds this limit during the adapter auto-detection, the file is considered 'already adapter-trimmed' and will not be adapter trimmed again. Quality trimming is carried out as usual (technically, the adapter sequence is set to `-a X`). This option was added so that pipelines that are being fed either already-trimmed or untrimmed data will do the right thing in both cases.

- Changed the trimming mode for paired-end `--rrbs` in conjunction with `--non_directional`: previously, Read 2 was only trimmed for `CGA` or `CAA` at the 5' end, but not trimmed for read-through contamination at the 3' end if no 5' contamination had been removed. This problem had been introduced in v0.4.3, but since non-directional RRBS is not very common it had not been spotted so far. 

- File names for single-end trimming are now changed correctly when both `--output_dir` and `--basename` were specified together (was working correctly for PE mode already) 

### Version 0.6.3 (Release on 27 06 2019)

- Also added the number of PolyA trimmed bases to the start of the read in the format `trimmed_bases:A:`

So an example trimmed read would look like this:
```
 @READ-ID:1:1102:22039:36996 1:N:0:CCTAATCC
GCCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAATAAAAACTTTATAAACACCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

@32:A:READ-ID:1:1102:22039:36996_1:N:0:CCTAATCC_PolyA:32
GCCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAATAAAAACTTTATAAACACC

```


### Version 0.6.2 (Release on 08 05 2019)

- Changed the version checking mechanism so that Trim Galore even works in single-core mode if the version of Cutadapt was 7 years old...

- Fixed setting `-j $cores` for Cutadapt versions 2.X or above.


### Version 0.6.1 (Release on 20 Mar 2019)

* Added a check for very old versions Cutadapt, so that single-core trimming still works with Cutadapt versions prior to 1.15. 

- Fixed the way single-core trimming was dealt with in paired-end mode (which was introduced by the above 'fix')

- the option `--basename preferred_name` should now correclty work when specified in conjunction with `--output_dir`

### Version 0.6.0 (Release on 1 Mar 2019)

* Added option `--hardtrim3 INT,` which allows you to hard-clip sequences from their 5' end. This option processes one or more files (plain FastQ or gzip compressed files) and produces hard-trimmed FastQ files ending in .{INT}bp_3prime.fq(.gz). We found this quite useful in a number of scenarios where we wanted to removed biased residues from the start of sequences. Here is an example :

```
before:         CCTAAGGAAACAAGTACACTCCACACATGCATAAAGGAAATCAAATGTTATTTTTAAGAAAATGGAAAAT
--hardtrim3 20:                                                   TTTTTAAGAAAATGGAAAAT
```

* Added new option `--basename <PREFERRED_NAME>` to use `PREFERRED_NAME` as the basename for output files, instead of deriving the filenames from the input files. Single-end data would be called `PREFERRED_NAME_trimmed.fq(.gz)`, or `PREFERRED_NAME_val_1.fq(.gz)` and `PREFERRED_NAME_val_2.fq(.gz)` for paired-end data. `--basename` only works when 1 file (single-end) or 2 files (paired-end) are specified, but not for longer lists (see #43).

* Added option `--2colour/--nextseq INT` whereby INT selects the quality cutoff that is normally set with `-q`, only that qualities of G bases are ignored. `-q` and `--2colour/--nextseq INT` are mutually exclusive (see [#41](https://github.com/FelixKrueger/TrimGalore/issues/41))

* Added check to see if Read 1 and Read 2 files were given as the very same file.

* If an output directory which was specified with `-o output_directory` did not exist, it will be created for you.

* The option `--max_n INT` now also works in single-end RRBS mode.

* Added multi-threading support with the new option `-j/--cores INT`; many thanks to Frankie James for initiating this. Multi-threading support works effectively if Cutadapt is run with Python 3, and if parallel gzip (`pigz`) is installed:

<img title="Multi-threading benchmark" style="float:right;margin:20px 20 20 600px" id="Multi-threading support" src="Docs/Images/pigz_bench.png" >

For Cutadapt to work with multiple cores, it requires Python 3 as well as parallel gzip (pigz) installed on the system. The version of Python used is detected from the shebang line of the Cutadapt executable (either 'cutadapt', or a specified path). If Python 2 is detected, `--cores` is set to 1 and multi-core processing will be disabled. If `pigz` cannot be detected on your system, Trim Galore reverts to using `gzip` compression. Please note however, that `gzip` compression will slow down multi-core processes so much that it is hardly worthwhile, please see: [here](https://github.com/FelixKrueger/TrimGalore/issues/16#issuecomment-458557103) for more info).

Actual core usage: It should be mentioned that the actual number of cores used is a little convoluted. Assuming that Python 3 is used and `pigz` is installed, `--cores 2` would use:

- 2 cores to read the input (probably not at a high usage though)
- 2 cores to write to the output (at moderately high usage)
- 2 cores for Cutadapt itself
- 2 additional cores for Cutadapt (not sure what they are used for)
- 1 core for Trim Galore itself

So this can be up to 9 cores, even though most of them won't be used at 100% for most of the time. Paired-end processing uses twice as many cores for the validation (= writing out) step as Trim Galore reads and writes from and to two files at the same time, respectively.

`--cores 4` would then be: 4 (read) + 4 (write) + 4 (Cutadapt) + 2 (extra Cutadapt) +     1 (Trim Galore) = ~15 cores in total.

From the graph above it seems that `--cores 4` could be a sweet spot, anything above appear to have diminishing returns.


### 28-06-18: Version 0.5.0 

* Adapters can now be specified as single bases with a multiplier in squiggly brackets, e.g. -a "A{10}" to trim poly-A tails

* Added option `--hardtrim5 INT` to enable hard-clipping from the 5' end. This option processes one or more files (plain FastQ or gzip compressed files) and produce hard-trimmed FastQ files ending in `.{INT}bp.fq(.gz)`. 

* Added option `--clock` to trim reads in a specific way that is currently used for the Mouse Epigenetic Clock (see here: [Multi-tissue DNA methylation age predictor in mouse, Stubbs et al., Genome Biology, 2017 18:68](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1203-5)). Following the trimming, Trim Galore exits.

In it's current implementation, the dual-UMI RRBS reads come in the following format:

```
Read 1  5' UUUUUUUU CAGTA FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF TACTG UUUUUUUU 3'
Read 2  3' UUUUUUUU GTCAT FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF ATGAC UUUUUUUU 5'
```

Where UUUUUUUU is a random 8-mer unique molecular identifier (UMI), CAGTA is a constant region,
and FFFFFFF... is the actual RRBS-Fragment to be sequenced. The UMIs for Read 1 (R1) and
Read 2 (R2), as well as the fixed sequences (F1 or F2), are written into the read ID and
removed from the actual sequence. Here is an example:

```
R1: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 1:N:0: CGATGTTT
    ATCTAGTTCAGTACGGTGTTTTCGAATTAGAAAAATATGTATAGAGGAAATAGATATAAAGGCGTATTCGTTATTG
R2: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 3:N:0: CGATGTTT
    CAATTTTGCAGTACAAAAATAATACCTCCTCTATTTATCCAAAATCACAAAAAACCACCCACTTAACTTTCCCTAA

R1: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 1:N:0: CGATGTTT:R1:ATCTAGTT:R2:CAATTTTG:F1:CAGT:F2:CAGT
                 CGGTGTTTTCGAATTAGAAAAATATGTATAGAGGAAATAGATATAAAGGCGTATTCGTTATTG
R2: @HWI-D00436:407:CCAETANXX:1:1101:4105:1905 3:N:0: CGATGTTT:R1:ATCTAGTT:R2:CAATTTTG:F1:CAGT:F2:CAGT
                 CAAAAATAATACCTCCTCTATTTATCCAAAATCACAAAAAACCACCCACTTAACTTTCCCTAA
```
Following clock trimming, the resulting files (.clock_UMI.R1.fq(.gz) and .clock_UMI.R2.fq(.gz))
should be adapter- and quality trimmed with Trim Galore as usual. In addition, reads need to be trimmed
by 15bp from their 3' end to get rid of potential UMI and fixed sequences. The command is:

`trim_galore --paired --three_prime_clip_R1 15 --three_prime_clip_R2 15 *.clock_UMI.R1.fq.gz *.clock_UMI.R2.fq.gz`

Following this, reads should be aligned with Bismark and deduplicated with UmiBam in `--dual_index` mode (see here: https://github.com/FelixKrueger/Umi-Grinder). UmiBam recognises the UMIs within this pattern: R1:(**ATCTAGTT**):R2:(**CAATTTTG**): as (UMI R1=**ATCTAGTT**) and (UMI R2=**CAATTTTG**).


### 13-11-17: Version 0.4.5 released

* Trim Galore now dies during the validation step when it encounters paired-end files that are not equal in length


### 28-03-17: Version 0.4.4 released

* Reinstated functionality of option `--rrbs` for single-end RRBS reads which had gone amiss in the previous release. What happened in detail was that RRBS trimming was de facto skipped if there was only a single file specified.

* Updated User Guide and Readme documents, added Installation instruction and Travis functionality - thanks Phil!


### 25-01-17: Version 0.4.3 released

* Changed the option `--rrbs` for paired-end libraries from removing 2 additional base pairs from the 3' end of both reads to trim 2 bp from the 3' end only for Read 1 and set `--clip_r2 2` for Read 2 instead. This is because Read 2 does not technically need 3' trimming since the end of Read 2 is not affected by the artificial methylation states introduced by the [end-repair] fill-in reaction. Instead, the first couple of positions of Read 2 suffer from the same fill-in problems as standard [paired-end libraries](https://sequencing.qcfail.com/articles/library-end-repair-reaction-introduces-methylation-biases-in-paired-end-pe-bisulfite-seq-applications/). Also see [this issue](https://github.com/FelixKrueger/TrimGalore/issues/3).

* Added a closing statement for the REPORT filehandle since it occasionally swallowed the last line...

* Setting `--length` now takes priority over the smallRNA adapter (which would set the length cutoff to 18 bp).


### 07-09-16: Version 0.4.2 released
* Replaced `zcat` with `gunzip -c` so that older versions of Mac OSX do not append a .Z to the end of the file and subsequently fail because the file is not present. Dah...
*	Added option `--max_n COUNT` to remove all reads (or read pairs) exceeding this limit of tolerated Ns. In a paired-end setting it is sufficient if one read exceeds this limit. Reads (or read pairs) are removed altogether and are not further trimmed or written to the unpaired output
* Enabled option `--trim-n` to remove Ns from both end of the reads. Does currently not work for RRBS-mode
* Added new option `--max_length` which reads that are longer than bp after trimming. This is only advised for smallRNA sequencing to remove non-small RNA sequences

### 12-11-15: Version 0.4.1 released: Essential update for smallRNA libraries!
*	Changed the Illumina small RNA sequence used for auto-detection to `TGGAATTCTCGG` (from formerly `ATGGAATTCTCG`). The reason for this is that smallRNA libraries have ssRNA adapters ligated to their -OH end, a signature of dicer cleavage, so there is no A-tailing involved. Thanks to Q. Gouil for bringing this to our attention
*	Changed the length cutoff for sequences to 16bp (down from 20bp) for smallRNA libraries before sequences get removed entirely. This is because some 20-23bp long smallRNAs species that had removed T, TG, or TGG etc. might just about pass the 20bp cutoff
*	Added a small description to the `--help` message for users of the NuGEN Ovation RRBS kit as to *NOT* use the `--rrbs` option (see `--help`)

### 06-05-15: Version 0.4.0 released
*	Unless instructed otherwise Trim Galore will now attempt to auto-detect the adapter which had been used for library construction (choosing from the Illumina universal, Nextera transposase and Illumina small RNA adapters). For this the first 1 million sequences of the first file specified are analysed. If no adapter can be detected within the first 1 million sequences Trim Galore defaults to --illumina. The auto-detection behaviour can be overruled by specifying an adapter sequence or using `--illumina`, `--nextera` or `--small_rna`
*	Added the new options `--illumina`, `--nextera` and `--small_rna` to use different default sequences for trimming (instead of `-a`): Illumina: `AGATCGGAAGAGC`; Small RNA: `TGGAATTCTCGG`; Nextera: `CTGTCTCTTATA`
*	Added a sanity check to the start of a Trim Galore run to see if the (first) FastQ file in question does contain information at all or appears to be in SOLiD colorspace format, and bails if either is true. Trim Galore does not support colorspace trimming, but users wishing to do this are kindly referred to using Cutadapt as a standalone program
*	Added a new option `--path_to_cutadapt /path/to/cudapt`. Unless this option is specified it is assumed that Cutadapt is in the PATH (equivalent to `--path_to_cutadapt cutadapt`). Also added a test to see if Cutadapt seems to be working before the actual trimming is launched
*	Fixed an open command for a certain type of RRBS processing (was open() instead of open3())

### 16-07-14: Version 0.3.7 released
*	Applied small change that makes paired-end mode work again (it was accidentally broken by changing @ARGV for @filenames when looping through the filenames...)

### 11-07-14: Version 0.3.6 released
*	Added the new options `--three_prime_clip_r1` and `--three_prime_clip_r2` to clip any number of bases from the 3' end after adapter/quality trimming has completed
* Added a check to see if Cutadapt exits fine. Else, Trim Galore will bail a well
*	The option `--stringency` needs to be spelled out now since using `-s` was ambiguous because of `--suppress_warn`

### late 2013: Version 0.3.5 released
*	Added the Trim Galore version number to the summary report

### 19-09-13: Version 0.3.4 released
*	Added single-end or paired-end mode to the summary report
*	In paired-end mode, the Read 1 summary report will no longer state that no sequence have been discarded due to trimming. This will be stated in the trimming report of Read 2 once the validation step has been completed

### 10-09-13: Version 0.3.3 released
*	Fixed a bug what was accidentally introduced which would add an additional empty line in single-end trimming mode

### 03-09-13: Version 0.3.2 released
*	Specifying `--clip_R1` or `--clip_R2` will no longer attempt to clip sequences that have been adapter- or quality-trimmed below the clipping threshold
*	Specifying an output directory with `--rrbs` mode should now correctly create temporary files

### 15-07-13: Version 0.3.1 released
*	The default length cutoff is now set at an earlier timepoint to avoid a clash in paired-end mode when `--retain_unpaired` and individual read lengths for read 1 and read 2 had been defined

### 15-07-13: Version 0.3.0 released
*	Added the options `--clip_R1` and `--clip_R2` to trim off a fixed amount of bases at from the 5' end of reads. This can be useful if the quality is unusually low at the start, or whenever there is an undesired bias at the start of reads. An example for this could be PBAT-Seq in general, or the start of read 2 for every bisulfite-Seq paired-end library where end repair procedure introduces unmethylated cytosines. For more information on this see the M-bias section of the Bismark User Guide.

### 10-04-13: Version 0.2.8 released
* Trim Galore will now compress output files with GZIP on the fly instead of compressing the trimmed file once trimming has completed. In the interest of time temporary files are not being compressed
*	Added a small sanity check to exit if no files were supplied for trimming. Thanks to P. for 'bringing this to my attention'

### 01-03-13: Version 0.2.7 released
*	Added a new option `--dont_gzip` that will force the output files not to be gzip compressed. This overrides both the `--gzip` option or a .gz line ending of the input file(s)

### 07-02-13: Version 0.2.6 released
*	Fixes some bugs which would not gzip or run FastQC correctly when the option `-o` had been specified
*	When `--fastqc` is specified in paired-end mode the intermediate files '_trimmed.fq' are no longer analysed (only files '_val_1' and '_val_2')

### 19-10-12: Version 0.2.5 released
*	Added option `-o/--output_directory` to redirect all output (including temporary files) to another folder (required for implementation into Galaxy)
*	Added option `--no_report_file` to suppress a report file
*	Added option `--suppress_warn` to suppress any output to STDOUT or STDERR

### 02-10-12: Version 0.2.4 released
*	Removed the shorthand `-l` from the description as it might conflict with the paired-end options `-r1/--length1` or `-r2/--length2`. Please use `--length` instead
*	Changed the reporting to show the true Phred score quality cutoff
*	Corrected a typo in stringency...

### 31-07-12: Version 0.2.3 released
* Added an option `-e ERROR RATE` that allows one to specify the maximum error rate for trimming manually (the default is 0.1)

### 09-05-12: Version 0.2.2 released
*	Added an option `-a2/--adapter2` so that one can specify individual adapter sequences for the two reads of paired-end files; hereby the sequence provided as `-a/--adapter` is used to trim read 1, and the sequence provided as `-a2/--adapter2` is used to trim read 2. This option requires `--paired` to be specified as well

### 20-04-12: Version 0.2.1 released
*	Trim Galore! now has an option `--paired` which has the same functionality as the validate_paired_ends script we offered previously. This option discards read pairs if one (or both) reads of a read pair became shorter than a given length cutoff
*	Reads of a read-pair that are longer than a given threshold but for which the partner read has become too short can optionally be written out to single-end files. This ensures that the information of a read pair is not lost entirely if only one read is of good quality
*	Paired-end reads may be truncated by a further 1 bp from their 3' end to avoid problems with invalid alignments with Bowtie 1 (which regards alignments that contain each other as invalid...)
*	The output may be gzip compressed (this happens automatically if the input files were gzipped (i.e. end in .gz))
*	The documentation was extended substantially. We also added some recommendations for RRBS libraries for MseI digested material (recognition motif TTAA)

### 21-03-12: Version 0.1.4 released
* Phred33 (Sanger) encoding is now the default quality scheme
*	Fixed a bug for Phred64 encoding that would occur if several files were specified at once

### 14-03-12: Version 0.1.3 released
*	Initial stand-alone release; all basic functions working
*	Added the option `--fastqc_args` to pass extra options to FastQC

