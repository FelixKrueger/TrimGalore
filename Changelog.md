# Trim Galore Changelog

## Changes since Trim Galore Release v0.4.2

### Trim Galore

* Changed the option `--rrbs` for paired-end libraries from removing 2 additional base pairs from the 3' end of both reads to trim 2 bp from the 3' end only for Read 1 and set `--clip_r2 2` for Read 2 instead. This is because Read 2 does not technically need 3' trimming since the end of Read 2 is not affected by the artificial methylation states introduced by the [end-repair] fill-in reaction. Instead, the first couple of positions of Read 2 suffer from the same fill-in problems as standard [paired-end libraries](https://sequencing.qcfail.com/articles/library-end-repair-reaction-introduces-methylation-biases-in-paired-end-pe-bisulfite-seq-applications/). Also see [this issue](https://github.com/FelixKrueger/TrimGalore/issues/3).

* Added a closing statement for the REPORT filehandle since it occasionally swallowed the last line...

* Setting `--length` now takes priority over the smallRNA adapter (which would set the length cutoff to 18 bp).


### 07-09-16: Version 0.4.2 released
* Replaced 'zcat' with 'gunzip -c' so that older versions of Mac OSX do not append a .Z to the end of the file and subsequently fail because the file is not present. Dah...
*	Added option '--max_n COUNT' to remove all reads (or read pairs) exceeding this limit of tolerated Ns. In a paired-end setting it is sufficient if one read exceeds this limit. Reads (or read pairs) are removed altogether and are not further trimmed or written to the unpaired output
* Enabled option '--trim-n' to remove Ns from both end of the reads. Does currently not work for RRBS-mode
* Added new option '--max_length ' which reads that are longer than bp after trimming. This is only advised for smallRNA sequencing to remove non-small RNA sequences

### 12-11-15: Version 0.4.1 released: Essential update for smallRNA libraries!
*	Changed the Illumina small RNA sequence used for auto-detection to 'TGGAATTCTCGG' (from formerly 'ATGGAATTCTCG'). The reason for this is that smallRNA libraries have ssRNA adapters ligated to their -OH end, a signature of dicer cleavage, so there is no A-tailing involved. Thanks to Q. Gouil for bringing this to our attention
*	Changed the length cutoff for sequences to 16bp (down from 20bp) for smallRNA libraries before sequences get removed entirely. This is because some 20-23bp long smallRNAs species that had removed T, TG, or TGG etc. might just about pass the 20bp cutoff
*	Added a small description to the --help message for users of the NuGEN Ovation RRBS kit as to *NOT* use the --rrbs option (see --help)

### 06-05-15: Version 0.4.0 released
*	Unless instructed otherwise Trim Galore will now attempt to auto-detect the adapter which had been used for library construction (choosing from the Illumina universal, Nextera transposase and Illumina small RNA adapters). For this the first 1 million sequences of the first file specified are analysed. If no adapter can be detected within the first 1 million sequences Trim Galore defaults to --illumina. The auto-detection behaviour can be overruled by specifying an adapter sequence or using --illumina, --nextera or --small_rna
*	Added the new options '--illumina', '--nextera' and '--small_rna' to use different default sequences for trimming (instead of -a): Illumina: AGATCGGAAGAGC; Small RNA: TGGAATTCTCGG; Nextera: CTGTCTCTTATA
*	Added a sanity check to the start of a Trim Galore run to see if the (first) FastQ file in question does contain information at all or appears to be in SOLiD colorspace format, and bails if either is true. Trim Galore does not support colorspace trimming, but users wishing to do this are kindly referred to using Cutadapt as a standalone program
*	Added a new option '--path_to_cutadapt /path/to/cudapt'. Unless this option is specified it is assumed that Cutadapt is in the PATH (equivalent to '--path_to_cutadapt cutadapt'). Also added a test to see if Cutadapt seems to be working before the actual trimming is launched
*	Fixed an open command for a certain type of RRBS processing (was open() instead of open3())

### 16-07-14: Version 0.3.7 released
*	Applied small change that makes paired-end mode work again (it was accidentally broken by changing @ARGV for @filenames when looping through the filenames...)

### 11-07-14: Version 0.3.6 released
*	Added the new options '--three_prime_clip_r1' and '--three_prime_clip_r2' to clip any number of bases from the 3' end after adapter/quality trimming has completed
* Added a check to see if Cutadapt exits fine. Else, Trim Galore will bail a well
*	The option '--stringency' needs to be spelled out now since using -s was ambiguous because of '--suppress_warn'

### late 2013: Version 0.3.5 released
*	Added the Trim Galore version number to the summary report

### 19-09-13: Version 0.3.4 released
*	Added single-end or paired-end mode to the summary report
*	In paired-end mode, the Read 1 summary report will no longer state that no sequence have been discarded due to trimming. This will be stated in the trimming report of Read 2 once the validation step has been completed

### 10-09-13: Version 0.3.3 released
*	Fixed a bug what was accidentally introduced which would add an additional empty line in single-end trimming mode

### 03-09-13: Version 0.3.2 released
*	Specifying --clip_R1 or --clip_R2 will no longer attempt to clip sequences that have been adapter- or quality-trimmed below the clipping threshold
*	Specifying an output directory with --rrbs mode should now correctly create temporary files

### 15-07-13: Version 0.3.1 released
*	The default length cutoff is now set at an earlier timepoint to avoid a clash in paired-end mode when '--retain_unpaired' and individual read lengths for read 1 and read 2 had been defined

### 15-07-13: Version 0.3.0 released
*	Added the options '--clip_R1' and '--clip_R2' to trim off a fixed amount of bases at from the 5' end of reads. This can be useful if the quality is unusually low at the start, or whenever there is an undesired bias at the start of reads. An example for this could be PBAT-Seq in general, or the start of read 2 for every bisulfite-Seq paired-end library where end repair procedure introduces unmethylated cytosines. For more information on this see the M-bias section of the Bismark User Guide

### 10-04-13: Version 0.2.8 released
* Trim Galore will now compress output files with GZIP on the fly instead of compressing the trimmed file once trimming has completed. In the interest of time temporary files are not being compressed
*	Added a small sanity check to exit if no files were supplied for trimming. Thanks to P. for 'bringing this to my attention'

### 01-03-13: Version 0.2.7 released
*	Added a new option '--dont_gzip' that will force the output files not to be gzip compressed. This overrides both the '--gzip' option or a .gz line ending of the input file(s)

### 07-02-13: Version 0.2.6 released
*	Fixes some bugs which would not gzip or run FastQC correctly when the option '-o' had been specified
*	When '--fastqc' is specified in paired-end mode the intermediate files '_trimmed.fq' are no longer analysed (only files '_val_1' and '_val_2')

### 19-10-12: Version 0.2.5 released
*	Added option '-o/--output_directory' to redirect all output (including temporary files) to another folder (required for implementation into Galaxy)
*	Added option '--no_report_file' to suppress a report file
*	Added option '--suppress_warn' to suppress any output to STDOUT or STDERR

### 02-10-12: Version 0.2.4 released
*	Removed the shorthand '-l' from the description as it might conflict with the paired-end options '-r1/--length1' or '-r2/--length2'. Please use '--length' instead
*	Changed the reporting to show the true Phred score quality cutoff
*	Corrected a typo in stringency...

### 31-07-12: Version 0.2.3 released
* Added an option `-e ERROR RATE` that allows one to specify the maximum error rate for trimming manually (the default is 0.1)

### 09-05-12: Version 0.2.2 released
*	Added an option `-a2/--adapter2` so that one can specify individual adapter sequences for the two reads of paired-end files; hereby the sequence provided as `-a/--adapter` is used to trim read 1, and the sequence provided as `-a2/--adapter2` is used to trim read 2. This option requires `--paired` to be specified as well

### 20-04-12: Version 0.2.1 released
*	Trim Galore! now has an option '--paired' which has the same functionality as the validate_paired_ends script we offered previously. This option discards read pairs if one (or both) reads of a read pair became shorter than a given length cutoff
*	Reads of a read-pair that are longer than a given threshold but for which the partner read has become too short can optionally be written out to single-end files. This ensures that the information of a read pair is not lost entirely if only one read is of good quality
*	Paired-end reads may be truncated by a further 1 bp from their 3' end to avoid problems with invalid alignments with Bowtie 1 (which regards alignments that contain each other as invalid...)
*	The output may be gzip compressed (this happens automatically if the input files were gzipped (i.e. end in .gz))
*	The documentation was extended substantially. We also added some recommendations for RRBS libraries for MseI digested material (recognition motif TTAA)

### 21-03-12: Version 0.1.4 released
* Phred33 (Sanger) encoding is now the default quality scheme
*	Fixed a bug for Phred64 encoding that would occur if several files were specified at once

### 14-03-12: Version 0.1.3 released
*	Initial stand-alone release; all basic functions working
*	Added the option 'fastqc_args' to pass extra options to FastQC

