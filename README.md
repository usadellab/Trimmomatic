# Trimmomatic
A Java-based processing and trimming tool for Illumina NGS sequencing data, developed by the [Usadel lab](http://www.usadellab.org/), capable of single-end and paired-end read handling.
# Note 
While the software is licensed under the GPL, the adapter sequences are *not* included in the GPL part, but owned by and used with permission of Illumina. Oligonucleotide sequences © 2023 Illumina, Inc. All rights reserved.
# Quick start
## Installation
The easiest option is to download a binary release zip, and unpack it somewhere convenient. You'll need to modify the example command lines below to reference the trimmomatic JAR file and the location of the adapter fasta files. 

## Trimmomatic is a de.NBI & ELIXIR Service
<a href="https://denbi.de"><img src="https://www.denbi.de/templates/nbimaster/img/denbi-logo-color.svg" width="30%"></a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href="https://elixir-europe.org"><img src="https://raw.githubusercontent.com/elixir-europe/rdmkit/master/assets/img/elixir_logo_inverted.svg" width="15%"></a><br>
This software is provided as a service by the German Network for Bioinformatics Infrastructure (de.NBI). As the German node for ELIXIR, de.NBI contributes services like this one to a pan-European infrastructure for life science data. Together, de.NBI and ELIXIR offer a coordinated portfolio of resources — including databases, software, training, and cloud computing — to make it easier for scientists in academia and industry to find and share data, exchange expertise, and establish best practices.


## Build from Source

To build Trimmomatic from source, you will need a Java Development Kit (JDK 8 or higher) and Apache Maven.

Clone the repository, change into the top-level directory, and build using the following Maven command:

```bash
mvn clean package
```

## Paired End:

With most new data sets you can use gentle quality trimming and adapter clipping.

You often don't need leading and traling clipping. Also in general setting the `keepBothReads` to `True` can be useful when working with paired end data, you will keep even redunfant information but this likely makes your pipelines more manageable. Note the additional `:2` in front of the `True` (for `keepBothReads`) - this is the minimum adapter length in palindrome mode, you can even set this to 1. (Default is a very conservative 8)

If you have questions please don't hesitate to contact us, this is not necessarily one size fits all. (e.g. RNAseq expression analysis vs DNA assembly).

```
java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
```

for reference only (less sensitive for adapters)

```java -jar trimmomatic-0.35.jar PE -phred33 \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

This will perform the following:

* Remove adapters (`ILLUMINACLIP:TruSeq3-PE.fa:2:30:10`)
* Remove leading low quality or N bases (below quality 3) (`LEADING:3`)
* Remove trailing low quality or N bases (below quality 3) (`TRAILING:3`)
* Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (`SLIDINGWINDOW:4:15`)
* Drop reads below 36 bases long (`MINLEN:36`)

## Single End:
To perform the same steps using a single-ended adapter file, run:
```
java -jar trimmomatic-0.35.jar SE -phred33 \
input.fq.gz \
output.fq.gz \
ILLUMINACLIP:TruSeq3-SE:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
``` 
# Description of Trimming Steps

Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data. The selection of trimming steps and their associated parameters are supplied on the command line.

The current trimming steps are:

* `ILLUMINACLIP`: Cut adapter and other illumina-specific sequences from the read.
* `SLIDINGWINDOW`: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
* `LEADING`: Cut bases off the start of a read, if below a threshold quality
* `TRAILING`: Cut bases off the end of a read, if below a threshold quality
* `CROP`: Cut the read to a specified length
* `HEADCROP`: Cut the specified number of bases from the start of the read
* `MINLEN`: Drop the read if it is below a specified length
* `TOPHRED33`: Convert quality scores to Phred-33
* `TOPHRED64`: Convert quality scores to Phred-64

It works with FASTQ (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used), either uncompressed or gzipp'ed FASTQ. Use of gzip format is determined based on the `.gz` extension.

For single-ended data, one input and one output file are specified, plus the processing steps. For paired-end data, two input files are specified, and 4 output files, 2 for the 'paired' output where both reads survived the processing, and 2 for corresponding 'unpaired' output where a read survived, but the partner read did not.
 
# Running Trimmomatic

Since version 0.27, trimmomatic can be executed using -jar. The 'old' method, using the explicit class, continues to work.
## Paired End Mode:

```
java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] \
[-trimlog <logFile>] [-summary <summaryFile>] [-basein <templateInputFile>] [-baseout <templateOutputFile>] \
[-validatePairs] [-compressLevel <level>] [-compressStream | -compressBlock] [-quiet] [-version] \
<input 1> <input 2> \
<paired output 1> <unpaired output 1> \
<paired output 2> <unpaired output 2> \
<step 1> # Additional steps added as needed
```

or

```
java -classpath <path to trimmomatic jar> org.usadellab.trimmomatic.TrimmomaticPE [-threads <threads>] [-phred33 | -phred64] \
[-trimlog <logFile>] [-summary <summaryFile>] [-basein <templateInputFile>] [-baseout <templateOutputFile>] \
[-validatePairs] [-compressLevel <level>] [-compressStream | -compressBlock] [-quiet] [-version] \
<input 1> <input 2> \
<paired output 1> <unpaired output 1> \
<paired output 2> <unpaired output 2> \
<step 1> # Additional steps added as needed
```
## Single End Mode:

```
java -jar <path to trimmomatic jar> SE [-threads <threads>] [-phred33 | -phred64] \
[-trimlog <logFile>] [-summary <summaryFile>] [-compressLevel <level>] [-compressStream | -compressBlock] [-quiet] [-version] \
<input> <output> \
<step 1> # Additional steps added as needed
```

or

```
java -classpath <path to trimmomatic jar> org.usadellab.trimmomatic.TrimmomaticSE [-threads <threads>] [-phred33 | -phred64] \
[-trimlog <logFile>] [-summary <summaryFile>] [-compressLevel <level>] [-compressStream | -compressBlock] [-quiet] [-version] \
<input> <output> \
<step 1> # Additional steps added as needed
```

**Note**: Compression method is inferred from output file suffixes (`.gz`, `.bz2`, `.zip`).

#Optional parameters:

* `-threads <threads>`: specifies the number of CPU threads to use for multi-threading [default = 0]. If set to 0, Trimmomatic will try to calculate the number of available threads and run with that.
* `-phred33` | `-phred64`: specifies the quality score encoding. `phred33` is the standard for modern Illumina data.
* `-trimlog <logFile>`: writes a detailed log file.
* `-summary <summaryFile>`: writes a summary of trimming results to a file.
* `-basein <templateInputFile>`: path to one of the **paired-end** input files; its mate is auto-detected.
* `-baseout <templateOutputFile>`: template path used to generate the four **paired-end** output files (using the `_1P`, `_1U`, `_2P`, `_2U` suffixes).
* `-validatePairs`: performs an extra validation step on **paired-end** reads before trimming to ensure read pairs are consistent.
* `-compressLevel <level>`: sets the compression level for BZIP2/GZ output files (1=fastest, 9=best compression).
* `-compressStream` | `compressBlock`: specifies the compression mode. Block compression is the default.
* `-quiet`: suppresses progress output to the console.
* `-version`: prints the Trimmomatic version number to the console.

#Step options:

Multiple steps can be specified as required, by using additional arguments at the end.

Most steps take one or more settings, delimited by `:`.

* `ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLengthPalindrome>:<keepBothReads>`
    * `fastaWithAdaptersEtc`: specifies the path to a fasta file containing all the adapters, PCR sequences etc. The naming of the various sequences within this file determines how they are used. See below.
    * `seedMismatches`: specifies the maximum mismatch count which will still allow a full match to be performed
    * `palindromeClipThreshold`: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
    * `simpleClipThreshold`: specifies how accurate the match between any adapter etc. sequence must be against a read.
    * `minAdapterLengthPalindrome`: (optional int) specifies the minimum adapter length in palindrome mode [default = 8].
    * `keepBothReads`: (optional boolean) specifies if both reads should be kept in palindrome mode even when redundant information is found (small inserts)[default = False]. Note: `minAdapterLengthPalindrome` needs to be set manually to be able to activate `keepBothReads`.

* `LEADING:<quality>`
    * `quality`: specifies the minimum quality required to keep a base.

* `TRAILING:<quality>`
    * `quality`: specifies the minimum quality required to keep a base.

* `HEADCROP:<length>`
    * `length`: the number of bases to remove from the start of the read.

* `TAILCROP:<length>`
    * `length`: the number of bases to remove from the end of the read.

* `CROP:<length>`
    * `length`: the number of bases to keep, from the start of the read.

* `SLIDINGWINDOW:<windowSize>:<requiredQuality>`
    * `windowSize`: specifies the number of bases to average across
    * `requiredQuality`: specifies the average quality required.
    
* `MAXINFO:<targetLength>:<strictness>`
    * `targetLength`: the ideal length for a read. Scoring system favors reads trimmed to a length near this value.
    * `strictness`: specifies the trade-off between length and quality. A higher strictness value places more weight on base quality and less on the length score. Range between 0.0 and 1.0.

* `MINLEN:<length>`
    * `length`: specifies the minimum length of reads to be kept.

* `MAXLEN:<length>`
    * `length`: specifies the maximum length of reads to be kept.

* `AVGQUAL:<quality>`
    * `quality`: specifies the minimum Phred quality score of a read to be kept.
    
* `BASECOUNT:<bases>:<minCount>:<maxCount>`
    * `bases`: specifies the string of characters to be counted (e.g., `N`, `GC`).
    * `minCount`: (optional) the minimum number of times the `bases` must appear for the read to be kept.
    * `maxCount`: (optional) the maximum number of times the `bases` are allowed to appear.
    
* `TOPHRED33`: converts Phred-64 encoded records to Phred-33.

* `TOPHRED64`: converts Phred-33 encoded records to Phred-64.
    

# Step Order

Trimming occurs in the order which the steps are specified on the command line. It is recommended that adapter clipping, if required, is done as early as possible in most cases.
 
# The Adapter Fasta

Illumina adapter and other technical sequences are copyrighted by Illumina,but we have been granted permission to distribute them with Trimmomatic. Suggested adapter sequences are provided for TruSeq2 (as used in GAII machines) and TruSeq3 (as used by HiSeq and MiSeq machines), for both single-end and paired-end mode. These sequences have not been extensively tested, and depending on specific issues which may occur in library preparation, other sequences may work better for a given dataset.

To make a custom version of fasta, you must first understand how it will be used. Trimmomatic uses two strategies for adapter trimming: Palindrome and Simple

With 'simple' trimming, each adapter sequence is tested against the reads, and if a sufficiently accurate match is detected, the read is clipped appropriately.

'Palindrome' trimming is specifically designed for the case of 'reading through' a short fragment into the adapter sequence on the other end. In this approach, the appropriate adapter sequences are 'in silico ligated' onto the start of the reads, and the combined adapter+read sequences, forward and reverse are aligned. If they align in a manner which indicates 'read-through', the forward read is clipped and the reverse read dropped (since it contains no new data).

Naming of the sequences indicates how they should be used. For 'Palindrome' clipping, the sequence names should both start with 'Prefix', and end in '/1' for the forward adapter and '/2' for the reverse adapter. All other sequences are checked using 'simple' mode. Sequences with names ending in '/1' or '/2' will be checked only against the forward or reverse read. Sequences not ending in '/1' or '/2' will be checked against both the forward and reverse read. If you want to check for the reverse-complement of a specific sequence, you need to specifically include the reverse-complemented form of the sequence as well, with another name.

The thresholds used are a simplified log-likelihood approach. Each matching base adds just over 0.6, while each mismatch reduces the alignment score by Q/10. Therefore, a perfect match of a 12 base sequence will score just over 7, while 25 bases are needed to score 15. As such we recommend values between 7 - 15 for this parameter. For palindromic matches, a longer alignment is possible - therefore this threshold can be higher, in the range of 30. The 'seed mismatch' parameter is used to make alignments more efficient, specifying the maximum base mismatch count in the 'seed' (16 bases). Typical values here are 1 or 2.
 
