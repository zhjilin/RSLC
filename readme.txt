Scripts to perform the UMI counting and SSMD score calculation

1.	Installation
	
1.1 Download the file:
1.2 unzip 
1.3 Required dependencies:

(i)	perl packages: IO::Uncompress::Gunzip
(ii)	R: reshape2, plyr

2.	Counting Random Sequence Labels (RSLs)/UMIs
  
  2.1 GuideUMI_p0.1.pl

This perl script can be used to count RSLs/UMIs from only one fastq/fastq.gz file at a time.

2.1.1 Usage:

perl GuideUMI_p0.1.pl --library lib.csv Sample1_reads.fastq.gz

This will give 2 tables named as Sample1.count Sample1.UMI

2.1.2 Help: perl GuideUMI_p0.1.pl --help 

2.1.3 Instructions for file formats

Library file format: 

GuideID,GuideSequence,TargetGene

		AATF_03,GGACCCTGAAGCGGACCCCG,AATF
		AATF_04,GATGAAGGGGAAGATGGGGA,AATF
		AATF_05,CTTCAGATGAGCATTAGCAG,AATF
	
Note:  Any other form of annotation is not allowed.


Input file format: 
Both fastq and fastq.gz files are allowed.
Note: File name should have only one underscore in it. Whatever preceeds underscore in the input file name will be used as outputfile name. Data deposited under PRJEB18436.

2.2 BatchRun-p0.1.pl

This is a wrapper of the above script that will execute UMI counting for more than 2 input files and merge individual count files.

ALERT: Don't count too many files at the same time (< cpu numbers), because it runs as many processes as your input fastq files. This can kill your computer; thus please think about what you are going to do first! More than that, each fastq file probably requires a lot of memory, it means the memory will be file numbers* memory usage!!!


2.2.1 Usage:

perl BatchRun-p0.1.pl --library lib.csv --step 12 --fastq Sample1_reads.fq,Samplel2_reads.fastq.gz,Sample3_read.fastq.gz


This script will yield a SampleX.count files and a SampleX.UMI files [X stands for 1,2 or 3] for each of the input files. Additionally, it will merge all the count files into one file named as summary_count.output.raw. Subsequent analysis relies on this merged table.

The --step argument can also be used to only merge count files if UMIs were individually counted on several files using the GuideUMI_p0.1.pl script described above (section 2.1).

2.2.2 Help: perl BatchRun-p0.1.pl –help


3. Data normalization and SSMD Calculation


3.1 normalization20170111.R

This is the core R script to calculate the SSMD score for just for one treatment vs control.

 
3.1.1 Usage: 

Rscript --vanilla normalization20170111.R inputfile outputprefix count_threshold (if one wants to filter reads below certain counts)

3.1.2 Example: 

Rscript --vanilla normalization20170111.R summary_count.output.raw output_prefix 1 

3.2  SSMD_calculator.sh

This can be used in a time course experiment. For example when one wants to compare several time points of the screen with one control (Day4 vs Day10, Day15, Day20, Day30)

This shell script is wrapper to parse the big table in order to get pairwise tables (traits against control) to calculate SSMD score. Five arguments must be provided.

To see the help information: 
./SSMD_calculator.sh

3.2.1 Usage:

./SSMD_calculator.sh Rscript Inputfile output_prefix count_threshold (If one wants to filter reads below certain counts)

3.2.2 example

./SSMD_calculator.sh normalization20170111.R summary_count.output.raw ABC 1 


4.UMI Binning

4.1
countTruncatedRSLs.pl 

This is a perl script for binning CRISPR-Cas9 RSL guides based on the common RSL prefix. 

4.1.1 Usage: 
perl countTruncatedRSLs.pl <trunclen> <mincount> <input_countfile.csv>

The script reads the RSL guide counts the from the input file
<input_countfile.csv> and writes to standard output. It then bins
together all RSL guides that have the same first <trunclen> bases and
writes the sum counts of the truncated RSLs to standard output. Only
those RSL guides are considered in the binning that have been observed
at least <mincount> times in at least one of the samples.

4.1.2 Input/Output format:

The input file should be a comma separated text file containing the
following columns

RSL.guide,guide.set,Control,Treatment

The first column contains the guide name and its RSL separated by an
underscore '_'. The second column contains the guide set name. The
last two columns contain the RSL guide counts in two samples (control and treatment). 

Example:

RSL.guide,guide.set,Control,Treatment
AATF_01_AAAAGC,AATF_01,37,2

The output format is the same as input format.


